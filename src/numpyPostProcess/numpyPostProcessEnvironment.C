/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  www.openfoam.com
    \\  /    A nd           | Copyright (C) 2026 Andre Weiner
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "numpyPostProcessEnvironment.H"

#include "addToRunTimeSelectionTable.H"
#include "fluidThermo.H"
#include "fvcFlux.H"
#include "fvc.H"
#include "gravityMeshObject.H"
#include "rhoThermo.H"
#include "singlePhaseTransportModel.H"
#include "surfaceFields.H"
#include "turbulentFluidThermoModel.H"
#include "turbulentTransportModel.H"
#include "volFields.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(numpyPostProcessEnvironment, 0);
    defineRunTimeSelectionTable(numpyPostProcessEnvironment, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

Foam::numpyPostProcessEnvironment::numpyPostProcessEnvironment
(
    fvMesh& mesh,
    const dictionary& dict,
    const wordHashSet& importedFields
)
:
    mesh_(mesh),
    importedFields_(importedFields),
    templateInstance_(dict.getOrDefault<word>("templateInstance", "0")),
    UPtr_()
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::volVectorField& Foam::numpyPostProcessEnvironment::velocity()
{
    if (volVectorField* U = mesh_.getObjectPtr<volVectorField>("U"))
    {
        return *U;
    }

    UPtr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "U",
                templateInstance_,
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        )
    );

    return UPtr_();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::numpyPostProcessEnvironment>
Foam::numpyPostProcessEnvironment::New
(
    fvMesh& mesh,
    const dictionary& dict,
    const wordHashSet& importedFields
)
{
    word environmentType(dict.getOrDefault<word>("type", "none"));

    if (environmentType == "simpleFoam" || environmentType == "pimpleFoam")
    {
        environmentType = "incompressible";
    }
    else if
    (
        environmentType == "rhoSimpleFoam"
     || environmentType == "rhoPimpleFoam"
    )
    {
        environmentType = "compressible";
    }
    else if (environmentType == "buoyantPimpleFoam")
    {
        environmentType = "buoyantCompressible";
    }

    auto* ctorPtr = dictionaryConstructorTable(environmentType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "environment",
            environmentType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    Info<< "Selecting NumPy post-processing environment "
        << environmentType << endl;

    return autoPtr<numpyPostProcessEnvironment>
    (
        ctorPtr(mesh, dict, importedFields)
    );
}


// * * * * * * * * * * * * * * Local Classes * * * * * * * * * * * * * * //

namespace Foam
{

class noneNumpyEnvironment
:
    public numpyPostProcessEnvironment
{
public:

    TypeName("none");

    noneNumpyEnvironment
    (
        fvMesh& mesh,
        const dictionary& dict,
        const wordHashSet& importedFields
    )
    :
        numpyPostProcessEnvironment(mesh, dict, importedFields)
    {}

    virtual void correct()
    {}
};


defineTypeNameAndDebug(noneNumpyEnvironment, 0);
addToRunTimeSelectionTable
(
    numpyPostProcessEnvironment,
    noneNumpyEnvironment,
    dictionary
);


class incompressibleNumpyEnvironment
:
    public numpyPostProcessEnvironment
{
    // Private Data

        autoPtr<surfaceScalarField> phiPtr_;
        autoPtr<singlePhaseTransportModel> transportPtr_;
        autoPtr<incompressible::turbulenceModel> turbulencePtr_;


public:

    TypeName("incompressible");

    incompressibleNumpyEnvironment
    (
        fvMesh& mesh,
        const dictionary& dict,
        const wordHashSet& importedFields
    )
    :
        numpyPostProcessEnvironment(mesh, dict, importedFields),
        phiPtr_(),
        transportPtr_(),
        turbulencePtr_()
    {
        volVectorField& U = velocity();

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::flux(U)
            )
        );
        transportPtr_.reset(new singlePhaseTransportModel(U, phiPtr_()));

        if
        (
            isFile
            (
                mesh_.time().globalPath()
               /mesh_.time().constant()/"turbulenceProperties"
            )
        )
        {
            turbulencePtr_ = incompressible::turbulenceModel::New
            (
                U,
                phiPtr_(),
                transportPtr_()
            );
        }
    }

    virtual void prepare()
    {
        const volVectorField& U = velocity();
        phiPtr_() = fvc::flux(U);
        transportPtr_->correct();
    }

    virtual void correct()
    {
        prepare();

        if (turbulencePtr_)
        {
            volScalarField& nut =
                mesh_.lookupObjectRef<volScalarField>("nut");

            for (auto& patchField : nut.boundaryFieldRef())
            {
                patchField.setUpdated(false);
            }

            turbulencePtr_->validate();
        }
    }
};


defineTypeNameAndDebug(incompressibleNumpyEnvironment, 0);
addToRunTimeSelectionTable
(
    numpyPostProcessEnvironment,
    incompressibleNumpyEnvironment,
    dictionary
);


class compressibleNumpyEnvironment
:
    public numpyPostProcessEnvironment
{
protected:

    // Protected Data

        autoPtr<fluidThermo> thermoPtr_;
        autoPtr<volScalarField> rhoPtr_;
        autoPtr<surfaceScalarField> phiPtr_;
        autoPtr<compressible::turbulenceModel> turbulencePtr_;


    // Protected Constructors

        compressibleNumpyEnvironment
        (
            fvMesh& mesh,
            const dictionary& dict,
            const wordHashSet& importedFields,
            const bool buoyant
        )
        :
            numpyPostProcessEnvironment(mesh, dict, importedFields),
            thermoPtr_(),
            rhoPtr_(),
            phiPtr_(),
            turbulencePtr_()
        {
            if (buoyant)
            {
                autoPtr<rhoThermo> thermo(rhoThermo::New(mesh_));
                thermoPtr_.reset(thermo.release());
            }
            else
            {
                thermoPtr_ = fluidThermo::New(mesh_);
            }

            if (!mesh_.foundObject<volScalarField>("rho"))
            {
                rhoPtr_.reset
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "rho",
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        thermoPtr_->rho()
                    )
                );
            }

            volVectorField& U = velocity();
            volScalarField& rho =
                mesh_.lookupObjectRef<volScalarField>("rho");

            phiPtr_.reset
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "phi",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    fvc::flux(rho*U)
                )
            );

            if
            (
                isFile
                (
                    mesh_.time().globalPath()
                   /mesh_.time().constant()/"turbulenceProperties"
                )
            )
            {
                turbulencePtr_ = compressible::turbulenceModel::New
                (
                    rho,
                    U,
                    phiPtr_(),
                    thermoPtr_()
                );
            }
        }


public:

    TypeName("compressible");

    compressibleNumpyEnvironment
    (
        fvMesh& mesh,
        const dictionary& dict,
        const wordHashSet& importedFields
    )
    :
        compressibleNumpyEnvironment
        (
            mesh,
            dict,
            importedFields,
            false
        )
    {}

    virtual void prepare()
    {
        fluidThermo& thermo = thermoPtr_();
        thermo.correct();

        volScalarField* rho = mesh_.getObjectPtr<volScalarField>("rho");

        if (!rho)
        {
            rhoPtr_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "rho",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    thermo.rho()
                )
            );
            rho = rhoPtr_.get();
        }
        else if (!importedFields_.contains("rho"))
        {
            *rho = thermo.rho();
        }

        const volVectorField& U = velocity();
        phiPtr_() = fvc::flux((*rho)*U);
    }

    virtual void correct()
    {
        prepare();

        if (turbulencePtr_)
        {
            volScalarField& nut =
                mesh_.lookupObjectRef<volScalarField>("nut");

            for (auto& patchField : nut.boundaryFieldRef())
            {
                patchField.setUpdated(false);
            }

            turbulencePtr_->validate();
        }
    }
};


defineTypeNameAndDebug(compressibleNumpyEnvironment, 0);
addToRunTimeSelectionTable
(
    numpyPostProcessEnvironment,
    compressibleNumpyEnvironment,
    dictionary
);


class buoyantCompressibleNumpyEnvironment
:
    public compressibleNumpyEnvironment
{
public:

    TypeName("buoyantCompressible");

    buoyantCompressibleNumpyEnvironment
    (
        fvMesh& mesh,
        const dictionary& dict,
        const wordHashSet& importedFields
    )
    :
        compressibleNumpyEnvironment(mesh, dict, importedFields, true)
    {
        (void)meshObjects::gravity::New(mesh_.time());
    }
};


defineTypeNameAndDebug(buoyantCompressibleNumpyEnvironment, 0);
addToRunTimeSelectionTable
(
    numpyPostProcessEnvironment,
    buoyantCompressibleNumpyEnvironment,
    dictionary
);

} // End namespace Foam


// ************************************************************************* //
