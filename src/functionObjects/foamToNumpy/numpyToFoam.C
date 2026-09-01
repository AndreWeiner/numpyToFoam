/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  www.openfoam.com
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tanuj Ravi
    Copyright (C) 2026 Keysight Technologies
    Copyright (C) 2026 Andre Weiner
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

#include "numpyToFoam.H"

#include "IOobject.H"
#include "ListOps.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(numpyToFoam, 0);
    addToRunTimeSelectionTable(functionObject, numpyToFoam, dictionary);
}
}


// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::functionObjects::numpyToFoam::numpyToFoam
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldNames_(),
    templateInstance_("0"),
    writeFields_(false),
    correctBoundaryConditions_(false),
    catalog_(),
    lastTime_(-VGREAT),
    scalarFields_(),
    vectorFields_(),
    sphericalTensorFields_(),
    symmTensorFields_(),
    tensorFields_()
{
    read(dict);
}


Foam::functionObjects::numpyToFoam::~numpyToFoam()
{}


Foam::functionObjects::numpyToFoam::numpyToFoam
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, obr, dict),
    fieldNames_(),
    templateInstance_("0"),
    writeFields_(false),
    correctBoundaryConditions_(false),
    catalog_(),
    lastTime_(-VGREAT),
    scalarFields_(),
    vectorFields_(),
    sphericalTensorFields_(),
    symmTensorFields_(),
    tensorFields_()
{
    read(dict);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

Foam::word Foam::functionObjects::numpyToFoam::fieldClass
(
    const word& fieldName
) const
{
    if (const regIOobject* objectPtr = mesh_.cfindObject<regIOobject>(fieldName))
    {
        return objectPtr->type();
    }

    IOobject templateIO
    (
        fieldName,
        templateInstance_,
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (!templateIO.typeHeaderOk<regIOobject>(false))
    {
        FatalErrorInFunction
            << "Cannot read a field template for " << fieldName << " from "
            << templateInstance_ << nl
            << "Imported fields require a template for dimensions and "
            << "boundary conditions unless a compatible field is already "
            << "registered."
            << exit(FatalError);
    }

    return templateIO.headerClassName();
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::numpyToFoam::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    wordList newFieldNames(dict.get<wordList>("fields"));
    inplaceUniqueSort(newFieldNames);

    const word newTemplateInstance
    (
        dict.getOrDefault<word>("templateInstance", "0")
    );
    const bool newWriteFields = dict.getOrDefault("writeFields", false);
    const bool newCorrectBoundaryConditions =
        dict.getOrDefault("correctBoundaryConditions", false);

    if
    (
        catalog_
     &&
        (
            newFieldNames != fieldNames_
         || newTemplateInstance != templateInstance_
        )
    )
    {
        FatalIOErrorInFunction(dict)
            << "Input fields and templates cannot change after construction"
            << exit(FatalIOError);
    }

    fieldNames_ = std::move(newFieldNames);
    templateInstance_ = newTemplateInstance;
    writeFields_ = newWriteFields;
    correctBoundaryConditions_ = newCorrectBoundaryConditions;

    if (!catalog_)
    {
        catalog_.reset
        (
            new numpyDetail::numpyInputCatalog(time_.globalPath(), dict)
        );
    }

    return true;
}


bool Foam::functionObjects::numpyToFoam::load(const scalar timeValue)
{
    if (mag(timeValue - lastTime_) <= SMALL)
    {
        return true;
    }

    if (timeValue < lastTime_)
    {
        FatalErrorInFunction
            << "Import time moved backwards from " << lastTime_
            << " to " << timeValue
            << exit(FatalError);
    }

    const numpyDetail::numpyInputCatalog::snapshot& sample =
        catalog_->find(timeValue);

    for (const word& fieldName : fieldNames_)
    {
        const word className(fieldClass(fieldName));

        if (className == volScalarField::typeName)
        {
            loadField(sample, fieldName, scalarFields_);
        }
        else if (className == volVectorField::typeName)
        {
            loadField(sample, fieldName, vectorFields_);
        }
        else if (className == volSphericalTensorField::typeName)
        {
            loadField(sample, fieldName, sphericalTensorFields_);
        }
        else if (className == volSymmTensorField::typeName)
        {
            loadField(sample, fieldName, symmTensorFields_);
        }
        else if (className == volTensorField::typeName)
        {
            loadField(sample, fieldName, tensorFields_);
        }
        else
        {
            FatalErrorInFunction
                << "Unsupported class " << className
                << " for imported field " << fieldName
                << exit(FatalError);
        }
    }

    lastTime_ = timeValue;
    Log << type() << ' ' << name() << ": imported time "
        << Time::timeName(timeValue) << endl;
    return true;
}


bool Foam::functionObjects::numpyToFoam::execute()
{
    return load(time_.value());
}


void Foam::functionObjects::numpyToFoam::correctBoundaryConditions()
{
    for (const word& fieldName : fieldNames_)
    {
        const word className(fieldClass(fieldName));

        if (className == volScalarField::typeName)
        {
            correctFieldBoundaryConditions
            (
                mesh_.lookupObjectRef<volScalarField>(fieldName)
            );
        }
        else if (className == volVectorField::typeName)
        {
            correctFieldBoundaryConditions
            (
                mesh_.lookupObjectRef<volVectorField>(fieldName)
            );
        }
        else if (className == volSphericalTensorField::typeName)
        {
            correctFieldBoundaryConditions
            (
                mesh_.lookupObjectRef<volSphericalTensorField>(fieldName)
            );
        }
        else if (className == volSymmTensorField::typeName)
        {
            correctFieldBoundaryConditions
            (
                mesh_.lookupObjectRef<volSymmTensorField>(fieldName)
            );
        }
        else if (className == volTensorField::typeName)
        {
            correctFieldBoundaryConditions
            (
                mesh_.lookupObjectRef<volTensorField>(fieldName)
            );
        }
    }
}


bool Foam::functionObjects::numpyToFoam::write()
{
    if (!writeFields_)
    {
        return true;
    }

    for (const word& fieldName : fieldNames_)
    {
        const word className(fieldClass(fieldName));
        bool written = false;

        if (className == volScalarField::typeName)
        {
            written = writeField<volScalarField>(fieldName);
        }
        else if (className == volVectorField::typeName)
        {
            written = writeField<volVectorField>(fieldName);
        }
        else if (className == volSphericalTensorField::typeName)
        {
            written = writeField<volSphericalTensorField>(fieldName);
        }
        else if (className == volSymmTensorField::typeName)
        {
            written = writeField<volSymmTensorField>(fieldName);
        }
        else if (className == volTensorField::typeName)
        {
            written = writeField<volTensorField>(fieldName);
        }

        if (!written)
        {
            FatalErrorInFunction
                << "Cannot write imported field " << fieldName
                << exit(FatalError);
        }
    }

    return true;
}


Foam::scalarList Foam::functionObjects::numpyToFoam::times() const
{
    return catalog_->times();
}


// ************************************************************************* //
