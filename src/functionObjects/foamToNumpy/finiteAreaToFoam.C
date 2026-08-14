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

#include "finiteAreaToFoam.H"

#include "HashPtrTable.H"
#include "IOobject.H"
#include "ListOps.H"
#include "Pstream.H"
#include "addToRunTimeSelectionTable.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "faMesh.H"
#include "numpyFileReader.H"
#include "numpyInputCatalog.H"

#include <memory>

// * * * * * * * * * * * * * * Local Classes * * * * * * * * * * * * * * //

class Foam::functionObjects::finiteAreaNumpyImporter
{
    const Time& time_;
    faMesh& mesh_;
    const wordList areaFields_;
    const wordList edgeFields_;
    const word templateInstance_;
    const bool writeFields_;
    const bool correctBoundaryConditions_;
    numpyDetail::numpyInputCatalog catalog_;
    scalar lastTime_;

    HashPtrTable<areaScalarField> areaScalarFields_;
    HashPtrTable<areaVectorField> areaVectorFields_;
    HashPtrTable<areaSphericalTensorField> areaSphericalTensorFields_;
    HashPtrTable<areaSymmTensorField> areaSymmTensorFields_;
    HashPtrTable<areaTensorField> areaTensorFields_;
    HashPtrTable<edgeScalarField> edgeScalarFields_;
    HashPtrTable<edgeVectorField> edgeVectorFields_;
    HashPtrTable<edgeSphericalTensorField> edgeSphericalTensorFields_;
    HashPtrTable<edgeSymmTensorField> edgeSymmTensorFields_;
    HashPtrTable<edgeTensorField> edgeTensorFields_;


    word fieldClass(const word& fieldName) const
    {
        if
        (
            const regIOobject* objectPtr =
                mesh_.thisDb().cfindObject<regIOobject>(fieldName)
        )
        {
            return objectPtr->type();
        }

        IOobject templateIO
        (
            fieldName,
            templateInstance_,
            mesh_.thisDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!templateIO.typeHeaderOk<regIOobject>(false))
        {
            FatalErrorInFunction
                << "Cannot read a finite-area template for " << fieldName
                << " from " << templateInstance_ << " on area "
                << mesh_.name() << exit(FatalError);
        }

        return templateIO.headerClassName();
    }


    template<class GeoField>
    GeoField& field
    (
        const word& fieldName,
        HashPtrTable<GeoField>& ownedFields
    )
    {
        if
        (
            GeoField* fieldPtr =
                mesh_.thisDb().getObjectPtr<GeoField>(fieldName)
        )
        {
            return *fieldPtr;
        }

        const GeoField templateField
        (
            IOobject
            (
                fieldName,
                templateInstance_,
                mesh_.thisDb(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_
        );

        GeoField* fieldPtr = new GeoField
        (
            IOobject
            (
                fieldName,
                time_.timeName(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            ),
            templateField
        );

        ownedFields.set(fieldName, fieldPtr);
        return *fieldPtr;
    }


    template<class GeoField>
    void loadField
    (
        const numpyDetail::numpyInputCatalog::snapshot& sample,
        const word& fieldName,
        HashPtrTable<GeoField>& ownedFields
    )
    {
        GeoField& target = field<GeoField>(fieldName, ownedFields);
        const fileName inputFile
        (
            catalog_.fieldPath(sample, fieldName, Pstream::myProcNo())
        );

        if (!isFile(inputFile))
        {
            FatalErrorInFunction
                << "Cannot find NumPy finite-area field file " << inputFile
                << exit(FatalError);
        }

        numpyDetail::numpyFileReader reader(inputFile);
        Field<typename GeoField::value_type> values;
        reader.readField(sample.index, target.size(), values);
        target.primitiveFieldRef() = values;
        if (correctBoundaryConditions_)
        {
            target.correctBoundaryConditions();
        }

        if (!target.nOldTimes())
        {
            (void)target.oldTime();
        }
    }


    void loadOne
    (
        const numpyDetail::numpyInputCatalog::snapshot& sample,
        const word& fieldName,
        const bool areaField
    )
    {
        const word className(fieldClass(fieldName));

        #define loadFieldType(FieldType, table)                               \
            if (className == FieldType::typeName)                            \
            {                                                                \
                loadField(sample, fieldName, table);                         \
                return;                                                      \
            }

        if (areaField)
        {
            loadFieldType(areaScalarField, areaScalarFields_);
            loadFieldType(areaVectorField, areaVectorFields_);
            loadFieldType
            (
                areaSphericalTensorField,
                areaSphericalTensorFields_
            );
            loadFieldType(areaSymmTensorField, areaSymmTensorFields_);
            loadFieldType(areaTensorField, areaTensorFields_);
        }
        else
        {
            loadFieldType(edgeScalarField, edgeScalarFields_);
            loadFieldType(edgeVectorField, edgeVectorFields_);
            loadFieldType
            (
                edgeSphericalTensorField,
                edgeSphericalTensorFields_
            );
            loadFieldType(edgeSymmTensorField, edgeSymmTensorFields_);
            loadFieldType(edgeTensorField, edgeTensorFields_);
        }

        #undef loadFieldType

        FatalErrorInFunction
            << "Unsupported " << (areaField ? "area" : "edge")
            << " field class " << className << " for " << fieldName
            << exit(FatalError);
    }


    template<class GeoField>
    bool writeField(const word& fieldName)
    {
        GeoField* fieldPtr =
            mesh_.thisDb().getObjectPtr<GeoField>(fieldName);

        if (!fieldPtr)
        {
            return false;
        }

        const word oldInstance(fieldPtr->instance());
        fieldPtr->instance() = time_.timeName();
        const bool written = fieldPtr->write();
        fieldPtr->instance() = oldInstance;
        return written;
    }


    void writeOne(const word& fieldName)
    {
        const word className(fieldClass(fieldName));
        bool written = false;

        #define writeFieldType(FieldType)                                    \
            if (!written && className == FieldType::typeName)                \
            {                                                                \
                written = writeField<FieldType>(fieldName);                  \
            }

        writeFieldType(areaScalarField);
        writeFieldType(areaVectorField);
        writeFieldType(areaSphericalTensorField);
        writeFieldType(areaSymmTensorField);
        writeFieldType(areaTensorField);
        writeFieldType(edgeScalarField);
        writeFieldType(edgeVectorField);
        writeFieldType(edgeSphericalTensorField);
        writeFieldType(edgeSymmTensorField);
        writeFieldType(edgeTensorField);

        #undef writeFieldType

        if (!written)
        {
            FatalErrorInFunction
                << "Cannot write imported finite-area field " << fieldName
                << exit(FatalError);
        }
    }


public:

    finiteAreaNumpyImporter
    (
        const Time& time,
        faMesh& mesh,
        const dictionary& dict,
        const fileName& areaInputPath
    )
    :
        time_(time),
        mesh_(mesh),
        areaFields_(dict.getOrDefault<wordList>("areaFields", wordList())),
        edgeFields_(dict.getOrDefault<wordList>("edgeFields", wordList())),
        templateInstance_
        (
            dict.getOrDefault<word>("templateInstance", "0")
        ),
        writeFields_(dict.getOrDefault("writeFields", false)),
        correctBoundaryConditions_
        (
            dict.getOrDefault("correctBoundaryConditions", false)
        ),
        catalog_
        (
            time.globalPath(),
            [&]()
            {
                dictionary areaDict(dict);
                areaDict.set("inputDir", areaInputPath);
                return areaDict;
            }()
        ),
        lastTime_(-VGREAT),
        areaScalarFields_(),
        areaVectorFields_(),
        areaSphericalTensorFields_(),
        areaSymmTensorFields_(),
        areaTensorFields_(),
        edgeScalarFields_(),
        edgeVectorFields_(),
        edgeSphericalTensorFields_(),
        edgeSymmTensorFields_(),
        edgeTensorFields_()
    {}


    bool load(const scalar timeValue)
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
            catalog_.find(timeValue);

        for (const word& fieldName : areaFields_)
        {
            loadOne(sample, fieldName, true);
        }
        for (const word& fieldName : edgeFields_)
        {
            loadOne(sample, fieldName, false);
        }

        lastTime_ = timeValue;
        Info<< "finiteAreaToFoam: imported area " << mesh_.name()
            << " at time " << Time::timeName(timeValue) << endl;
        return true;
    }


    void write()
    {
        if (!writeFields_)
        {
            return;
        }

        for (const word& fieldName : areaFields_)
        {
            writeOne(fieldName);
        }
        for (const word& fieldName : edgeFields_)
        {
            writeOne(fieldName);
        }
    }


    scalarList times() const
    {
        return catalog_.times();
    }
};


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(finiteAreaToFoam, 0);
    addToRunTimeSelectionTable(functionObject, finiteAreaToFoam, dictionary);
}
}


// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::functionObjects::finiteAreaToFoam::finiteAreaToFoam
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    areaNames_(),
    ownedMeshes_(),
    importers_()
{
    read(dict);
}


Foam::functionObjects::finiteAreaToFoam::finiteAreaToFoam
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, obr, dict),
    areaNames_(),
    ownedMeshes_(),
    importers_()
{
    read(dict);
}


Foam::functionObjects::finiteAreaToFoam::~finiteAreaToFoam()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void Foam::functionObjects::finiteAreaToFoam::initialise
(
    const dictionary& dict
)
{
    const objectRegistry* registry = faMesh::registry(mesh_);

    if (areaNames_.empty() && registry)
    {
        areaNames_ = registry->sortedNames<faMesh>();
    }

    if (areaNames_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "No finite-area mesh is registered and no 'area' or 'areas' "
            << "entry was specified"
            << exit(FatalIOError);
    }

    fileName basePath(dict.get<fileName>("inputDir"));
    if (!basePath.isAbsolute())
    {
        basePath = time_.globalPath()/basePath;
    }
    basePath.clean();

    ownedMeshes_.resize(areaNames_.size());
    importers_.resize(areaNames_.size());

    forAll(areaNames_, areai)
    {
        const word& areaName = areaNames_[areai];
        faMesh* areaMeshPtr = registry
          ? const_cast<faMesh*>(registry->cfindObject<faMesh>(areaName))
          : nullptr;

        if (!areaMeshPtr)
        {
            autoPtr<faMesh> meshPtr(faMesh::TryNew(areaName, mesh_));
            if (!meshPtr)
            {
                FatalIOErrorInFunction(dict)
                    << "Cannot find or read finite-area mesh " << areaName
                    << exit(FatalIOError);
            }
            ownedMeshes_.set(areai, meshPtr.release());
            areaMeshPtr = &ownedMeshes_[areai];
        }

        importers_.set
        (
            areai,
            new finiteAreaNumpyImporter
            (
                time_,
                *areaMeshPtr,
                dict,
                basePath/areaName
            )
        );
    }
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::finiteAreaToFoam::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (!importers_.empty())
    {
        return true;
    }

    dict.readIfPresent("areas", areaNames_);
    if (areaNames_.empty())
    {
        word areaName;
        if (dict.readIfPresent("area", areaName))
        {
            areaNames_.resize(1);
            areaNames_.front() = areaName;
        }
    }

    const wordList areaFields
    (
        dict.getOrDefault<wordList>("areaFields", wordList())
    );
    const wordList edgeFields
    (
        dict.getOrDefault<wordList>("edgeFields", wordList())
    );

    if (areaFields.empty() && edgeFields.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Specify at least one of 'areaFields' or 'edgeFields'"
            << exit(FatalIOError);
    }

    initialise(dict);
    return true;
}


bool Foam::functionObjects::finiteAreaToFoam::load(const scalar timeValue)
{
    for (finiteAreaNumpyImporter& importer : importers_)
    {
        importer.load(timeValue);
    }
    return true;
}


bool Foam::functionObjects::finiteAreaToFoam::execute()
{
    return load(time_.value());
}


bool Foam::functionObjects::finiteAreaToFoam::write()
{
    for (finiteAreaNumpyImporter& importer : importers_)
    {
        importer.write();
    }
    return true;
}


Foam::scalarList Foam::functionObjects::finiteAreaToFoam::times() const
{
    const scalarList values(importers_.front().times());

    for (label areai = 1; areai < importers_.size(); ++areai)
    {
        const scalarList areaTimes(importers_[areai].times());

        if (areaTimes.size() != values.size())
        {
            FatalErrorInFunction
                << "Finite-area inputs have different time sets"
                << exit(FatalError);
        }

        forAll(values, timei)
        {
            if
            (
                mag(areaTimes[timei] - values[timei])
              > SMALL*max(scalar(1), mag(values[timei]))
            )
            {
                FatalErrorInFunction
                    << "Finite-area inputs have different time sets"
                    << exit(FatalError);
            }
        }
    }

    return values;
}


// ************************************************************************* //
