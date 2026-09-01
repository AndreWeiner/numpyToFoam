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

#include "finiteAreaToNumpy.H"

#include "HashPtrTable.H"
#include "HashSet.H"
#include "IOobjectList.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "faMesh.H"
#include "mapPolyMesh.H"
#include "numpyFileWriter.H"

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <vector>

// * * * * * * * * * * * * * * Local Classes * * * * * * * * * * * * * * //

class Foam::functionObjects::finiteAreaNumpyWriter
{
    struct fieldInfo
    {
        word name;
        word className;
        label localSize;
        direction nComponents;
    };

    const word functionName_;
    const Time& time_;
    const faMesh& mesh_;
    const wordRes areaSelection_;
    const wordRes edgeSelection_;
    const numpyDetail::dataType dtype_;
    const label batchSize_;
    const bool writeTimes_;
    const bool writeAreaCentres_;
    const bool writeFaceAreas_;
    const bool writeEdgeCentres_;
    const fileName outputPath_;

    fileName segmentPath_;
    fileName batchPath_;
    bool segmentReady_;
    bool batchOpen_;
    bool geometryDirty_;
    label batchIndex_;
    label batchCount_;
    label meshRevision_;
    scalar firstBatchTime_;
    scalar lastBatchTime_;
    scalar lastOutputTime_;
    std::vector<fieldInfo> batchFields_;
    HashPtrTable<numpyDetail::numpyFileWriter> fieldWriters_;
    std::unique_ptr<numpyDetail::numpyFileWriter> timesWriter_;


    template<class GeoField>
    void appendFieldInfo
    (
        std::vector<fieldInfo>& info,
        HashSet<word>& selected,
        const wordRes& selection
    ) const
    {
        HashSet<word> candidates;

        for
        (
            const word& fieldName
          : mesh_.thisDb().sortedNames<GeoField>(selection)
        )
        {
            candidates.insert(fieldName);
        }

        if (functionObject::postProcess)
        {
            const IOobjectList objects(mesh_.thisDb(), time_.timeName());

            for
            (
                const word& fieldName
              : objects.names<GeoField>(selection, true)
            )
            {
                candidates.insert(fieldName);
            }
        }

        for (const word& fieldName : candidates.sortedToc())
        {
            if (!selected.insert(fieldName))
            {
                continue;
            }

            const GeoField* fieldPtr =
                mesh_.thisDb().cfindObject<GeoField>(fieldName);

            if (fieldPtr)
            {
                info.push_back
                ({
                    fieldName,
                    GeoField::typeName,
                    fieldPtr->size(),
                    pTraits<typename GeoField::value_type>::nComponents
                });
            }
            else
            {
                const GeoField field
                (
                    IOobject
                    (
                        fieldName,
                        time_.timeName(),
                        mesh_.thisDb(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                );

                info.push_back
                ({
                    fieldName,
                    GeoField::typeName,
                    field.size(),
                    pTraits<typename GeoField::value_type>::nComponents
                });
            }
        }
    }


    std::vector<fieldInfo> selectedFields() const
    {
        std::vector<fieldInfo> info;
        HashSet<word> selected;

        appendFieldInfo<areaScalarField>(info, selected, areaSelection_);
        appendFieldInfo<areaVectorField>(info, selected, areaSelection_);
        appendFieldInfo<areaSphericalTensorField>
        (
            info,
            selected,
            areaSelection_
        );
        appendFieldInfo<areaSymmTensorField>(info, selected, areaSelection_);
        appendFieldInfo<areaTensorField>(info, selected, areaSelection_);

        appendFieldInfo<edgeScalarField>(info, selected, edgeSelection_);
        appendFieldInfo<edgeVectorField>(info, selected, edgeSelection_);
        appendFieldInfo<edgeSphericalTensorField>
        (
            info,
            selected,
            edgeSelection_
        );
        appendFieldInfo<edgeSymmTensorField>(info, selected, edgeSelection_);
        appendFieldInfo<edgeTensorField>(info, selected, edgeSelection_);

        std::sort
        (
            info.begin(),
            info.end(),
            [](const fieldInfo& a, const fieldInfo& b)
            {
                return a.name < b.name;
            }
        );

        if (info.empty())
        {
            FatalErrorInFunction
                << "No supported finite-area fields match areaFields "
                << areaSelection_ << " or edgeFields " << edgeSelection_
                << " on area " << mesh_.name() << nl
                << "Available objects: " << mesh_.thisDb().sortedToc()
                << exit(FatalError);
        }

        return info;
    }


    template<class GeoField>
    bool appendField(const word& fieldName)
    {
        const GeoField* fieldPtr =
            mesh_.thisDb().cfindObject<GeoField>(fieldName);

        if (!fieldPtr)
        {
            if (!functionObject::postProcess)
            {
                return false;
            }

            IOobject io
            (
                fieldName,
                time_.timeName(),
                mesh_.thisDb(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            if (!io.typeHeaderOk<GeoField>(false))
            {
                return false;
            }

            const GeoField field(io, mesh_);
            fieldWriters_[fieldName]->appendField(field.primitiveField());
            return true;
        }

        fieldWriters_[fieldName]->appendField(fieldPtr->primitiveField());
        return true;
    }


    static word indexedName(const char* prefix, const label index)
    {
        std::ostringstream os;
        os << prefix << std::setw(6) << std::setfill('0') << index;
        return word(os.str());
    }


    void ensureSegment()
    {
        if (segmentReady_)
        {
            return;
        }

        fileName segmentName
        (
            functionObject::postProcess
          ? time_.timeName()
          : Time::timeName(time_.startTime().value())
        );

        if (Pstream::master())
        {
            const fileName initialName(segmentName);
            label suffix = 0;

            while (isDir(outputPath_/segmentName))
            {
                segmentName = initialName + '_' + Foam::name(++suffix);
            }
        }

        Pstream::broadcast(segmentName);
        segmentPath_ = outputPath_/segmentName;

        if (Pstream::master())
        {
            mkDir(segmentPath_);
            OFstream os(segmentPath_/"segmentInfo");
            IOobject::writeBanner(os);
            os  << "FoamFile" << nl
                << token::BEGIN_BLOCK << nl
                << "    version     2.0;" << nl
                << "    format      ascii;" << nl
                << "    class       dictionary;" << nl
                << "    object      segmentInfo;" << nl
                << token::END_BLOCK << nl;
            IOobject::writeDivider(os) << nl;

            os  << "meshType       finiteArea;" << nl
                << "area            " << mesh_.name() << token::END_STATEMENT
                << nl
                << "dataType        " << numpyDetail::dataTypeName(dtype_)
                << token::END_STATEMENT << nl
                << "order           fortran;" << nl
                << "batchSize       " << batchSize_ << token::END_STATEMENT
                << nl;
        }

        UPstream::barrier(UPstream::worldComm);
        segmentReady_ = true;
        Info<< "finiteAreaToNumpy " << functionName_ << ": writing area "
            << mesh_.name() << " to " << segmentPath_ << endl;
    }


    void writeState(const bool sealed) const
    {
        if (!Pstream::master() || !batchOpen_)
        {
            return;
        }

        const fileName stateFile(batchPath_/"state");
        const fileName temporaryFile(stateFile + ".tmp");
        {
            OFstream os(temporaryFile);

            IOobject::writeBanner(os);
            os  << "FoamFile" << nl
                << token::BEGIN_BLOCK << nl
                << "    version     2.0;" << nl
                << "    format      ascii;" << nl
                << "    class       dictionary;" << nl
                << "    object      state;" << nl
                << token::END_BLOCK << nl;
            IOobject::writeDivider(os) << nl;

            os  << "batch           " << batchIndex_
                << token::END_STATEMENT << nl
                << "meshRevision    " << meshRevision_
                << token::END_STATEMENT << nl
                << "count           " << batchCount_
                << token::END_STATEMENT << nl
                << "sealed          " << Switch(sealed)
                << token::END_STATEMENT << nl;

            if (batchCount_)
            {
                os  << "firstTime       " << firstBatchTime_
                    << token::END_STATEMENT << nl
                    << "lastTime        " << lastBatchTime_
                    << token::END_STATEMENT << nl;
            }

            os << "fieldClasses" << nl << token::BEGIN_BLOCK << nl;
            for (const fieldInfo& info : batchFields_)
            {
                os  << "    " << info.name << ' ' << info.className
                    << token::END_STATEMENT << nl;
            }
            os << token::END_BLOCK << nl;
        }

        if (!mv(temporaryFile, stateFile))
        {
            FatalErrorInFunction
                << "Cannot commit batch state " << stateFile
                << exit(FatalError);
        }
    }


    void beginBatch(const std::vector<fieldInfo>& fields)
    {
        ensureSegment();
        batchPath_ = segmentPath_/indexedName("batch_", batchIndex_);

        if (Pstream::master())
        {
            mkDir(batchPath_);
        }
        UPstream::barrier(UPstream::worldComm);

        batchFields_ = fields;
        fieldWriters_.clear();

        for (const fieldInfo& info : batchFields_)
        {
            std::vector<label> dimensions(1, info.localSize);
            if (info.nComponents > 1)
            {
                dimensions.push_back(info.nComponents);
            }

            fieldWriters_.set
            (
                info.name,
                new numpyDetail::numpyFileWriter
                (
                    batchPath_
                   /(info.name + "_proc_"
                   + Foam::name(Pstream::myProcNo()) + ".npy"),
                    dimensions,
                    dtype_
                )
            );
        }

        if (writeTimes_ && Pstream::master())
        {
            timesWriter_.reset
            (
                new numpyDetail::numpyFileWriter
                (
                    batchPath_/"times.npy",
                    std::vector<label>(),
                    numpyDetail::dataType::FLOAT64
                )
            );
        }

        batchCount_ = 0;
        firstBatchTime_ = VGREAT;
        lastBatchTime_ = -VGREAT;
        batchOpen_ = true;
        writeState(false);
    }


    void writeGeometry()
    {
        if
        (
            !geometryDirty_
         || (!writeAreaCentres_ && !writeFaceAreas_ && !writeEdgeCentres_)
        )
        {
            geometryDirty_ = false;
            return;
        }

        ensureSegment();
        const fileName geometryPath
        (
            segmentPath_/indexedName("geometry_", meshRevision_)
        );

        if (Pstream::master())
        {
            mkDir(geometryPath);
        }
        UPstream::barrier(UPstream::worldComm);

        const word suffix
        (
            "_proc_" + Foam::name(Pstream::myProcNo()) + ".npy"
        );

        if (writeAreaCentres_)
        {
            const vectorField& centres = mesh_.areaCentres().primitiveField();
            numpyDetail::numpyFileWriter writer
            (
                geometryPath/("areaCentres" + suffix),
                {centres.size(), vector::nComponents},
                dtype_
            );
            writer.appendField(centres);
        }

        if (writeFaceAreas_)
        {
            const scalarField& areas = mesh_.S().field();
            numpyDetail::numpyFileWriter writer
            (
                geometryPath/("faceAreas" + suffix),
                {areas.size()},
                dtype_
            );
            writer.appendField(areas);
        }

        if (writeEdgeCentres_)
        {
            const vectorField& centres = mesh_.edgeCentres().primitiveField();
            numpyDetail::numpyFileWriter writer
            (
                geometryPath/("edgeCentres" + suffix),
                {centres.size(), vector::nComponents},
                dtype_
            );
            writer.appendField(centres);
        }

        geometryDirty_ = false;
    }


public:

    finiteAreaNumpyWriter
    (
        const word& functionName,
        const Time& time,
        const faMesh& mesh,
        const wordRes& areaSelection,
        const wordRes& edgeSelection,
        const dictionary& dict,
        const fileName& outputPath
    )
    :
        functionName_(functionName),
        time_(time),
        mesh_(mesh),
        areaSelection_(areaSelection),
        edgeSelection_(edgeSelection),
        dtype_
        (
            numpyDetail::parseDataType
            (
                dict.getOrDefault<word>("dataType", "float64")
            )
        ),
        batchSize_(dict.getOrDefault<label>("batchSize", 100)),
        writeTimes_(dict.getOrDefault("writeTimes", true)),
        writeAreaCentres_(dict.getOrDefault("writeAreaCentres", false)),
        writeFaceAreas_(dict.getOrDefault("writeFaceAreas", false)),
        writeEdgeCentres_(dict.getOrDefault("writeEdgeCentres", false)),
        outputPath_(outputPath),
        segmentPath_(),
        batchPath_(),
        segmentReady_(false),
        batchOpen_(false),
        geometryDirty_(true),
        batchIndex_(0),
        batchCount_(0),
        meshRevision_(0),
        firstBatchTime_(VGREAT),
        lastBatchTime_(-VGREAT),
        lastOutputTime_(-VGREAT),
        batchFields_(),
        fieldWriters_(),
        timesWriter_()
    {
        if (batchSize_ <= 0)
        {
            FatalIOErrorInFunction(dict)
                << "batchSize must be greater than zero"
                << exit(FatalIOError);
        }
    }


    ~finiteAreaNumpyWriter()
    {
        end();
    }


    bool write()
    {
        const scalar outputTime = time_.value();

        if (mag(outputTime - lastOutputTime_) <= SMALL)
        {
            return true;
        }
        if (outputTime < lastOutputTime_)
        {
            FatalErrorInFunction
                << "Output time moved backwards from " << lastOutputTime_
                << " to " << outputTime
                << exit(FatalError);
        }

        const std::vector<fieldInfo> fields(selectedFields());

        if (!batchOpen_)
        {
            beginBatch(fields);
        }
        else if (fields.size() != batchFields_.size())
        {
            FatalErrorInFunction
                << "Finite-area field set changed within batch "
                << batchIndex_ << exit(FatalError);
        }

        for (std::size_t i = 0; i < fields.size(); ++i)
        {
            if
            (
                fields[i].name != batchFields_[i].name
             || fields[i].className != batchFields_[i].className
             || fields[i].localSize != batchFields_[i].localSize
            )
            {
                FatalErrorInFunction
                    << "Finite-area field schema changed for "
                    << fields[i].name << exit(FatalError);
            }
        }

        writeGeometry();

        for (const fieldInfo& info : batchFields_)
        {
            bool written = false;

            #define tryFieldType(FieldType)                                  \
                if (!written && info.className == FieldType::typeName)       \
                {                                                            \
                    written = appendField<FieldType>(info.name);             \
                }

            tryFieldType(areaScalarField);
            tryFieldType(areaVectorField);
            tryFieldType(areaSphericalTensorField);
            tryFieldType(areaSymmTensorField);
            tryFieldType(areaTensorField);
            tryFieldType(edgeScalarField);
            tryFieldType(edgeVectorField);
            tryFieldType(edgeSphericalTensorField);
            tryFieldType(edgeSymmTensorField);
            tryFieldType(edgeTensorField);

            #undef tryFieldType

            if (!written)
            {
                FatalErrorInFunction
                    << "Cannot access finite-area field " << info.name
                    << exit(FatalError);
            }
        }

        UPstream::barrier(UPstream::worldComm);
        if (timesWriter_)
        {
            timesWriter_->appendValue(outputTime);
        }
        UPstream::barrier(UPstream::worldComm);

        if (!batchCount_)
        {
            firstBatchTime_ = outputTime;
        }
        ++batchCount_;
        lastBatchTime_ = outputTime;
        lastOutputTime_ = outputTime;
        writeState(false);

        if (batchCount_ == batchSize_)
        {
            end();
        }
        return true;
    }


    void end()
    {
        if (!batchOpen_)
        {
            return;
        }

        UPstream::barrier(UPstream::worldComm);
        writeState(true);
        fieldWriters_.clear();
        timesWriter_.reset();
        batchFields_.clear();
        batchOpen_ = false;
        batchCount_ = 0;
        ++batchIndex_;
    }


    void meshChanged()
    {
        end();
        ++meshRevision_;
        geometryDirty_ = true;
    }
};


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(finiteAreaToNumpy, 0);
    addToRunTimeSelectionTable(functionObject, finiteAreaToNumpy, dictionary);
}
}


// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::functionObjects::finiteAreaToNumpy::finiteAreaToNumpy
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    areaNames_(),
    areaFieldSelection_(),
    edgeFieldSelection_(),
    ownedMeshes_(),
    writers_()
{
    read(dict);
}


Foam::functionObjects::finiteAreaToNumpy::~finiteAreaToNumpy()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void Foam::functionObjects::finiteAreaToNumpy::initialise
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

    ownedMeshes_.resize(areaNames_.size());
    writers_.resize(areaNames_.size());

    fileName basePath;
    if (dict.readIfPresent("outputDir", basePath))
    {
        if (!basePath.isAbsolute())
        {
            basePath = time_.globalPath()/basePath;
        }
    }
    else
    {
        basePath = time_.globalPath()/functionObject::outputPrefix;
        if (mesh_.name() != polyMesh::defaultRegion)
        {
            basePath /= mesh_.name();
        }
        basePath /= name();
    }
    basePath.clean();

    forAll(areaNames_, areai)
    {
        const word& areaName = areaNames_[areai];
        const faMesh* areaMeshPtr =
            registry ? registry->cfindObject<faMesh>(areaName) : nullptr;

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

        writers_.set
        (
            areai,
            new finiteAreaNumpyWriter
            (
                name(),
                time_,
                *areaMeshPtr,
                areaFieldSelection_,
                edgeFieldSelection_,
                dict,
                basePath/areaName
            )
        );
    }
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::finiteAreaToNumpy::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (dict.found("storageOrder"))
    {
        FatalIOErrorInFunction(dict)
            << "finiteAreaToNumpy always writes Fortran-order arrays"
            << exit(FatalIOError);
    }

    if (!writers_.empty())
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

    dict.readIfPresent("areaFields", areaFieldSelection_);
    dict.readIfPresent("edgeFields", edgeFieldSelection_);
    areaFieldSelection_.uniq();
    edgeFieldSelection_.uniq();

    if (areaFieldSelection_.empty() && edgeFieldSelection_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Specify at least one of 'areaFields' or 'edgeFields'"
            << exit(FatalIOError);
    }

    initialise(dict);
    return true;
}


bool Foam::functionObjects::finiteAreaToNumpy::execute()
{
    return true;
}


bool Foam::functionObjects::finiteAreaToNumpy::write()
{
    for (finiteAreaNumpyWriter& writer : writers_)
    {
        writer.write();
    }
    return true;
}


bool Foam::functionObjects::finiteAreaToNumpy::end()
{
    for (finiteAreaNumpyWriter& writer : writers_)
    {
        writer.end();
    }
    return true;
}


void Foam::functionObjects::finiteAreaToNumpy::updateMesh(const mapPolyMesh& mpm)
{
    if (&mpm.mesh() == &mesh_)
    {
        for (finiteAreaNumpyWriter& writer : writers_)
        {
            writer.meshChanged();
        }
    }
}


void Foam::functionObjects::finiteAreaToNumpy::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        for (finiteAreaNumpyWriter& writer : writers_)
        {
            writer.meshChanged();
        }
    }
}


void Foam::functionObjects::finiteAreaToNumpy::readUpdate
(
    const polyMesh::readUpdateState state
)
{
    if (state != polyMesh::UNCHANGED)
    {
        for (finiteAreaNumpyWriter& writer : writers_)
        {
            writer.meshChanged();
        }
    }
}


// ************************************************************************* //
