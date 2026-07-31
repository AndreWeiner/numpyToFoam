/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tanuj Ravi
    Copyright (C) 2026 Keysight Technologies
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    numpyToFoam

Description

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "fvCFD.H"
#include "readData.H"

#include <fstream>
#include <string>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// True if the meta information corresponds to a supported shape
bool is_supported_meta(const NpyMeta& meta)
{
    if (meta.shape.size() == 3)
    {
        // vector-space types
        return
        (
            (meta.shape[1] == 3)
         || (meta.shape[1] == 6)
         || (meta.shape[1] == 9)
        );
    }
    else
    {
        // scalar types
        return (meta.shape.size() == 2);
    }
}

// Helper function - read snapshot and update/write volume field
template<class Type>
void update_from_snapshot
(
    const NpyMeta& meta,
    const label timeIndex,
    GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    #ifdef FULLDEBUG
    if (label len = meta.shape[0]; nCells != fld.size())
    {
        FatalErrorInFunction
            << "Size mismatch. snapshot-size=" << len
            << " field-size=" << fld.size() << " for field "
            << fld.name() << nl
            << Foam::exit(FatalError);
    }
    #endif

    Field<Type> snapshot;
    readSnapshotTemplate(meta, timeIndex, snapshot);

    SubList<Type>(fld.primitiveFieldRef(), fld.size()) =
        SubList<Type>(snapshot, fld.size());

    fld.correctBoundaryConditions();
    fld.write();
}

// Helper function - read or create field
template<class GeoField>
void init_field
(
    const word& fieldName,
    const Time& runTime,
    const fvMesh& mesh,
    HashPtrTable<GeoField>& fields
)
{
    using value_type = typename GeoField::value_type;

    fileName field0Path = runTime.path() / "0" / fieldName;

    autoPtr<GeoField> tmpl;

    if (Foam::isFile(field0Path))
    {
        Info<< "Reading existing 0/" << fieldName
            << " [" << pTraits<value_type>::typeName << "]" << nl;

        tmpl = autoPtr<GeoField>::New
        (
            IOobject
            (
                fieldName,
                "0",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false  // no register
            ),
            mesh
        );
    }
    else
    {
        Info<< "Creating field " << fieldName << " (no 0/ file)"
            << " [" << pTraits<value_type>::typeName << "]" << nl;

        tmpl = autoPtr<GeoField>::New
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false  // no register
            ),
            mesh,
            dimensioned<value_type>(Zero)
        );
    }

    fields.insert
    (
        fieldName,
        autoPtr<GeoField>::New
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true  // Register
            ),
            tmpl()
        )
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption("dict", "file", "Alternative numpyToFoamDict");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << nl;
    runTime.printExecutionTime(Info);
    const word dictName("numpyToFoamDict");

    fileName dictPath;

    if (args.readIfPresent("dict", dictPath))
    {
        // Dictionary specified on the command-line ...

        if (isDir(dictPath))
        {
            dictPath /= dictName;
        }
    }
    else
    {
        // Assume dictionary is to be found in the system directory

        dictPath = runTime.system() / dictName;
    }

    IOobject DictIO(dictPath, runTime, IOobject::MUST_READ, IOobject::NO_WRITE,
                    false);

    if (!DictIO.typeHeaderOk<IOdictionary>(true))
    {
        FatalErrorInFunction << DictIO.objectPath() << nl << exit(FatalError);
    }

    Info << "Reading numpyToFoam settings from " << DictIO.objectRelPath()
         << endl;

    const IOdictionary dict(DictIO);

    wordList fieldNames(dict.get<wordList>("fields"));
    const dictionary &timeDict = dict.subDict("time");
    const scalar t_start = timeDict.get<scalar>("startTime");
    const scalar t_end = timeDict.get<scalar>("endTime");
    const scalar dt = timeDict.get<scalar>("deltaT");
    label procNo = Pstream::myProcNo();

    fileName dataDir(dict.getOrDefault<fileName>("dataDir", "data"));
    if (!dataDir.isAbsolute())
    {
        dataDir = runTime.rootPath() / dataDir;
    }

    fileNameList files = Foam::readDir(dataDir, fileName::FILE);

    label nSteps = floor((t_end - t_start) / dt + 0.5) + 1;
    scalarList timeValues(nSteps);
    for (label i = 0; i < nSteps; ++i)
    {
        timeValues[i] = t_start + i * dt;
    }
    timeValues[nSteps - 1] = t_end;

    // Warn once per field if any of its output files already exist
    for (const word& fieldName : fieldNames)
    {
        bool fieldExists = false;

        forAll(timeValues, timeIndex)
        {
            runTime.setTime(timeValues[timeIndex], timeIndex);
            if (isFile(runTime.path() / runTime.timeName() / fieldName))
            {
                fieldExists = true;
                break;
            }
        }

        if (fieldExists && Pstream::master())
        {
            WarningInFunction << "Field '" << fieldName
                              << "' already exists in one or more"
                              << " time directories — will be overwritten." << nl;
        }
    }

    HashTable<NpyMeta> metaByField; // fieldName -> meta

    bool hasRowMajor = false;
    bool hasSinglePrecision = false;

    for (const word& fieldName : fieldNames)
    {
        word fname = fieldName + "_proc_" + Foam::name(procNo) + ".npy";

        fileName fullPath = dataDir / fieldName / fname;

        if (!isFile(fullPath))
        {
            FatalErrorInFunction << "Missing file: " << fullPath << exit(FatalError);
        }

        NpyMeta meta = readNpyMeta(fullPath);

        if (!meta.fortranOrder)
        {
            hasRowMajor = true;
        }
        if (!meta.is_f8)
        {
            hasSinglePrecision = true;
        }

        // TBD: filter/remove unsupported fields here?

        metaByField.insert(fieldName, meta);

        // ---- PRINT INFO ----
        Info << "Rank " << Pstream::myProcNo() << " | Loaded file: " << meta.file
             << " | Shape: (";

        for (std::size_t d = 0; d < meta.shape.size(); ++d)
        {
            Info << meta.shape[d];
            if (d + 1 < meta.shape.size())
            {
                Info << " x ";
            }
        }

        Info << ")" << " | Type: " << (meta.is_f8 ? "float64" : "float32") << nl;
    }

    Info << nl;

    if (hasRowMajor && Pstream::master())
    {
        WarningInFunction << "Detected row-major (C-order) numpy arrays." << nl
                          << "    Reading snapshots will be inefficient due to "
                             "non-contiguous access."
                          << nl
                          << "    Consider saving data in column-major "
                             "(Fortran-order) for better performance."
                          << nl << endl;
    }
    if (hasSinglePrecision && Pstream::master())
    {
        WarningInFunction
            << "Fields stored as float32 (single precision)." << nl
            << "    It will be converted to OpenFOAM double precision during write."
            << nl << endl;
    }
    runTime.setTime(timeValues[0], 0);

    HashPtrTable<volScalarField> scalarFields;
    HashPtrTable<volVectorField> vectorFields;
    HashPtrTable<volSymmTensorField> symmTensorFields;
    HashPtrTable<volTensorField> tensorFields;

    for (const word& fieldName : fieldNames)
    {
        const NpyMeta &meta = metaByField[fieldName];

        if (meta.shape.size() == 2)
        {
            // scalar
            init_field(fieldName, runTime, mesh, scalarFields);
        }
        else if (meta.shape.size() == 3 && meta.shape[1] == 3)
        {
            // vector
            init_field(fieldName, runTime, mesh, vectorFields);
        }
        else if (meta.shape.size() == 3 && meta.shape[1] == 6)
        {
            // symmTensor
            init_field(fieldName, runTime, mesh, symmTensorFields);
        }
        else if (meta.shape.size() == 3 && meta.shape[1] == 9)
        {
            // tensor
            init_field(fieldName, runTime, mesh, tensorFields);
        }
        else
        {
            WarningInFunction
                << "Unhandled dims " << meta.shape.size()
                << " for " << meta.file << exit(FatalError);
        }
    }
    Info<< nl;

    // autoPtr<functionObjectList> functionsPtr;

    forAll(timeValues, timeIndex)
    {
        runTime.setTime(timeValues[timeIndex], timeIndex);

        for (const word& fieldName : fieldNames)
        {
            const NpyMeta &meta = metaByField[fieldName];

            // Now you can read based on meta.shape.size():
            if (meta.shape.size() == 2)
            {
                // scalar
                update_from_snapshot(meta, timeIndex, *(scalarFields[fieldName]));
            }
            else if (meta.shape.size() == 3 && meta.shape[1] == 3)
            {
                // vector
                update_from_snapshot(meta, timeIndex, *(vectorFields[fieldName]));
            }
            else if (meta.shape.size() == 3 && meta.shape[1] == 6)
            {
                // symmTensor
                update_from_snapshot(meta, timeIndex, *(symmTensorFields[fieldName]));
            }
            else if (meta.shape.size() == 3 && meta.shape[1] == 9)
            {
                // tensor
                update_from_snapshot(meta, timeIndex, *(tensorFields[fieldName]));
            }
            else
            {
                FatalErrorInFunction
                    << "Unsupported dims " << meta.shape.size()
                    << " for " << meta.file << exit(FatalError);
            }
        }

        Info<< "Finished writing fields for time = " << runTime.timeName() << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************** //
