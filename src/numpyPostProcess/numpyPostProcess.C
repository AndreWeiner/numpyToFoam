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

Application
    numpyPostProcess

Description
    Execute OpenFOAM function objects directly on fields imported from
    foamToNumpy output, without intermediate native field files.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "finiteAreaToFoam.H"
#include "functionObjectList.H"
#include "numpyPostProcessEnvironment.H"
#include "numpyToFoam.H"
#include "timeSelector.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Execute function objects on fields imported directly from NumPy"
    );
    argList::addOption
    (
        "dict",
        "file",
        "Alternative numpyPostProcessDict"
    );
    timeSelector::addOptions(false, false);
    #include "addRegionOption.H"

    functionObject::postProcess = true;

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word dictName
    (
        args.getOrDefault<word>("dict", "numpyPostProcessDict")
    );
    IOdictionary driverDict
    (
        IOobject
        (
            dictName,
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (!driverDict.found("input") && !driverDict.found("finiteAreaInput"))
    {
        FatalIOErrorInFunction(driverDict)
            << "Specify 'input', 'finiteAreaInput', or both"
            << exit(FatalIOError);
    }

    wordList importedNames;
    if (driverDict.found("input"))
    {
        importedNames =
            driverDict.subDict("input").get<wordList>("fields");
    }
    const wordHashSet importedFields(importedNames);

    // Declaration order is intentional. The finite-area importer must be
    // destroyed before environment-owned area meshes, whereas the environment
    // must be destroyed before importer-owned volume fields used by its models.
    dictionary environmentDict
    (
        driverDict.subOrEmptyDict("environment")
    );

    if (driverDict.found("input"))
    {
        const dictionary& inputDict = driverDict.subDict("input");
        environmentDict.set
        (
            "templateInstance",
            string
            (
                inputDict.getOrDefault<word>("templateInstance", "0")
            )
        );
    }

    autoPtr<functionObjects::numpyToFoam> importerPtr;
    autoPtr<numpyPostProcessEnvironment> environmentPtr
    (
        numpyPostProcessEnvironment::New
        (
            mesh,
            environmentDict,
            importedFields
        )
    );
    autoPtr<functionObjects::finiteAreaToFoam> areaImporterPtr;

    if (driverDict.found("input"))
    {
        dictionary importerDict(driverDict.subDict("input"));
        importerDict.set("correctBoundaryConditions", false);

        importerPtr.reset
        (
            new functionObjects::numpyToFoam
            (
                "numpyImport",
                mesh,
                importerDict
            )
        );
    }

    if (driverDict.found("finiteAreaInput"))
    {
        areaImporterPtr.reset
        (
            new functionObjects::finiteAreaToFoam
            (
                "finiteAreaImport",
                mesh,
                driverDict.subDict("finiteAreaInput")
            )
        );
    }

    scalarList inputTimes
    (
        importerPtr ? importerPtr->times() : areaImporterPtr->times()
    );

    if (importerPtr && areaImporterPtr)
    {
        const scalarList areaTimes(areaImporterPtr->times());

        if (areaTimes.size() != inputTimes.size())
        {
            FatalErrorInFunction
                << "Volume and finite-area NumPy inputs have different "
                << "time sets" << exit(FatalError);
        }

        forAll(inputTimes, timei)
        {
            if
            (
                mag(areaTimes[timei] - inputTimes[timei])
              > SMALL*max(scalar(1), mag(inputTimes[timei]))
            )
            {
                FatalErrorInFunction
                    << "Volume and finite-area NumPy inputs have different "
                    << "time sets" << exit(FatalError);
            }
        }
    }
    instantList availableTimes(inputTimes.size());

    forAll(inputTimes, timei)
    {
        availableTimes[timei] = instant(inputTimes[timei]);
    }

    const instantList selectedTimes
    (
        timeSelector::select(availableTimes, args, runTime.constant())
    );

    if (selectedTimes.empty())
    {
        FatalErrorInFunction
            << "No NumPy times selected from " << inputTimes
            << exit(FatalError);
    }

    autoPtr<functionObjectList> functionsPtr;

    forAll(selectedTimes, timei)
    {
        runTime.setTime(selectedTimes[timei], timei);
        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (importerPtr)
        {
            importerPtr->load(runTime.value());
        }

        environmentPtr->prepare();

        if (importerPtr)
        {
            importerPtr->correctBoundaryConditions();
        }

        environmentPtr->correct();

        if (areaImporterPtr)
        {
            areaImporterPtr->load(runTime.value());
        }

        if (importerPtr)
        {
            importerPtr->write();
        }

        if (!functionsPtr)
        {
            functionsPtr.reset(new functionObjectList(runTime, driverDict));
            functionsPtr->start();
        }

        functionsPtr->execute();

        if (timei == selectedTimes.size()-1)
        {
            functionsPtr->end();
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
