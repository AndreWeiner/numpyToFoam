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

Class
    Foam::functionObjects::finiteAreaToFoam

Description
    Imports finite-area face and edge fields from finiteAreaToNumpy output.

SourceFiles
    finiteAreaToFoam.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_finiteAreaToFoam_H
#define Foam_functionObjects_finiteAreaToFoam_H

#include "PtrList.H"
#include "fvMesh.H"
#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faMesh;

namespace functionObjects
{

class finiteAreaNumpyImporter;

/*---------------------------------------------------------------------------*\
                   Class finiteAreaToFoam Declaration
\*---------------------------------------------------------------------------*/

class finiteAreaToFoam
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Finite-area region names
        wordList areaNames_;

        //- Finite-area meshes owned for offline operation
        PtrList<faMesh> ownedMeshes_;

        //- Per-area importers
        PtrList<finiteAreaNumpyImporter> importers_;


    // Private Member Functions

        //- Build mesh references and importers
        void initialise(const dictionary& dict);

        //- No copy construct
        finiteAreaToFoam(const finiteAreaToFoam&) = delete;

        //- No copy assignment
        void operator=(const finiteAreaToFoam&) = delete;


public:

    //- Runtime type information
    TypeName("finiteAreaToFoam");


    // Constructors

        //- Construct from Time and dictionary
        finiteAreaToFoam
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct for an object registry
        finiteAreaToFoam
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~finiteAreaToFoam();


    // Member Functions

        //- Read the function-object dictionary
        virtual bool read(const dictionary& dict);

        //- Import at the current time
        virtual bool execute();

        //- Optionally write imported finite-area fields
        virtual bool write();

        //- Explicitly import a time
        bool load(scalar timeValue);

        //- Available committed input times
        scalarList times() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
