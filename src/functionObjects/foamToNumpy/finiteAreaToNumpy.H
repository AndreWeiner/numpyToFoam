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
    Foam::functionObjects::finiteAreaToNumpy

Description
    Writes internal finite-area face and edge fields to batched,
    Fortran-order NumPy files.

SourceFiles
    finiteAreaToNumpy.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_finiteAreaToNumpy_H
#define Foam_functionObjects_finiteAreaToNumpy_H

#include "PtrList.H"
#include "fvMesh.H"
#include "fvMeshFunctionObject.H"
#include "polyMesh.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faMesh;

namespace functionObjects
{

class finiteAreaNumpyWriter;

/*---------------------------------------------------------------------------*\
                  Class finiteAreaToNumpy Declaration
\*---------------------------------------------------------------------------*/

class finiteAreaToNumpy
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Finite-area region names
        wordList areaNames_;

        //- Selected area fields
        wordRes areaFieldSelection_;

        //- Selected edge fields
        wordRes edgeFieldSelection_;

        //- Finite-area meshes owned for offline operation
        PtrList<faMesh> ownedMeshes_;

        //- Per-area output writers
        PtrList<finiteAreaNumpyWriter> writers_;


    // Private Member Functions

        //- Build mesh references and writers
        void initialise(const dictionary& dict);

        //- No copy construct
        finiteAreaToNumpy(const finiteAreaToNumpy&) = delete;

        //- No copy assignment
        void operator=(const finiteAreaToNumpy&) = delete;


public:

    //- Runtime type information
    TypeName("finiteAreaToNumpy");


    // Constructors

        //- Construct from Time and dictionary
        finiteAreaToNumpy
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~finiteAreaToNumpy();


    // Member Functions

        //- Read the function-object dictionary
        virtual bool read(const dictionary& dict);

        //- Execute is a no-op
        virtual bool execute();

        //- Append all selected finite-area fields
        virtual bool write();

        //- Seal partially filled batches
        virtual bool end();

        //- Handle parent mesh topology changes
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Handle parent mesh point motion
        virtual void movePoints(const polyMesh& mesh);

        //- Handle parent mesh changes read from disk
        virtual void readUpdate(const polyMesh::readUpdateState state);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
