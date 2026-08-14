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

Class
    Foam::numpyPostProcessEnvironment

Description
    Runtime-selectable owner of solver-like objects required by function
    objects operating on fields imported by numpyPostProcess.

SourceFiles
    numpyPostProcessEnvironment.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_numpyPostProcessEnvironment_H
#define Foam_numpyPostProcessEnvironment_H

#include "autoPtr.H"
#include "dictionary.H"
#include "HashSet.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
              Class numpyPostProcessEnvironment Declaration
\*---------------------------------------------------------------------------*/

class numpyPostProcessEnvironment
{
protected:

    // Protected Data

        //- Volume mesh
        fvMesh& mesh_;

        //- Names of fields supplied by NumPy
        const wordHashSet importedFields_;

        //- Instance used to initialise solver-owned fields
        const word templateInstance_;

        //- Velocity field when it is not owned by another model
        autoPtr<volVectorField> UPtr_;


    // Protected Member Functions

        //- Find or read the velocity field before constructing solver models
        volVectorField& velocity();


public:

    //- Runtime type information
    TypeName("numpyPostProcessEnvironment");


    //- Runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        numpyPostProcessEnvironment,
        dictionary,
        (
            fvMesh& mesh,
            const dictionary& dict,
            const wordHashSet& importedFields
        ),
        (mesh, dict, importedFields)
    );


    // Constructors

        //- Construct from mesh and dictionary
        numpyPostProcessEnvironment
        (
            fvMesh& mesh,
            const dictionary& dict,
            const wordHashSet& importedFields
        );


    //- Destructor
    virtual ~numpyPostProcessEnvironment() = default;


    // Selectors

        //- Construct the selected environment
        static autoPtr<numpyPostProcessEnvironment> New
        (
            fvMesh& mesh,
            const dictionary& dict,
            const wordHashSet& importedFields
        );


    // Member Functions

        //- Create objects needed while correcting imported patch fields
        virtual void prepare()
        {}

        //- Refresh or construct objects after importing a snapshot
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
