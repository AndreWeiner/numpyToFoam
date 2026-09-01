/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tanuj Ravi
    Copyright (C) 2026 Keysight Technologies
    Copyright (C) 2026 Andre Weiner
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

Class
    Foam::functionObjects::numpyDetail::numpyFileWriter

Description
    Incremental Fortran-order NumPy file writer. The final array dimension
    is the appendable output-time dimension.

    The shape uses constant-width integer fields. Rewriting the header after
    an append therefore never moves the binary payload. Each completed append
    leaves a directly loadable NumPy file.

SourceFiles
    numpyFileWriter.C
    numpyFileWriter.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_numpyFileWriter_H
#define Foam_functionObjects_numpyFileWriter_H

#include "Field.H"
#include "fileName.H"
#include "foamNumpyCoreAdapter.H"

#include <cstddef>
#include <fstream>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace numpyDetail
{

/*---------------------------------------------------------------------------*\
                       Class numpyFileWriter Declaration
\*---------------------------------------------------------------------------*/

class numpyFileWriter
{
    // Private Data

        //- Output file
        fileName path_;

        //- Array dimensions preceding the appendable dimension
        std::vector<label> leadingDimensions_;

        //- Output data type
        dataType dtype_;

        //- Binary stream opened for payload and header updates
        std::fstream stream_;

        //- Byte offset of the binary payload
        std::size_t dataStart_;

        //- Number of committed snapshots
        label count_;

        //- Number of values in one snapshot
        std::size_t valuesPerSnapshot_;


    // Private Member Functions

        //- Write one numeric value in the selected precision
        void writeValue(scalar value);

        //- Rewrite the fixed-width header for the current count
        void rewriteHeader();

        //- No copy construct
        numpyFileWriter(const numpyFileWriter&) = delete;

        //- No copy assignment
        void operator=(const numpyFileWriter&) = delete;


public:

    // Constructors

        //- Construct a new appendable NumPy file
        numpyFileWriter
        (
            const fileName& path,
            const std::vector<label>& leadingDimensions,
            dataType dtype
        );


    //- Destructor
    ~numpyFileWriter() = default;


    // Member Functions

        //- Append an OpenFOAM field with components in Fortran order
        template<class Type>
        void appendField(const UList<Type>& field);

        //- Append a single scalar value
        void appendValue(scalar value);

        //- Number of committed snapshots
        label count() const noexcept
        {
            return count_;
        }

        //- Flush the underlying stream
        void flush();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace numpyDetail
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "numpyFileWriter.txx"

#endif

// ************************************************************************* //
