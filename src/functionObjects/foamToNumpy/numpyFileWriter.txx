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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "pTraits.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::numpyDetail::numpyFileWriter::appendField
(
    const UList<Type>& field
)
{
    const label nComponents = pTraits<Type>::nComponents;
    const std::size_t expected =
        static_cast<std::size_t>(field.size())*nComponents;

    if (expected != valuesPerSnapshot_)
    {
        FatalErrorInFunction
            << "Snapshot size changed for " << path_ << nl
            << "Expected " << valuesPerSnapshot_
            << " values but received " << expected
            << exit(FatalError);
    }

    const std::size_t byteOffset =
        dataStart_
      + static_cast<std::size_t>(count_)*valuesPerSnapshot_
       *dataTypeSize(dtype_);

    stream_.seekp(byteOffset, std::ios::beg);

    for (direction cmpti = 0; cmpti < nComponents; ++cmpti)
    {
        for (const Type& value : field)
        {
            writeValue(component(value, cmpti));
        }
    }

    stream_.flush();

    if (!stream_.good())
    {
        FatalErrorInFunction
            << "Failed writing NumPy payload to " << path_
            << exit(FatalError);
    }

    ++count_;
    rewriteHeader();
}


// ************************************************************************* //
