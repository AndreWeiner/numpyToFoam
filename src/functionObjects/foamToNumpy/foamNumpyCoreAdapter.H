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

Namespace
    Foam::functionObjects::numpyDetail

Description
    Version-independent adapter for the NumPy header support used by the
    foamToNumpy function object.

SourceFiles
    foamNumpyCoreAdapter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_foamNumpyCoreAdapter_H
#define Foam_functionObjects_foamNumpyCoreAdapter_H

#include "label.H"

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace numpyDetail
{

//- Numeric output types supported by foamToNumpy
enum class dataType
{
    FLOAT32,
    FLOAT64
};

//- Parse a function-object dataType entry
dataType parseDataType(const word& value);

//- Return the dictionary spelling for a data type
const char* dataTypeName(dataType dtype) noexcept;

//- Return the byte width of a data type
std::size_t dataTypeSize(dataType dtype) noexcept;

//- Return the NumPy dtype descriptor for the native machine endian
std::string dtypeDescription(dataType dtype);

//- True when the adapter was compiled with OpenFOAM's native numpyCore
bool hasNativeCore() noexcept;

//- Build a constant-width shape with an appendable final time dimension
std::string appendableShape
(
    const std::vector<label>& leadingDimensions,
    label count
);

//- Write a Fortran-order NumPy header and return its total byte count
std::size_t writeHeader
(
    std::ostream& os,
    dataType dtype,
    const std::string& shape
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace numpyDetail
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
