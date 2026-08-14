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

#include "foamNumpyCoreAdapter.H"
#include "foamNumpyCoreCompat.H"
#include "error.H"

#if OPENFOAM >= 2607
    // Hide the optional header name from older wmkdepend implementations.
    #define FOAM_TO_NUMPY_CORE_HEADER "foamNumpyCore.H"
    #include FOAM_TO_NUMPY_CORE_HEADER
    #define FOAM_TO_NUMPY_HAS_NATIVE_CORE 1
#else
    #define FOAM_TO_NUMPY_HAS_NATIVE_CORE 0
#endif

#include <iomanip>
#include <sstream>

// * * * * * * * * * * * * * * Free Functions  * * * * * * * * * * * * * * //

Foam::functionObjects::numpyDetail::dataType
Foam::functionObjects::numpyDetail::parseDataType(const word& value)
{
    if (value == "float32")
    {
        return dataType::FLOAT32;
    }
    else if (value == "float64")
    {
        return dataType::FLOAT64;
    }

    FatalErrorInFunction
        << "Unsupported dataType " << value << nl
        << "Valid data types are float32 and float64"
        << exit(FatalError);

    return dataType::FLOAT64;
}


const char* Foam::functionObjects::numpyDetail::dataTypeName
(
    dataType dtype
) noexcept
{
    return (dtype == dataType::FLOAT32 ? "float32" : "float64");
}


std::size_t Foam::functionObjects::numpyDetail::dataTypeSize
(
    dataType dtype
) noexcept
{
    return (dtype == dataType::FLOAT32 ? sizeof(float) : sizeof(double));
}


std::string Foam::functionObjects::numpyDetail::dtypeDescription
(
    dataType dtype
)
{
#if FOAM_TO_NUMPY_HAS_NATIVE_CORE
    const auto nativeType =
    (
        dtype == dataType::FLOAT32
      ? Foam::numpy::dataType::type_float32
      : Foam::numpy::dataType::type_float64
    );

    return Foam::numpy::dtypeDescr(nativeType);
#else
    #ifdef WM_BIG_ENDIAN
    constexpr char endian = '>';
    #else
    constexpr char endian = '<';
    #endif

    return std::string(1, endian)
        + (dtype == dataType::FLOAT32 ? "f4" : "f8");
#endif
}


bool Foam::functionObjects::numpyDetail::hasNativeCore() noexcept
{
    return FOAM_TO_NUMPY_HAS_NATIVE_CORE;
}


std::string Foam::functionObjects::numpyDetail::appendableShape
(
    const std::vector<label>& leadingDimensions,
    label count
)
{
    constexpr int width = 20;
    std::ostringstream os;
    os << '(';

    for (const label size : leadingDimensions)
    {
        os << std::setw(width) << size << ',';
    }

    os << std::setw(width) << count << ',' << ')';
    return os.str();
}


std::size_t Foam::functionObjects::numpyDetail::writeHeader
(
    std::ostream& os,
    dataType dtype,
    const std::string& shape
)
{
    const std::string description(dtypeDescription(dtype));

#if FOAM_TO_NUMPY_HAS_NATIVE_CORE
    return Foam::fileFormats::numpyCore::writeHeader
    (
        os,
        description,
        shape,
        true
    );
#else
    return foamNumpyCoreCompat::writeHeader
    (
        os,
        description,
        shape,
        true
    );
#endif
}


// ************************************************************************* //
