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

#include "foamNumpyCoreCompat.H"

#include <cstdint>
#include <ostream>
#include <sstream>

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

std::size_t
Foam::functionObjects::numpyDetail::foamNumpyCoreCompat::writeHeader
(
    std::ostream& os,
    const std::string& dtype,
    const std::string& shape,
    bool fortranOrder
)
{
    constexpr char magic[] = "\x93NUMPY";
    constexpr std::size_t preamble = 10;
    constexpr std::size_t alignment = 16;

    std::ostringstream dict;
    dict
        << "{'descr': '" << dtype
        << "', 'fortran_order': " << (fortranOrder ? "True" : "False")
        << ", 'shape': " << shape << '}';

    std::string header(dict.str());
    const std::size_t total = preamble + header.size() + 1;
    const std::size_t remainder = total % alignment;

    if (remainder)
    {
        header.append(alignment - remainder, ' ');
    }

    header.push_back('\n');

    if (header.size() > UINT16_MAX)
    {
        return 0;
    }

    const std::uint16_t length =
        static_cast<std::uint16_t>(header.size());

    os.write(magic, 6);
    os.put(char(0x01));
    os.put(char(0x00));
    os.put(static_cast<char>(length & 0xFF));
    os.put(static_cast<char>((length >> 8) & 0xFF));
    os.write(header.data(), header.size());

    return preamble + header.size();
}


// ************************************************************************* //
