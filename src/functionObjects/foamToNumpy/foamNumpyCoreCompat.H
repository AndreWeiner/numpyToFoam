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
    Foam::functionObjects::numpyDetail::foamNumpyCoreCompat

Description
    Compatibility implementation of the NumPy header routines introduced
    by Foam::fileFormats::numpyCore after OpenFOAM-v2606.

    This class is private to the foamToNumpy function-object library. It is
    not intended to provide or emulate a public OpenFOAM API.

SourceFiles
    foamNumpyCoreCompat.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_foamNumpyCoreCompat_H
#define Foam_functionObjects_foamNumpyCoreCompat_H

#include <cstddef>
#include <iosfwd>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace numpyDetail
{

/*---------------------------------------------------------------------------*\
                     Class foamNumpyCoreCompat Declaration
\*---------------------------------------------------------------------------*/

class foamNumpyCoreCompat
{
public:

    //- Write a version 1.0 NumPy header and return its total byte count
    static std::size_t writeHeader
    (
        std::ostream& os,
        const std::string& dtype,
        const std::string& shape,
        bool fortranOrder
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace numpyDetail
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
