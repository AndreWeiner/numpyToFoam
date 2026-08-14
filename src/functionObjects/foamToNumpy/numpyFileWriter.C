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

#include "numpyFileWriter.H"
#include "OSspecific.H"
#include "error.H"

#include <ios>

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::numpyDetail::numpyFileWriter::numpyFileWriter
(
    const fileName& path,
    const std::vector<label>& leadingDimensions,
    dataType dtype
)
:
    path_(path),
    leadingDimensions_(leadingDimensions),
    dtype_(dtype),
    stream_(),
    dataStart_(0),
    count_(0),
    valuesPerSnapshot_(1)
{
    for (const label size : leadingDimensions_)
    {
        if (size < 0)
        {
            FatalErrorInFunction
                << "Negative NumPy dimension " << size << " for " << path_
                << exit(FatalError);
        }

        valuesPerSnapshot_ *= static_cast<std::size_t>(size);
    }

    mkDir(path_.path());

    stream_.open
    (
        path_.c_str(),
        std::ios::binary | std::ios::in | std::ios::out | std::ios::trunc
    );

    if (!stream_.good())
    {
        FatalErrorInFunction
            << "Cannot open NumPy output file " << path_
            << exit(FatalError);
    }

    dataStart_ = writeHeader
    (
        stream_,
        dtype_,
        appendableShape(leadingDimensions_, count_)
    );

    if (!dataStart_ || !stream_.good())
    {
        FatalErrorInFunction
            << "Cannot write NumPy header to " << path_
            << exit(FatalError);
    }

    stream_.flush();
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::numpyDetail::numpyFileWriter::writeValue
(
    scalar value
)
{
    if (dtype_ == dataType::FLOAT32)
    {
        const float output = static_cast<float>(value);
        stream_.write
        (
            reinterpret_cast<const char*>(&output),
            sizeof(output)
        );
    }
    else
    {
        const double output = static_cast<double>(value);
        stream_.write
        (
            reinterpret_cast<const char*>(&output),
            sizeof(output)
        );
    }
}


void Foam::functionObjects::numpyDetail::numpyFileWriter::rewriteHeader()
{
    stream_.seekp(0, std::ios::beg);

    const std::size_t newDataStart = writeHeader
    (
        stream_,
        dtype_,
        appendableShape(leadingDimensions_, count_)
    );

    stream_.flush();

    if (newDataStart != dataStart_ || !stream_.good())
    {
        FatalErrorInFunction
            << "NumPy header size changed while updating " << path_ << nl
            << "Old payload offset: " << dataStart_ << nl
            << "New payload offset: " << newDataStart
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void Foam::functionObjects::numpyDetail::numpyFileWriter::appendValue
(
    scalar value
)
{
    if (valuesPerSnapshot_ != 1)
    {
        FatalErrorInFunction
            << "Cannot append a scalar to a snapshot containing "
            << valuesPerSnapshot_ << " values in " << path_
            << exit(FatalError);
    }

    const std::size_t byteOffset =
        dataStart_
      + static_cast<std::size_t>(count_)*dataTypeSize(dtype_);

    stream_.seekp(byteOffset, std::ios::beg);
    writeValue(value);
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


void Foam::functionObjects::numpyDetail::numpyFileWriter::flush()
{
    stream_.flush();
}


// ************************************************************************* //
