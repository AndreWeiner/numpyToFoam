// Copyright (C) 2026 Tanuj Ravi
// Copyright (C) 2026 Keysight Technologies
// SPDX-License-Identifier: (GPL-3.0-or-later)

#include "fvCFD.H"
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

struct NpyMeta
{
  fileName path;
  word file;
  word field;       // e.g. "U"
  label procNo{-1}; // e.g. 0
  std::vector<size_t> shape;
  bool is_f8{false};
  bool is_f4{false};
  bool fortranOrder{false};
  std::streamoff dataOffset{0}; // where the raw array bytes start
};
static NpyMeta readNpyMeta(const fileName &filename)
{
  std::ifstream f(filename.c_str(), std::ios::binary);
  if (!f)
    throw std::runtime_error("Cannot open " + std::string(filename));

  char magic[6];
  f.read(magic, 6);

  char v[2];
  f.read(v, 2);

  unsigned short header_len;
  f.read(reinterpret_cast<char *>(&header_len), 2);

  std::string header(header_len, ' ');
  f.read(&header[0], header_len);
  NpyMeta meta;
  meta.path = filename;
  meta.file = filename.name();
  word base = filename.lessExt().name(); // e.g. "U_proc_0"

  label pos = base.rfind('_');                   // last underscore
  meta.procNo = readLabel(base.substr(pos + 1)); // "0"
  word left = base.substr(0, pos);               // "U_proc"
  label pos2 = left.rfind('_');                  // underscore before proc
  meta.field = left.substr(0, pos2);             // "U"

  meta.is_f8 = header.find("<f8") != std::string::npos;
  meta.is_f4 = header.find("<f4") != std::string::npos;

  if (!meta.is_f8 && !meta.is_f4)
    throw std::runtime_error("Unsupported dtype in " + std::string(filename));
  // Parse Fortran order
  if (header.find("fortran_order") != std::string::npos)
  {
    if (header.find("True") != std::string::npos)
      meta.fortranOrder = true;
    else if (header.find("False") != std::string::npos)
      meta.fortranOrder = false;
    else
      throw std::runtime_error("Could not parse fortran_order in " +
                               std::string(filename));
  }
  else
  {
    throw std::runtime_error("Missing fortran_order in header of " +
                             std::string(filename));
  }

  auto p1 = header.find("(");
  auto p2 = header.find(")");
  std::string shape_str = header.substr(p1 + 1, p2 - p1 - 1);

  size_t start = 0;
  while (true)
  {
    auto comma = shape_str.find(",", start);
    std::string token = shape_str.substr(start, comma - start);

    // trim spaces
    token.erase(0, token.find_first_not_of(" \t"));
    token.erase(token.find_last_not_of(" \t") + 1);

    if (!token.empty())
      meta.shape.push_back(std::stoul(token));

    if (comma == std::string::npos)
      break;
    start = comma + 1;
  }

  meta.dataOffset = f.tellg(); // current position is start of raw data
  return meta;
}

static inline size_t prod(const std::vector<size_t> &v)
{
  size_t p = 1;
  for (auto x : v)
    p *= x;
  return p;
}

static inline void require(bool cond, const std::string &msg)
{
  if (!cond)
    throw std::runtime_error(msg);
}

static inline double readAtIndex(std::ifstream &f, const bool is_f8,
                                 const std::streamoff bytePos)
{
  f.seekg(bytePos, std::ios::beg);

  if (is_f8)
  {
    double v;
    f.read(reinterpret_cast<char *>(&v), sizeof(double));
    return v;
  }
  else
  {
    float v;
    f.read(reinterpret_cast<char *>(&v), sizeof(float));
    return static_cast<double>(v);
  }
}


// Templated read of scalar or vector-space fields
template<class Type>
static void readSnapshot
(
    const NpyMeta& meta,
    const label timeIndex,
    Field<Type>& snapshot
)
{
  constexpr bool is_scalar_type = std::is_arithmetic_v<Type>;
  constexpr direction nCmpt = pTraits<Type>::nComponents;
  constexpr int expectedShape = (is_scalar_type) ? 2 : 3;

  // scalar shape : (nCells, nTimes)
  // vectorspace  : (nCells, nComponents, nTimes)

  if (meta.shape.size() != expectedShape)
  {
      throw std::runtime_error
      (
          std::string
          (
              "Expected " + std::string(pTraits<Type>::typeName)
            + " npy with " + std::to_string(expectedShape) + " dims"
          )
      );
  }

  const label nCells = meta.shape[0];
  const label nCmpt_shape = (is_scalar_type ? 1 : meta.shape[1]);
  const label nTimes = (is_scalar_type ? meta.shape[1] : meta.shape[2]);

  // require(nCmpt_shape == nCmpt)
  if (!is_scalar_type && nCmpt_shape != nCmpt)
  {
      throw std::runtime_error
      (
          std::string
          (
              std::string(pTraits<Type>::typeName)
            + " middle dim must be " + std::to_string(nCmpt)
          )
      );
  }
  require(timeIndex >= 0 && timeIndex < nTimes, "timeIndex out of range");

  snapshot.resize(nCells);

  std::ifstream f(meta.path.c_str(), std::ios::binary);
  require(bool(f), "Cannot open " + std::string(meta.path));

  const size_t elemSize = meta.is_f8 ? sizeof(double) : sizeof(float);

  for (label celli = 0; celli < nCells; ++celli)
  {
    auto& elemValue = snapshot[celli];
    std::streamoff idx;

    if constexpr (is_scalar_type)
    {
      if (!meta.fortranOrder)
      {
        // C-order: (nCells, nCmpt, nTimes)
        // idx = ((cell*nCmpt + cmpt)*nTimes + time)
        idx =
        (
            (static_cast<std::streamoff>(celli)) * nTimes
          + timeIndex
        );
      }
      else
      {
        // Fortran-order: (nCells, nCmpt, nTimes)
        // idx = cell + nCells*(comp + nCmpt*time)
        idx =
        (
            static_cast<std::streamoff>(celli)
          + static_cast<std::streamoff>(nCells) *
            (timeIndex)
        );
      }

      double value =
            readAtIndex(f, meta.is_f8, meta.dataOffset + idx * elemSize);

      elemValue = value;
    }
    else
    {
      for (direction cmpt = 0; cmpt < nCmpt; ++cmpt)
      {
        if (!meta.fortranOrder)
        {
          // C-order: (nCells, nCmpt, nTimes)
          // idx = ((cell*nCmpt + cmpt)*nTimes + time)
          idx =
          (
              (static_cast<std::streamoff>(celli) * nCmpt + cmpt) * nTimes
            + timeIndex
          );
        }
        else
        {
          // Fortran-order: (nCells, nCmpt, nTimes)
          // idx = cell + nCells*(comp + nCmpt*time)
          idx =
          (
              static_cast<std::streamoff>(celli)
            + static_cast<std::streamoff>(nCells) *
              (cmpt + static_cast<std::streamoff>(nCmpt) * timeIndex)
          );
        }

        double value =
            readAtIndex(f, meta.is_f8, meta.dataOffset + idx * elemSize);

        elemValue[cmpt] = value;
      }
    }
  }
}


// Templated read of scalar or vector-space fields
template<class Type>
static tmp<Field<Type>> readSnapshot
(
    const NpyMeta& meta,
    const label timeIndex
)
{
  auto tresult = tmp<Field<Type>>::New();
  readSnapshot(meta, timeIndex, tresult.ref());
  return tresult;
}


// ************************************************************************* //
