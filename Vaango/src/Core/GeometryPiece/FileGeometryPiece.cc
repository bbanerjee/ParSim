/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/GeometryPiece/FileGeometryPiece.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/Parallel.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/Endian.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

namespace Uintah {

const std::string FileGeometryPiece::TYPE_NAME = "file";

static std::string
numbered_str(const std::string& s, int is)
{
  std::ostringstream b;
  b << s << is;
  return b.str();
}

FileGeometryPiece::FileGeometryPiece(ProblemSpecP& ps)
{
  d_name = "Unnamed " + TYPE_NAME + " from PS";
  ps->require("name", d_file_name);

  // Set up the vars contained in the file that can be used by
  // the particles
  d_scalar_vars = { "p.volume", "p.temperature", "p.color" };
  d_vector_vars = { "p.externalforce", "p.velocity", "p.fiberdir", "p.rvec1",
                    "p.rvec2",         "p.rvec3",    "p.area" };
  d_tensor_vars = { "p.size" };

  proc0cout << "File Geometry Piece: reading points and";
  for (ProblemSpecP varblock = ps->findBlock("var"); varblock;
       varblock              = varblock->findNextBlock("var")) {
    std::string next_var_name("");
    varblock->get(next_var_name);

    if (std::find(d_scalar_vars.begin(), d_scalar_vars.end(), next_var_name) !=
        d_scalar_vars.end()) {
      d_scalar_vars.push_back(next_var_name);
      proc0cout << " and " << next_var_name;
    } else if (std::find(d_vector_vars.begin(),
                         d_vector_vars.end(),
                         next_var_name) != d_vector_vars.end()) {
      d_vector_vars.push_back(next_var_name);
      proc0cout << " and " << next_var_name;
    } else if (std::find(d_tensor_vars.begin(),
                         d_tensor_vars.end(),
                         next_var_name) != d_tensor_vars.end()) {
      d_tensor_vars.push_back(next_var_name);
      proc0cout << " and " << next_var_name;
    } else {
      proc0cout << " .. unknown variable. Skipping.";
    }
  }
  proc0cout << "\n";

  ps->getWithDefault("format", d_file_format, "text");
  if (d_file_format == "bin") {
    d_file_format = isLittleEndian() ? "lsb" : "msb";
  }

  ps->getWithDefault("usePFS", d_usePFS, true);
  Point min(1e30, 1e30, 1e30), max(-1e30, -1e30, -1e30);

  if (d_usePFS) {
    // We must first read in the min and max from file.0 so
    // that we can determine the BoundingBox for the geometry
    std::string file_name = numbered_str(d_file_name + ".", 0);
    std::ifstream source(file_name.c_str());
    if (!source) {
      std::ostringstream err;
      err << "**ERROR**: Coud not open MPM geometry file " << file_name
          << "\n         Failed to find points file.";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }

    checkFileType(source, d_file_format, file_name);
    read_bbox(source, min, max);
    source.close();

  } else {
    // If we do not use PFS then we should read the entire points file now.
    std::ifstream source(d_file_name.c_str());
    if (!source) {
      std::ostringstream err;
      err << "**ERROR**: Coud not open MPM geometry file " << d_file_name
          << "\n         Failed to find points file.";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }

    checkFileType(source, d_file_format, d_file_name);

    // While the file is read the max and min are updated.
    while (source) {
      read_line(source, min, max);
    }
    source.close();
  }

  Vector fudge(1.e-5, 1.e-5, 1.e-5);
  min   = min - fudge;
  max   = max + fudge;
  d_box = Box(min, max);
}

void
FileGeometryPiece::checkFileType(std::ifstream& source,
                                 std::string& fileType,
                                 std::string& file_name)
{
  int c;
  while ((c = source.get()) != EOF && c <= 127) {
    ;
  }

  if (c == EOF) {
    // the file is ascii
    if (fileType != "text") {
      std::ostringstream err;
      err << "ERROR: opening MPM geometry file (" << file_name + ")\n"
          << "In the ups file you've specified that the file format is bin or "
          << fileType << "\n"
          << "However this is a text file.\n";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }
  } else {
    // the file is binary
    if (fileType != "bin" && fileType != "lsb" && fileType != "msb") {
      std::ostringstream err;
      err << "ERROR: opening MPM geometry file (" << file_name + ")\n"
          << "In the ups file you've specified that the file format is bin or "
          << fileType << "\n"
          << "However this is a binary file.  Please correct this "
             "inconsistency.\n";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }
  }

  // Reset read pointer to beginning of file
  source.clear();
  source.seekg(0, std::ios::beg);

  // Check that file can now be read.
  ASSERT(!source.eof());
}

FileGeometryPiece::FileGeometryPiece(const string& /*file_name*/)
{
  d_name = "Unnamed " + TYPE_NAME + " from file_name";
}

void
FileGeometryPiece::outputHelper(ProblemSpecP& ps) const
{
  ps->appendElement("name", d_file_name);
  ps->appendElement("format", d_file_format);
  for (auto var : d_scalar_vars) {
    ps->appendElement("var", var);
  }
  for (auto var : d_vector_vars) {
    ps->appendElement("var", var);
  }
  for (auto var : d_tensor_vars) {
    ps->appendElement("var", var);
  }
}

GeometryPieceP
FileGeometryPiece::clone() const
{
  return std::make_shared<FileGeometryPiece>(*this);
}

bool
FileGeometryPiece::inside(const Point& p) const
{
  if (p == Max(p, d_box.lower()) && p == Min(p, d_box.upper())) {
    return true;
  } else {
    return false;
  }
}

Box
FileGeometryPiece::getBoundingBox() const
{
  return d_box;
}

void
FileGeometryPiece::read_bbox(std::istream& source, Point& min, Point& max) const
{
  if (d_file_format == "text") {
    source >> min(0) >> min(1) >> min(2) >> max(0) >> max(1) >> max(2);

  } else {
    // FIXME: never changes, should save this !
    const bool iamlittle = isLittleEndian();
    const bool needflip  = (iamlittle && (d_file_format == "msb")) ||
                          (!iamlittle && (d_file_format == "lsb"));
    double t;
    source.read((char*)&t, sizeof(double));
    if (needflip) {
      swapbytes(t);
    }
    min(0) = t;
    source.read((char*)&t, sizeof(double));
    if (needflip) {
      swapbytes(t);
    }
    min(1) = t;
    source.read((char*)&t, sizeof(double));
    if (needflip) {
      swapbytes(t);
    }
    min(2) = t;
    source.read((char*)&t, sizeof(double));
    if (needflip) {
      swapbytes(t);
    }
    max(0) = t;
    source.read((char*)&t, sizeof(double));
    if (needflip) {
      swapbytes(t);
    }
    max(1) = t;
    source.read((char*)&t, sizeof(double));
    if (needflip) {
      swapbytes(t);
    }
    max(2) = t;
  }
}
//______________________________________________________________________
//
bool
FileGeometryPiece::read_line(std::istream& is, Point& xmin, Point& xmax)
{
  Point position{ 0.0, 0.0, 0.0 };

  // CPTI and CPDI pass the size matrix columns containing rvec1, rvec2, rvec3
  Matrix3 size{ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

  if (d_file_format == "text") {

    //__________________________________
    //  TEXT FILE
    read_line_text(is, position, size);

  } else if (d_file_format == "lsb" || d_file_format == "msb") {
    //__________________________________
    //  BINARY FILE
    // read unformatted binary numbers
    read_line_binary(is, position, size);
  }

  // CPDI and CPTI check for negative volumes due to Rvector order
  double vol = size.Determinant();
  if (vol < 0) {
    // switch r2 and r3 in size to get a positive volume
    Matrix3 tmpsize(size(0, 0),
                    size(0, 2),
                    size(0, 1),
                    size(1, 0),
                    size(1, 2),
                    size(1, 1),
                    size(2, 0),
                    size(2, 2),
                    size(2, 1));
    d_tensors.at("p.size").push_back(tmpsize);
  } else {
    // CPTI and CPDI populate psize matrix with Rvectors in columns
    // normalized by the grid spacing for interpolators in ParticleCreator.cc
    d_tensors.at("p.size").push_back(size);
  }

  xmin = Min(xmin, position);
  xmax = Max(xmax, position);
  return true;
}

bool
FileGeometryPiece::read_line_text(std::istream& is,
                                  Point& position,
                                  Matrix3& size)
{
  double x1 = 0.0, x2 = 0.0, x3 = 0.0;
  double v1, v2, v3;


  // line always starts with coordinates
  is >> x1 >> x2 >> x3;
  if (is.eof()) {
    return false; // out of points
  }
  // Particle coordinates
  position = Point(x1, x2, x3);
  d_points.push_back(position);

  for (auto var : d_scalar_vars) {
    if (is >> v1) {
      d_scalars.at(var).push_back(v1);
    }
    if (!is) {
      std::ostringstream err;
      err << "**ERROR** Failed while reading point text point file \n"
          << "Position: " << Point(x1, x2, x3) << "\n";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }
  }

  for (auto var : d_vector_vars) {
    if (is >> v1 >> v2 >> v3) {
      d_vectors.at(var).push_back(Vector(v1, v2, v3));
      if (var == "p.rvec1") {
        size(0, 0) = v1;
        size(1, 0) = v2;
        size(2, 0) = v3;
      } else if (var == "p.rvec2") {
        size(0, 1) = v1;
        size(1, 1) = v2;
        size(2, 1) = v3;
      } else if (var == "p.rvec3") {
        size(0, 2) = v1;
        size(1, 2) = v2;
        size(2, 2) = v3;
      }
    }
    if (!is) {
      std::ostringstream err;
      err << "**ERROR** Failed while reading point text point file \n"
          << "Position: " << Point(x1, x2, x3) << "\n";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }
  }

  return true;
}

bool
FileGeometryPiece::read_line_binary(std::istream& is,
                                    Point& position,
                                    Matrix3& size)
{
  double x1 = 0.0, x2 = 0.0, x3 = 0.0;
  double v[3];

  // never changes, should save this !
  const bool iamlittle = isLittleEndian();
  const bool needflip  = (iamlittle && (d_file_format == "msb")) ||
                        (!iamlittle && (d_file_format == "lsb"));

  is.read((char*)&x1, sizeof(double));

  if (!is) {
    return false; // out of points
  }

  is.read((char*)&x2, sizeof(double));
  is.read((char*)&x3, sizeof(double));
  if (needflip) {
    swapbytes(x1);
    swapbytes(x2);
    swapbytes(x3);
  }
  position = Point(x1, x2, x3);
  d_points.push_back(position);

  for (auto var : d_scalar_vars) {
    if (is.read((char*)&v[0], sizeof(double))) {
      if (needflip) {
        swapbytes(v[0]);
      }
      d_scalars.at(var).push_back(v[0]);
    }
    if (!is) {
      std::ostringstream err;
      err << "**ERROR** Failed while reading point text point file \n"
          << "Position: " << Point(x1, x2, x3) << "\n";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }
  }

  for (auto var : d_vector_vars) {
    if (is.read((char*)&v[0], sizeof(double) * 3)) {
      if (needflip) {
        swapbytes(v[0]);
        swapbytes(v[1]);
        swapbytes(v[2]);
      }
      d_vectors.at(var).push_back(Vector(v[0], v[1], v[2]));
      if (var == "p.rvec1") {
        size(0, 0) = v[0];
        size(1, 0) = v[1];
        size(2, 0) = v[2];
      } else if (var == "p.rvec2") {
        size(0, 1) = v[0];
        size(1, 1) = v[1];
        size(2, 1) = v[2];
      } else if (var == "p.rvec3") {
        size(0, 2) = v[0];
        size(1, 2) = v[1];
        size(2, 2) = v[2];
      }
    }
    if (!is) {
      std::ostringstream err;
      err << "**ERROR** Failed while reading point text point file \n"
          << "Position: " << Point(x1, x2, x3) << "\n";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }
  }

  return true;
}

void
FileGeometryPiece::readPoints(int patchID)
{
  if (d_usePFS) {
    std::ifstream source;

    Point minpt(1e30, 1e30, 1e30);
    Point maxpt(-1e30, -1e30, -1e30);

    std::string file_name;
    char fnum[5];

    sprintf(fnum, ".%d", patchID);
    file_name = d_file_name + fnum;

    source.open(file_name.c_str());

    checkFileType(source, d_file_format, file_name);

    // ignore the first line of the file;
    // this has already been processed
    Point notUsed;
    read_bbox(source, notUsed, notUsed);

    while (source) {
      read_line(source, minpt, maxpt);
    }
  }
}

unsigned int
FileGeometryPiece::createPoints()
{
  std::cerr << "You should be reading points .. not creating them" << std::endl;
  return 0;
}

} // end namespace Uintah