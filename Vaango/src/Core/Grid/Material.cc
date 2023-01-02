/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

#include <Core/Grid/Material.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <iostream>
#include <sstream>
#include <string>

namespace Uintah {

Material::Material() {}

Material::Material(ProblemSpecP& ps)
{
  d_mat_subset = nullptr;
  if (ps->getAttribute("name", d_name)) {
    d_have_name = true;
  } else {
    d_have_name = false;
  }
}

Material::~Material()
{
  if (d_mat_subset && d_mat_subset->removeReference()) {
    delete d_mat_subset;
  }
}

ProblemSpecP
Material::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP mat = 0;
  if (d_have_name) {
    mat = ps->appendChild("material");
    mat->setAttribute("name", d_name);
  } else {
    mat = ps->appendChild("material");
  }

  std::stringstream strstream;
  strstream << getDWIndex();
  std::string index_val = strstream.str();
  mat->setAttribute("index", index_val);
  return mat;
}

int
Material::getDWIndex() const
{
  // Return this material's index into the data warehouse
  return d_dwindex;
}

void
Material::setDWIndex(int idx)
{
  d_dwindex = idx;

  ASSERT(!d_mat_subset);
  d_mat_subset = scinew MaterialSubset();

  d_mat_subset->addReference();
  d_mat_subset->add(d_dwindex);
}

} // end namespace Uintah