/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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

#include <CCA/Components/MPM/Crack/PartialEllipticCrack.h>
#include <iostream>

using std::cout;
using std::endl;
using std::vector;

using namespace Uintah;

PartialEllipticCrack::PartialEllipticCrack(ProblemSpecP& ps)
{
  readCrack(ps);
}

PartialEllipticCrack::~PartialEllipticCrack()
{
}



void PartialEllipticCrack::readCrack(ProblemSpecP& pellipse_ps)
{
  // Center,two points on major and minor axes
  Point p;
 
  pellipse_ps->require("center",p);       vertices.push_back(p);
  pellipse_ps->require("point_axis1",p);  vertices.push_back(p);
  pellipse_ps->require("point_axis2",p);  vertices.push_back(p);
  
  // Extent of the partial ellipse (quarter or half)
  pellipse_ps->require("extent",Extent);
    
  // Resolution on circumference
  NCells=1; 
  pellipse_ps->require("resolution_circumference",NCells);
  
  // Crack front segment ID, -1 by default which means all segments are crack front
  CrkFrtSegID=-1; 
  pellipse_ps->get("crack_front_segment_ID",CrkFrtSegID);

}


void PartialEllipticCrack::outputInitialCrackPlane(int i)
{
  
  std::cout << "  * Partial ellipse " << i+1 << " (" << Extent
       << "): meshed by " << NCells
       << " cells on the circumference." << std::endl;
  if(CrkFrtSegID==-1)
    std::cout << "    crack front: on the ellipse circumference" << std::endl;
  else
    std::cout << "    crack front segment ID: " << CrkFrtSegID
         << std::endl;
  std::cout << "    center: " << vertices[0] << std::endl;
  std::cout << "    end point on axis1: " << vertices[1] << std::endl;
  std::cout << "    end point on axis2: " << vertices[2] << std::endl;

}

void PartialEllipticCrack::discretize([[maybe_unused]] int& nstart0,[[maybe_unused]] vector<Point>& cx, 
                           [[maybe_unused]] std::vector<IntVector>& ce,[[maybe_unused]] vector<int>& SegNodes)
{
}

