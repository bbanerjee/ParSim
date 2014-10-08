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


#include <testprograms/Regridders/BNRRegridder.h>

class LBNRRegridder : BNRRegridder
{
  public:
    LBNRRegridder(double tol, IntVector rr) : BNRRegridder(tol,rr){};

    void regrid(const std::vector<std::list<IntVector> >&lflags, std::vector<Region> &patches);

  private:
  
};


void LBNRRegridder::regrid(const std::vector<std::list<IntVector> > &lflags, std::vector<Region> &patches)
{
  patches.resize(0);

  for(size_t p=0;p<lflags.size();p++)
  {
    if(lflags[p].size()>0)
    {
      std::list<IntVector> flags_tmp(lflags[p].begin(),lflags[p].end());
      brsplit(flags_tmp, patches); 
    }
  }
}
    
