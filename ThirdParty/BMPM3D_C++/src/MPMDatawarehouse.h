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

#ifndef __BRMPM_MPMDATAWAREHOUSE__
#define __BRMPM_MPMDATAWAREHOUSE__

#include <MPMDataTypes.h>
#include <MPMPatchP.h>
//#include <Output.h>
//#include <OutputVTK.h>
//#include <MPMTime.h>
// #include <MPMShapeFunction.h>

#include <map>

namespace BrMPM {

  class MPMDatawarehouse {

  public:
    MPMDatawarehouse();
//     MPMDatawarehouse(Uintah::ProblemSpecP& ps);
    ~MPMDatawarehouse();

    /*
    void saveData(double dt, MaterialSPArray& matlist);

    void dumpData(double dt, MaterialSPArray& matlist);

    bool checkSave(double dt);
    */

    // Add a variable to the data warehouse
    template<typename T>
    void add(const std::string& label, int dwi, T& val);

    // Zero the variable data in the data warehouse
    void zero(const std::string& label, int dwi);

    // Get the particle data
    template<typename T>
    void get(const std::string& label, int dwi, T& val);

    // Update the datawarehouse with new the particle data
    template<typename T>
    void put(const std::string& label, int dwi, T& val);

    // Add particles to the datawarehouse
    void addParticles(const int& dwi, Point3DParticleData& pX, DoubleParticleData& pVol,
                      Vector3DParticleData& pN, DoubleParticleData& density,
                      const int& shSize);

    void createGrid(int dwi, MPMPatchP& patch);

    void zeroGrid(int dwi);

  private:

    // Datawarehouse data map
    // The string contains the join of the variable name and the 
    // material index (dwi)
    int d_id;
    std::map<std::string, MPMVar> d_var;

//     OutputVTK d_out;
//     MPMTime d_time;
//     MPMShapeFunction d_shapefunction;

  }; // end class

}

 // end namespace

#endif
