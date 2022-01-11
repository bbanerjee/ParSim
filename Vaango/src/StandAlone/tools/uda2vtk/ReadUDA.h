/*
 * The MIT License
 *
 * Copyright (c) 2014-2022 Parresia Research Limited, New Zealand
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

#include <StandAlone/tools/uda2vis/udaData.h>
#include <Core/DataArchive/DataArchive.h>

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <algorithm>

namespace Vaango {

  class ReadUDA {

  private:

    std::string d_input_uda_name;
    DataArchive* d_archive;
    std::vector<double> d_cycleTimes;

    bool d_useExtraCells;
    int d_currTimeStep;

    std::map<std::string, void_ref_ptr> d_mesh_domains;
    std::map<std::string, void_ref_ptr> d_mesh_boundaries;

  private:

    void  GetLevelAndLocalPatchNumber(int, int&, int&);
    int   GetGlobalDomainNumber(int, int);
    void  CalculateDomainNesting(int, const std::string&);
        
    void  CheckNaNs(int num, double *data, int level, int patch);

    virtual double GetTime( int timestep );
    virtual int    GetNTimesteps( void );

    virtual void   ActivateTimestep( int timestep ); 

    virtual vtkDataSet*   GetMesh( int timestate, int domain, const char * meshname );
    virtual vtkDataArray* GetVar(  int timestate, int domain, const char * varname );
    virtual vtkDataArray* GetVectorVar( int timestate, int domain, const char * varname );

  public:

    ReadUDA();
    ReadUDA(const std::string& input_uda_name);
    ~ReadUDA();

    /////////////////////////////////////////////////////////////////////
    // Open a data archive.
    Uintah::DataArchive* openDataArchive(const std::string& input_uda_name);

    /////////////////////////////////////////////////////////////////////
    // Close a data archive
    void closeDataArchive(DataArchive *archive);

    /////////////////////////////////////////////////////////////////////
    // Get the grid for the current timestep, so we don't have to query
    // it over and over.  
    GridP getGrid(DataArchive *archive, int timeStepNo);

    /////////////////////////////////////////////////////////////////////
    // Get the time for each cycle.
    std::vector<double> getCycleTimes(DataArchive *archive);

    /////////////////////////////////////////////////////////////////////
    // Get all the information that may be needed for the current timestep,
    // including variable/material info, and level/patch info
    TimeStepInfo* getTimeStepInfo(DataArchive *archive, GridP grid, 
                                  int timestep, bool useExtraCells);

    /////////////////////////////////////////////////////////////////////
    // Read the grid data for the given index range
    template<template <typename> class VAR, typename T>
    static GridDataRaw* readGridData(DataArchive *archive,
                                     const Patch *patch,
                                     const LevelP level,
                                     string variable_name,
                                     int material,
                                     int timestep,
                                     int low[3],
                                     int high[3]);

    template<template<typename> class VAR>
    GridDataRaw* getGridDataMainType(DataArchive *archive,
                                     const Patch *patch,
                                     const LevelP level,
                                     string variable_name,
                                     int material,
                                     int timestep,
                                     int low[3],
                                     int high[3],
                                     const Uintah::TypeDescription *subtype);

    GridDataRaw* getGridData(DataArchive *archive,
                             GridP grid,
                             int level_i,
                             int patch_i,
                             string variable_name,
                             int material,
                             int timestep,
                             int low[3],
                             int high[3]);

    /////////////////////////////////////////////////////////////////////
    // Read all the particle data for a given patch.
    template<typename T>
    ParticleDataRaw* readParticleData(DataArchive *archive,
                                      const Patch *patch,
                                      string variable_name,
                                      int material,
                                      int timestep);

    ParticleDataRaw* getParticleData(DataArchive *archive,
                                     GridP grid,
                                     int level_i,
                                     int patch_i,
                                     string variable_name,
                                     int material,
                                     int timestep);

    /////////////////////////////////////////////////////////////////////
    // Utility functions for copying data from Uintah structures into
    // simple arrays.
    void copyIntVector(int to[3], const IntVector &from) {
      to[0]=from[0];  to[1]=from[1];  to[2]=from[2];
    }

    void copyVector(double to[3], const Vector &from) {
      to[0]=from[0];  to[1]=from[1];  to[2]=from[2];
    }

    void copyVector(double to[3], const Point &from) {
      to[0]=from.x();  to[1]=from.y();  to[2]=from.z();
    }

    /////////////////////////////////////////////////////////////////////
    // Utility functions for serializing Uintah data structures into
    // a simple array 
    template <typename T>
    int numComponents() {
      return 1;
    }

    template <>
    int numComponents<Vector>() {
      return 3;
    }

    template <>
    int numComponents<Point>() {
      return 3;
    }

    template <>
    int numComponents<Matrix3>() {
      return 9;
    }

    template <typename T>
    void copyComponents(double *dest, const T &src) {
      (*dest) = (double)src;
    }

    template <>
    void copyComponents<Vector>(double *dest, const Vector &src) {
      dest[0] = (double)src[0];
      dest[1] = (double)src[1];
      dest[2] = (double)src[2];
    }

    template <>
    void copyComponents<Point>(double *dest, const Point &src) {
      dest[0] = (double)src.x();
      dest[1] = (double)src.y();
      dest[2] = (double)src.z();
    }

    template <>
    void copyComponents<Matrix3>(double *dest, const Matrix3 &src) {
      for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
          dest[i*3+j] = (double)src(i,j);
        }
      }
    }
  }; // end class ReadUDA
} // end namespace Vaango

