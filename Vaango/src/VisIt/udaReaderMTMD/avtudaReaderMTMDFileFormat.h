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

/*****************************************************************************
 *
 * Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
 * Produced at the Lawrence Livermore National Laboratory
 * LLNL-CODE-400142
 * All rights reserved.
 *
 * This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
 * full copyright notice is contained in the file COPYRIGHT located at the root
 * of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
 *
 * Redistribution  and  use  in  source  and  binary  forms,  with  or  without
 * modification, are permitted provided that the following conditions are met:
 *
 *  - Redistributions of  source code must  retain the above  copyright notice,
 *    this list of conditions and the disclaimer below.
 *  - Redistributions in binary form must reproduce the above copyright notice,
 *    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
 *    documentation and/or other materials provided with the distribution.
 *  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
 * ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
 * LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
 * DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
 * CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
 * LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
 * OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 *
 *****************************************************************************/

// ************************************************************************* //
//                            avtudaReaderMTMDFileFormat.h                   //
// ************************************************************************* //

#ifndef AVT_udaReaderMTMD_FILE_FORMAT_H
#define AVT_udaReaderMTMD_FILE_FORMAT_H

#include <StandAlone/tools/uda2vis/udaData.h>

#include <avtMTMDFileFormat.h>

#include <vector>
#include <map>
#include <string>

// ****************************************************************************
//  Class: avtudaReaderMTMDFileFormat
//
//  Purpose:
//      Reads in udaReaderMTMD files as a plugin to VisIt.
//
//  Programmer: sshankar -- generated by xml2avt
//  Creation:   Tue May 13 19:02:26 PST 2008
//
// ****************************************************************************

// the DataArchive and GridP types are considered opaque from this library, we dlsym()
// functions from uintah libs to allocate, free, and perform operations on them.
class DataArchive;
class GridP;
class DBOptionsAttributes;

class avtudaReaderMTMDFileFormat : public avtMTMDFileFormat
{
public:
  avtudaReaderMTMDFileFormat( const char * filename, DBOptionsAttributes* attrs);
  virtual           ~avtudaReaderMTMDFileFormat();

  virtual double        GetTime( int timestep );
  virtual int           GetNTimesteps( void );

  virtual const char    *GetType( void )   { return "udaReaderMTMD"; };
  virtual void          ActivateTimestep( int timestep ); 

  virtual vtkDataSet    *GetMesh( int timestate, int domain, const char * meshname );
  virtual vtkDataArray  *GetVar(  int timestate, int domain, const char * varname );
  virtual vtkDataArray  *GetVectorVar( int timestate, int domain, const char * varname );


protected:

  virtual void     PopulateDatabaseMetaData(avtDatabaseMetaData *, int);
  void             ReadMetaData(avtDatabaseMetaData *, int);

  virtual void     *GetAuxiliaryData(const char *var, int,
                                     const char *type, void *args,
                                     DestructorFunction &);

  void             GetLevelAndLocalPatchNumber(int, int&, int&);
  int              GetGlobalDomainNumber(int, int);
  void             CalculateDomainNesting(int, const std::string&);
        
  virtual bool     HasInvariantMetaData(void) const { return false; };
  virtual bool     HasInvariantSIL(void) const { return false; };

  void             AddExpressionsToMetadata(avtDatabaseMetaData *md);
  void             CheckNaNs(int num,double *data,int level,int patch);

  // DATA MEMBERS
  bool useExtraCells;
  int currTimestep;

  //VisIt meshes (see https://visitbugs.ornl.gov/issues/52)
  std::map<std::string, void_ref_ptr> mesh_domains;
  std::map<std::string, void_ref_ptr> mesh_boundaries;

  // data that is not dependent on time
  DataArchive *archive;
  std::vector<double> cycleTimes;

  // data that is dependent on time
  GridP *grid;
  TimestepInfo *stepInfo;
        
  // interface to the uda2vis library
  void  * libHandle;

  DataArchive*     (*openDataArchive)(const std::string&);
  void             (*closeDataArchive)(DataArchive*);

  GridP*           (*getGrid)(DataArchive*, int);
  void             (*releaseGrid)(GridP*);

  std::vector<double>   (*getCycleTimes)(DataArchive*);
  TimestepInfo*    (*getTimestepInfo)(DataArchive*, GridP*, int, bool);

  GridDataRaw*     (*getGridData)(DataArchive*, GridP*, int, int, std::string, int, int, int[3], int[3]);
  ParticleDataRaw* (*getParticleData)(DataArchive*, GridP*, int, int, std::string, int, int);
};

#endif
