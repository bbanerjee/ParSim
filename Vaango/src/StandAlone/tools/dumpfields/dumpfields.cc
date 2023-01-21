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

/*
 *  dumpfields.cc: Print out a uintah data archive
 *
 *  The fault of 
 *   Andrew D. Brydon
 *   Los Alamos National Laboratory
 *   Mar 2004
 *
 *  Based on puda, written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   February 2000
 *
 */

#include <cassert>

#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Math/Matrix3.h>
#include <Core/DataArchive/DataArchive.h>
#include <Core/Grid/Variables/ShareAssignParticleVariable.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Math/MinMax.h>
#include <Core/Util/Endian.h>
#include <Core/Geometry/Point.h>
#include <Core/Grid/Box.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Geometry/Vector.h>
#include <Core/OS/Dir.h>
//#include <Core/Containers/Array3.h>

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cerrno>
#include <algorithm>

#include "utils.h"
#include "Args.h"
#include "FieldSelection.h"
#include "ScalarDiags.h"
#include "VectorDiags.h"
#include "TensorDiags.h"

#include "TextDumper.h"
#include "InfoDumper.h"
#include "HistogramDumper.h"
/*
#include "EnsightDumper.h"
#include "DXDumper.h"
*/


using namespace std;
using namespace Uintah;


// -----------------------------------------------------------------------------

// store tuple of (variable, it's type)
typedef pair<std::string, const Uintah::TypeDescription*> typed_varname;		  

static 
void usage(const string& badarg, const string& progname)
{
  if(badarg != "")
    std::cerr <<  "Error parsing argument: " << badarg << std::endl;
  std::cerr <<  "Usage: " << progname << " [options] <archive file>\n\n";
  std::cerr <<  "Valid options are:\n";
  std::cerr <<  "  -basename [bnm]            alternate output basename\n";
  std::cerr <<  "  -format [fmt]              output format, one of (text,histogram,info)\n"; // ensight, dx
  std::cerr <<  "                             default is info\n";
  std::cerr <<  "  selection options:" << std::endl;
  std::cerr <<  FieldSelection::options() << std::endl;
  std::cerr <<  "  time options:" << std::endl;
  std::cerr <<  "      -datasetlow  [int]    output data sets starting at [int]\n";
  std::cerr <<  "      -datasethigh [int]    output data sets up to [int]\n";
  std::cerr <<  "      -datasetinc  [int]    output every [int] data sets\n";
  std::cerr <<  "  info options:" << std::endl;
  std::cerr <<  InfoOpts::options() << std::endl;
  std::cerr <<  "  text options:" << std::endl;
  std::cerr <<  TextOpts::options() << std::endl;
  // std::cerr <<  "  ensight options:" << std::endl;
  // std::cerr <<  EnsightOpts::options() << std::endl;
  std::cerr <<  "  histogram options:" << std::endl;
  std::cerr <<  HistogramOpts::options() << std::endl;
  std::cerr <<  "  help options:" << std::endl;
  std::cerr <<  "      -help                  this help\n";
  std::cerr <<  "      -showdiags             print available diagnostic names\n";
  std::cerr <<  "      -showtensorops         print known tensor transformations\n";
  std::cerr <<  "      -showfields            print available field names (requires archive name)\n";
  std::cerr <<  std::endl;
  exit(EXIT_SUCCESS);
}

int
main(int argc, char** argv)
{
  try {
    /*
     * Parse arguments
     */
    Args args(argc, argv);
  
    // throw help early
    if(args.getLogical("help") || args.getLogical("h"))
      usage("", args.progname());
    
    // global options
    //bool do_verbose = args.getLogical("verbose");
  
    // time stepping
    int time_step_lower = args.getInteger("datasetlow",  0);
    int time_step_upper = args.getInteger("datasethigh", INT_MAX);
    int time_step_inc   = args.getInteger("datasetinc",  1);

    // general writing options
    string fmt          = args.getString("format",   "info");
    string basedir      = args.getString("basename", "");
    
    if(args.getLogical("showdiags")) {
      std::cout << "Valid diagnostics: " << std::endl;
      describeScalarDiags(cout);
      describeVectorDiags(cout);
      describeTensorDiags(cout);
      std::cout << std::endl;
      exit(EXIT_SUCCESS);
    }
    
    if(args.getLogical("showtensorops")) {
      std::cout << "Valid tensor operations: " << std::endl;
      describeTensorDiags(cout);
      std::cout << std::endl;
      exit(EXIT_SUCCESS);
    }
    
    string filebase = args.trailing();
    if(filebase=="")
      usage("", args.progname());
    
    if(basedir=="")
      basedir = filebase.substr(0, filebase.find('.'));
    
    std::cout << "filebase: " << filebase << std::endl;
    DataArchive* da = scinew DataArchive(filebase);
    
    // load list of possible variables from the data archive
    std::vector<std::string> allvars;
    std::vector<const Uintah::TypeDescription*> alltypes;
    da->queryVariables(allvars, alltypes);
    ASSERTEQ(allvars.size(), alltypes.size());
    
    if(args.getLogical("showfields")) {
      std::cout << "Valid field names are: " << std::endl;
      for(vector<std::string>::const_iterator vit(allvars.begin());vit!=allvars.end();vit++) {
        if(*vit != "p.x") 
          std::cout << "   " << *vit << std::endl;
      }
      std::cout << std::endl;
      exit(EXIT_SUCCESS);
    }
    
    // select appropriate fields, materials and diagnostics
    FieldSelection fldselection(args, allvars);
    
    // build a specific dumper
    FieldDumper * dumper = 0;
    if(fmt=="text") {
      dumper = scinew TextDumper(da, basedir, args, fldselection);
      /* untested 
    } else if(fmt=="ensight") {
      dumper = scinew EnsightDumper(da, basedir, args, fldselection);
      */
    } else if(fmt=="histogram" || fmt=="hist") {
      dumper = new HistogramDumper(da, basedir, args, fldselection);
      /* untested
    } else if(fmt=="dx" || fmt=="opendx") {
      dumper = scinew DXDumper(da, basedir, binary, onedim);
      */
    } else if(fmt=="info") {
      dumper = new InfoDumper(da, basedir, args, fldselection);
    } else {
      std::cerr <<  "Failed to find match to format '" + fmt + "'" << std::endl;
      usage("", argv[0]);
    }
    
    if(args.hasUnused()) {
      std::cerr <<  "Unused options detected" << std::endl;
      std::vector<std::string> extraargs = args.unusedArgs();
      for(vector<std::string>::const_iterator ait(extraargs.begin());ait!=extraargs.end();ait++)
        {
          std::cerr <<  "    " << *ait << std::endl;
        }
      usage("", argv[0]);
    }
    
    // load list of possible indices and times
    std::vector<int>    index;
    std::vector<double> times;
    da->queryTimesteps(index, times);
    ASSERTEQ(index.size(), times.size());
    std::cout << "There are " << index.size() << " timesteps:\n";
    
    if(time_step_lower<0)                  time_step_lower = 0;
    if(time_step_upper>=(int)index.size()) time_step_upper = (int)index.size()-1;
    if(time_step_inc<=0)                   time_step_inc   = 1;
    
    // build list of (variable, type tuples) for any fields in use
     std::list<typed_varname> dumpvars;
    int nvars = (int)allvars.size();
    for(int i=0;i<nvars;i++) {
      if( fldselection.wantField(allvars[i]) )
        dumpvars.push_back( typed_varname(allvars[i], alltypes[i]) );
    }
    
    for(list<typed_varname>::const_iterator vit(dumpvars.begin());vit!=dumpvars.end();vit++) {
      const string fieldname = vit->first;
      const Uintah::TypeDescription* td  = vit->second;
      dumper->addField(fieldname, td);
    }
    
    // loop over the times
    for(int i=time_step_lower;i<=time_step_upper;i+=time_step_inc) {
      std::cout << index[i] << ": " << times[i] << std::endl;
        
      FieldDumper::Step * step_dumper = dumper->addStep(index[i], times[i], i);
        
      step_dumper->storeGrid();
        
      for(list<typed_varname>::const_iterator vit(dumpvars.begin());vit!=dumpvars.end();vit++) {
        const string fieldname = vit->first;
        if(fieldname=="p.x") continue; // dont work with point field
        
        const Uintah::TypeDescription* td      = vit->second;
        
        step_dumper->storeField(fieldname, td);
      }
      std::cout << std::endl;
	
      dumper->finishStep(step_dumper);
      
      // FIXME: 
      // delete step_dumper;
    }
      
    // delete dumper;
    
  } catch (ProblemSetupException & e) {
    std::cerr <<  std::endl;
    std::cerr <<  "----------------------------------------------------------------------" << std::endl;
    std::cerr <<  std::endl;
    std::cerr <<  "ERROR: " << e.message() << std::endl;
    std::cerr <<  std::endl;
    std::cerr <<  "----------------------------------------------------------------------" << std::endl;
    usage("", argv[0]);
    exit(EXIT_FAILURE);
  } catch (Exception& e) {
    std::cerr <<  "Caught exception: " << e.message() << std::endl;
    abort();
  } catch(...){
    std::cerr <<  "Caught unknown exception\n";
    abort();
  }
  exit(EXIT_SUCCESS);
}
