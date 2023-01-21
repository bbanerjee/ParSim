/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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
 *  timeextract.cc:  extract data at a point over time
 *
 *  Written by:
 *   James L. Bigler
 *   Bryan J. Worthen
 *   Department of Computer Science
 *   University of Utah
 *   March 2003
 *
 */

#include <Core/DataArchive/DataArchive.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Math/Matrix3.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/ShareAssignParticleVariable.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>

#include <Core/Math/MinMax.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/OS/Dir.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <algorithm>


using namespace std;
using namespace Uintah;

bool verbose = false;
bool quiet = false;

void
usage(const std::string& badarg, const std::string& progname)
{
    if(badarg != "")
        std::cerr <<  "Error parsing argument: " << badarg << std::endl;
    std::cerr <<  "Usage: " << progname << " [options] "
         << "-uda <archive file>\n\n";
    std::cerr <<  "Valid options are:\n";
    std::cerr <<  "  -h,      --help\n";
    std::cerr <<  "  -v,      --variable          <variable name>\n";
    std::cerr <<  "  -m,      --material          <material number> [defaults to 0]\n";
    //    std::cerr <<  "  -binary (prints out the data in binary)\n";
    std::cerr <<  "  -tlow,   --timesteplow       [int] (only outputs timestep from int) [defaults to 0]\n";
    std::cerr <<  "  -thigh,  --timestephigh      [int] (only outputs timesteps up to int) [defaults to last timestep]\n";
    std::cerr <<  "  -i,      --index             <i> <j> <k> [intx] cell index [defaults to 0,0,0]\n";
    std::cerr <<  "  -p,      --point             <x> <y> <z> [doubles] point location in physical coordinates \n";
    std::cerr <<  "  -l,      --level             [int] (level index to query range from) [defaults to 0]\n";
    std::cerr <<  "  -o,      --out               <outputfilename> [defaults to stdout]\n";
    std::cerr <<  "  -vv,     --verbose           (prints status of output)\n";
    std::cerr <<  "  -q,      --quiet             (only print data values)\n";
    std::cerr <<  "  -noxml,  --xml-cache-off (turn off XML caching in DataArchive)\n";
    exit(1);
}

// arguments are the dataarchive, the successive arguments are the same as 
// the arguments to archive->query for data values.  Then comes a type 
// dexcription of the variable being queried, and last is an output stream.

template<class T>
void
printData(DataArchive* archive, string& variable_name,
          int material, IntVector& var_id, int levelIndex,
          unsigned long time_step_lower, unsigned long time_step_upper,
          unsigned long output_precision, ostream& out) 
{
  std::vector<int> index;
  std::vector<double> times;

  // query time info from dataarchive
  archive->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());
  if (!quiet) std::cout << "There are " << index.size() << " timesteps:\n";
      
  //------------------------------
  // figure out the lower and upper bounds on the timesteps
  if (time_step_lower >= times.size()) {
    std::cerr <<  "timesteplow must be between 0 and " << times.size()-1 << std::endl;
    exit(1);
  }

  // set default max time value
  if (time_step_upper == (unsigned long)-1) {
    if (verbose)
      std::cout <<"Initializing time_step_upper to "<<times.size()-1<<"\n";
    time_step_upper = times.size() - 1;
  }

  if (time_step_upper >= times.size() || time_step_upper < time_step_lower) {
    std::cerr <<  "timestephigh("<<time_step_upper<<") must be greater than " << time_step_lower 
         << " and less than " << times.size()-1 << std::endl;
    exit(1);
  }
  
  if (!quiet){
    std::cout << "outputting for times["<<time_step_lower<<"] = " << times[time_step_lower]<<" to times["<<time_step_upper<<"] = "<<times[time_step_upper] << std::endl;
  }
  
  // set defaults for output stream
  out.setf(ios::scientific,ios::floatfield);
  out.precision(output_precision);
  
  // for each type available, we need to query the values for the time range, 
  // variable name, and material
  std::vector<T> values;
  try {
    archive->query(values, variable_name, material, var_id, times[time_step_lower], times[time_step_upper], levelIndex);
  } catch (const VariableNotFoundInGrid& exception) {
    std::cerr <<  "Caught VariableNotFoundInGrid Exception: " << exception.message() << std::endl;
    exit(1);
  }

  // Get the cell position
  GridP grid = archive->queryGrid(time_step_lower);
  const LevelP level = grid->getLevel(levelIndex);  
  Point pos(0.0, 0.0, 0.0);
  if (level) {                                       
    pos = level->getCellPosition(var_id);                    
  }

  // Print out data
  for(unsigned int i = 0; i < values.size(); i++) {
    out << times[time_step_lower + i] << "  " << values[i] 
        << " " << pos.x() << " " << pos.y() << " " << pos.z() << std::endl;
  }
} 

int
main(int argc, char** argv)
{
  /*
   * Default values
   */
  //bool do_binary=false;
  bool findCellIndex=false;

  unsigned long time_step_lower = 0;
  // default to be last timestep, but can be set to 0
  unsigned long time_step_upper = (unsigned long)-1;
  unsigned long output_precision = 16;
  string input_uda_name;
  string output_file_name("-");
  IntVector var_id(0,0,0);
  Point var_pt(0,0,0);
  string variable_name;
  int levelIndex = 0;

  // Now the material index is kind of a hard thing.  There is no way
  // to reliably determine a default material.  Materials are defined
  // on the patch for each varialbe, so this subset of materials could
  // change over patches.  We can guess, that there will be a material
  // 0.  This shouldn't cause the program to crash.  It will spit out
  // an exception and exit gracefully.
  int material = 0;

  // Some datasets run out of memory storing the XML data.  This turns
  // off storing that data.
  bool storeXML = true;
  
  /*
   * Parse arguments
   */
  for(int i=1;i<argc;i++){
    string s=argv[i];
    if(s == "-v" || s == "--variable") {
      variable_name = string(argv[++i]);
    } else if (s == "-m" || s == "--material") {
      material = atoi(argv[++i]);
    } else if (s == "-vv" || s == "--verbose") {
      verbose = true;
    } else if (s == "-q" || s == "--quiet") {
      quiet = true;
    } else if (s == "-tlow" || s == "--timesteplow") {
      time_step_lower = strtoul(argv[++i],(char**)nullptr,10);
    } else if (s == "-thigh" || s == "--timestephigh") {
      time_step_upper = strtoul(argv[++i],(char**)nullptr,10);
    } else if (s == "-pr" || s == "--precision") {
      output_precision = strtoul(argv[++i],(char**)nullptr,10);
      if (output_precision > 32) {
        std::cout << "Output precision cannot be larger than 32. Setting precision to 32 \n";        
        output_precision = 32;
      }
      if (output_precision < 1 ) {
        std::cout << "Output precision cannot be less than 1. Setting precision to 16 \n";
        output_precision = 16;
      }      
    } else if (s == "-i" || s == "--index") {
      int x = atoi(argv[++i]);
      int y = atoi(argv[++i]);
      int z = atoi(argv[++i]);
      var_id = IntVector(x,y,z);
    } else if (s == "-p" || s == "--point") {
      double x = atof(argv[++i]);
      double y = atof(argv[++i]);
      double z = atof(argv[++i]);
      var_pt = Point(x,y,z);
      findCellIndex = true;
    } else if (s == "-l" || s == "--level") {
      levelIndex = atoi(argv[++i]);
    } else if( (s == "-h") || (s == "--help") ) {
      usage( "", argv[0] );
    } else if (s == "-uda") {
      input_uda_name = string(argv[++i]);
    } else if (s == "-o" || s == "--out") {
      output_file_name = string(argv[++i]);
//    } else if(s == "-binary") {
//      do_binary=true;
    } else if(s == "-noxml" || s == "--xml-cache-off") {
      storeXML = false;
    } else {
      usage(s, argv[0]);
    }
  }
  
  if(input_uda_name == ""){
    std::cerr <<  "No archive file specified\n";
    usage("", argv[0]);
  }

  try {
    DataArchive* archive = scinew DataArchive(input_uda_name);

    if (!storeXML){
      archive->turnOffXMLCaching();
    }
    
    std::vector<std::string> vars;
    std::vector<const Uintah::TypeDescription*> types;

    archive->queryVariables(vars, types);
    ASSERTEQ(vars.size(), types.size());
    if (verbose) std::cout << "There are " << vars.size() << " variables:\n";
    bool var_found = false;
    unsigned int var_index = 0;
    for (;var_index < vars.size(); var_index++) {
      if (variable_name == vars[var_index]) {
        var_found = true;
        break;
      }
    }
    
    //__________________________________
    //  Bullet proofing
    if (!var_found) {
      std::cerr <<  "\n\n";
      std::cerr <<  "ERROR: Variable \"" << variable_name << "\" was not found.\n";
      std::cerr <<  "If a variable name was not specified try -v [name].\n";
      std::cerr <<  "Possible variable names are:\n";
      var_index = 0;
      for (;var_index < vars.size(); var_index++) {
        std::cout << "vars[" << var_index << "] = " << vars[var_index] << std::endl;
      }
      std::cerr <<  "\n";
      std::cerr <<  "Goodbye!!\n";
      std::cerr <<  "\n";
      exit(-1);
      //      var = vars[0];
    }

    //__________________________________
    //  compute the cell index from the point
    if(findCellIndex){
      std::vector<int> index;
      std::vector<double> times;
      archive->queryTimesteps(index, times);
      ASSERTEQ(index.size(), times.size());

      GridP grid = archive->queryGrid(time_step_lower);
      const LevelP level = grid->getLevel(levelIndex);  
      if (level){                                       
        var_id=level->getCellIndex(var_pt);                    
      }
    }
      
    if (!quiet){
      std::cout << vars[var_index] << ": " << types[var_index]->getName() << " being extracted for material "<<material<<" at index "<<var_id<<endl;
    }
    
    // get type and subtype of data
    const Uintah::TypeDescription* td = types[var_index];
    const Uintah::TypeDescription* subtype = td->getSubType();

    if( subtype == nullptr ) {
      std::cout << "\n";
      std::cout << "An ERROR occurred.  Subtype is nullptr.  Most likely this means that the automatic\n";
      std::cout << "type instantiation is not working... Are you running on a strange architecture?\n";
      std::cout << "Types should be constructed when global static variables of each type are instantiated\n";
      std::cout << "automatically when the program loads.  The registering of the types occurs in:\n";
      std::cout << "src/Core/Disclosure/TypeDescription.cc in register_type() (called from the\n";
      std::cout << "TypeDescription() constructor(s).  However, I'm not quite sure where the variables\n";
      std::cout << "are initially (or in this case not initially) instantiated...  Need to track\n";
      std::cout << "that down and force them to be created... Dd.\n";
      std::cout << "\n";
      exit( 1 );
    }

    // Open output file, call printData with it's ofstream
    // if no output file, call with cout
    ostream *output_stream = &cout;
    if (output_file_name != "-") {
      if (verbose) std::cout << "Opening \""<<output_file_name<<"\" for writing.\n";
      ofstream *output = new ofstream();
      output->open(output_file_name.c_str());
      if (!(*output)) {
        // Error!!
        std::cerr <<  "Could not open "<<output_file_name<<" for writing.\n";
        exit(1);
      }
      output_stream = output;
    } else {
      //output_stream = cout;
    }

  //__________________________________
  //  Now print out the data  
  switch (subtype->getType()) {
  case Uintah::TypeDescription::Type::double_type:
    printData<double>(archive, variable_name, material, var_id, levelIndex,
                      time_step_lower, time_step_upper, output_precision, *output_stream);
    break;
  case Uintah::TypeDescription::Type::float_type:
    printData<float>(archive, variable_name, material, var_id, levelIndex,
                      time_step_lower, time_step_upper, output_precision, *output_stream);
    break;
  case Uintah::TypeDescription::Type::int_type:
    printData<int>(archive, variable_name, material, var_id, levelIndex,
                   time_step_lower, time_step_upper, output_precision, *output_stream);
    break;
  case Uintah::TypeDescription::Type::Vector:
    printData<Vector>(archive, variable_name, material, var_id, levelIndex,
                   time_step_lower, time_step_upper, output_precision, *output_stream);
    break;
  case Uintah::TypeDescription::Type::Matrix3:
  case Uintah::TypeDescription::Type::bool_type:
  case Uintah::TypeDescription::Type::short_int_type:
  case Uintah::TypeDescription::Type::long_type:
  case Uintah::TypeDescription::Type::long64_type:
    std::cerr <<  "Subtype is not implemented\n";
    exit(1);
    break;
  default:
    std::cerr <<  "Unknown subtype\n";
    exit(1);
  }

  // Delete the output file if it was created.
  if (output_file_name != "-") {
    delete((ofstream*)output_stream);
  }

  } catch (Exception& e) {
    std::cerr <<  "Caught exception: " << e.message() << std::endl;
    exit(1);
  } catch(...){
    std::cerr <<  "Caught unknown exception\n";
    exit(1);
  }
} // end main

