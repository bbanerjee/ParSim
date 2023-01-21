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


#include <StandAlone/tools/puda/asci.h>

#include <StandAlone/tools/puda/util.h>

#include <Core/DataArchive/DataArchive.h>
#include <Core/Grid/GridP.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;
using namespace Uintah;


void
Uintah::asci( DataArchive *   da,
              const bool      tslow_set,
              const bool      tsup_set,
              unsigned long & time_step_lower,
              unsigned long & time_step_upper )
{
  std::vector<std::string> vars;
  std::vector<const Uintah::TypeDescription*> types;

  da->queryVariables(vars, types);
  ASSERTEQ(vars.size(), types.size());
  int freq = 1; int ts=1;
      
  std::vector<int> index;
  std::vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());
  if (index.size() == 1) {
    std::cout << "There is only 1 timestep:\n"; 
  }
  else {
    std::cout << "There are " << index.size() << " timesteps:\n";
  }
      
  findTimestep_loopLimits(tslow_set, tsup_set,times, time_step_lower, time_step_upper);
      
  // Loop over time
  for( unsigned long t = time_step_lower; t <= time_step_upper; t++ ) {
    double time = times[t];
    int partnum = 1;
    int num_of_particles = 0;
    std::cout << "timestep " << ts << " inprogress... ";
	
    if( ( ts % freq) == 0 ) {
   		
      // dumps header and variable info to file
      //int variable_count =0;
       std::ostringstream fnum;
      string filename;
      int stepnum=ts/freq;
      fnum << setw(4) << setfill('0') << stepnum;
      string partroot("partout");
      filename = partroot+ fnum.str();
      ofstream partfile(filename.c_str());

      partfile << "TITLE = \"Time Step # " << time <<"\"," << std::endl;
                
      // Code to print out a list of Variables
      partfile << "VARIABLES = ";
	
      GridP grid = da->queryGrid(t);
      int l=0;
      LevelP level = grid->getLevel(l);
      Level::const_patchIterator iter = level->patchesBegin();
      const Patch* patch = *iter;
		
		
      // for loop over variables for name printing
      for(unsigned int v=0;v<vars.size();v++){
        std::string var = vars[v];
	       
        ConsecutiveRangeSet matls= da->queryMaterials(var, patch, t);
        // loop over materials
        for( ConsecutiveRangeSet::iterator matlIter = matls.begin(); matlIter != matls.end(); matlIter++ ) {
          int matl = *matlIter;
          const Uintah::TypeDescription* td = types[v];
          const Uintah::TypeDescription* subtype = td->getSubType();
          switch(td->getType()){
	        
            // The following only accesses particle data
          case Uintah::TypeDescription::Type::ParticleVariable:
            switch(subtype->getType()){
            case Uintah::TypeDescription::Type::double_type:
              {
                ParticleVariable<double> value;
                da->query(value, var, matl, patch, t);
                ParticleSubset* pset = value.getParticleSubset();
		      
                if(pset->numParticles() > 0){
                  ParticleSubset::iterator iter = pset->begin();
			
                  if(matl == 0){
                    partfile << ", \"" << var << "\"";}
                  for(;iter != pset->end(); iter++){
                    num_of_particles++;
                  }
                }
                partnum=num_of_particles;
              }
            break;
            case Uintah::TypeDescription::Type::float_type:
              {
                ParticleVariable<float> value;
                da->query(value, var, matl, patch, t);
                ParticleSubset* pset = value.getParticleSubset();
		      
                if(pset->numParticles() > 0){
                  ParticleSubset::iterator iter = pset->begin();
			
                  if(matl == 0){
                    partfile << ", \"" << var << "\"";}
                  for(;iter != pset->end(); iter++){
                    num_of_particles++;
                  }
                }
                partnum=num_of_particles;
              }
            break;
            case Uintah::TypeDescription::Type::Point:
              {
                ParticleVariable<Point> value;
                da->query(value, var, matl, patch, t);
                ParticleSubset* pset = value.getParticleSubset();
		      
                if(pset->numParticles() > 0 && (matl == 0)){
                  partfile << ", \"" << var << ".x\"" << ", \"" << var <<
                    ".y\"" << ", \"" <<var << ".z\"";
                }
              }
            break;
            case Uintah::TypeDescription::Type::Vector:
              {
                ParticleVariable<Vector> value;
                da->query(value, var, matl, patch, t);
                ParticleSubset* pset = value.getParticleSubset();
                //cout << td->getName() << " over " << pset->numParticles() << " particles\n";
                if(pset->numParticles() > 0 && (matl == 0)){
                  partfile << ", \"" << var << ".x\"" << ", \"" << var <<
                    ".y\"" << ", \"" << var << ".z\"";
                }
              }
            break;
            case Uintah::TypeDescription::Type::Matrix3:
              {
                ParticleVariable<Matrix3> value;
                da->query(value, var, matl, patch, t);
                ParticleSubset* pset = value.getParticleSubset();
                //cout << td->getName() << " over " << pset->numParticles() << " particles\n";
                if(pset->numParticles() > 0 && (matl == 0)){
                  partfile << ", \"" << var << ".1.1\"" << ", \"" << var << ".1.2\"" << ", \"" << var << ".1.3\""
                           << ", \"" << var << ".2.1\"" << ", \"" << var << ".2.2\"" << ", \"" << var << ".2.3\""
                           << ", \"" << var << ".3.1\"" << ", \"" << var << ".3.2\"" << ", \"" << var << ".3.3\"";
                }
              }
            break;
            default:
              std::cerr <<  "Particle Variable of unknown type: " << subtype->getName() << std::endl;
              break;
            }
            break;
          default:
            // Dd: Is this an error!?
            break;
          } // end switch( td->getType() )
		 
        } // end of for loop over materials

        // resets counter of number of particles, so it doesn't count for multiple
        // variables of the same type
        num_of_particles = 0;
	       
      } // end of for loop over variables
		
      partfile << std::endl << "ZONE I=" << partnum << ", F=BLOCK" << std::endl;	
		
      // Loop to print values for specific timestep
      // Because header has already been printed
		
      //variable initialization
      grid = da->queryGrid(t);
      level = grid->getLevel(l);
      iter = level->patchesBegin();
      patch = *iter;
	
      // loop over variables for printing values
      for(unsigned int v=0;v<vars.size();v++){
        std::string var = vars[v];
		
        ConsecutiveRangeSet matls=da->queryMaterials(var, patch, t);
        // loop over materials
        for(ConsecutiveRangeSet::iterator matlIter = matls.begin();
            matlIter != matls.end(); matlIter++){
          int matl = *matlIter;
          const Uintah::TypeDescription* td = types[v];
          const Uintah::TypeDescription* subtype = td->getSubType();
	        
          // the following only accesses particle data
          switch(td->getType()){
          case Uintah::TypeDescription::Type::ParticleVariable:
            switch(subtype->getType()){
            case Uintah::TypeDescription::Type::double_type:
              {
                ParticleVariable<double> value;
                da->query(value, var, matl, patch, t);
                ParticleSubset* pset = value.getParticleSubset();
                if(pset->numParticles() > 0){
                  ParticleSubset::iterator iter = pset->begin();
                  for(;iter != pset->end(); iter++){
                    partfile << value[*iter] << " " << std::endl;
                  }
                  partfile << std::endl;
                }
              }
            break;
            case Uintah::TypeDescription::Type::float_type:
              {
                ParticleVariable<float> value;
                da->query(value, var, matl, patch, t);
                ParticleSubset* pset = value.getParticleSubset();
                if(pset->numParticles() > 0){
                  ParticleSubset::iterator iter = pset->begin();
                  for(;iter != pset->end(); iter++){
                    partfile << value[*iter] << " " << std::endl;
                  }
                  partfile << std::endl;
                }
              }
            break;
            case Uintah::TypeDescription::Type::Point:
              {
                ParticleVariable<Point> value;
                da->query(value, var, matl, patch, t);
                ParticleSubset* pset = value.getParticleSubset();
                if(pset->numParticles() > 0){
                  ParticleSubset::iterator iter = pset->begin();
                  for(;iter != pset->end(); iter++){
                    partfile << value[*iter].x() << " " << std::endl;
                  }
                  partfile << std::endl;
                  iter = pset->begin();
                  for(;iter != pset->end(); iter++){
                    partfile << value[*iter].y() << " " << std::endl;
                  }
                  partfile << std::endl;
                  iter = pset->begin();
                  for(;iter != pset->end(); iter++){
                    partfile << value[*iter].z() << " " << std::endl;
                  }  
                  partfile << std::endl;  
                }
              }
            break;
            case Uintah::TypeDescription::Type::Vector:
              {
                ParticleVariable<Vector> value;
                da->query(value, var, matl, patch, t);
                ParticleSubset* pset = value.getParticleSubset();
                if(pset->numParticles() > 0){
                  ParticleSubset::iterator iter = pset->begin();
                  for(;iter != pset->end(); iter++){
                    partfile << value[*iter].x() << " " << std::endl;
                  }
                  partfile << std::endl;
                  iter = pset->begin();
                  for(;iter != pset->end(); iter++){
                    partfile << value[*iter].y() << " " << std::endl;
                  }
                  partfile << std::endl;
                  iter = pset->begin();
                  for(;iter != pset->end(); iter++){
                    partfile << value[*iter].z() << " " << std::endl;
                  }  
                  partfile << std::endl; 
                }
              }
            break;
            case Uintah::TypeDescription::Type::Matrix3:
              {
                ParticleVariable<Matrix3> value;
                da->query(value, var, matl, patch, t);
                ParticleSubset* pset = value.getParticleSubset();
                if(pset->numParticles() > 0){
                  ParticleSubset::iterator iter = pset->begin();
                  for(;iter != pset->end(); iter++){
                    partfile << (value[*iter])(0,0) << " " << std::endl;
                  }
                  partfile << std::endl;
                  iter = pset->begin();
                  for(;iter !=pset->end(); iter++){
                    partfile << (value[*iter])(0,1) << " " << std::endl;
                  }
                  partfile << std::endl;
                  iter = pset->begin();
                  for(;iter !=pset->end(); iter++){
                    partfile << (value[*iter])(0,2) << " " << std::endl;
                  }
                  partfile << std::endl;
                  iter = pset->begin();
                  for(;iter !=pset->end(); iter++){
                    partfile << (value[*iter])(1,0) << " " << std::endl;
                  }
                  partfile << std::endl;
                  iter = pset->begin();
                  for(;iter !=pset->end(); iter++){
                    partfile << (value[*iter])(1,1) << " " << std::endl;
                  }
                  partfile << std::endl;
                  iter = pset->begin();
                  for(;iter !=pset->end(); iter++){
                    partfile << (value[*iter])(1,2) << " " << std::endl;
                  }
                  partfile << std::endl;
                  iter = pset->begin();
                  for(;iter !=pset->end(); iter++){
                    partfile << (value[*iter])(2,0) << " " << std::endl;
                  }
                  partfile << std::endl;
                  iter = pset->begin();
                  for(;iter !=pset->end(); iter++){
                    partfile << (value[*iter])(2,1) << " " << std::endl;
                  }
                  partfile << std::endl;
                  iter = pset->begin();
                  for(;iter !=pset->end(); iter++){
                    partfile << (value[*iter])(2,2) << " " << std::endl;
                  }
                  partfile << std::endl;
                }
              }
            break;
            default:
              std::cerr <<  "Particle Variable of unknown type: " << subtype->getName() << std::endl;
              break;
            }
            break;
          default:
            // Dd: Is this an error?
            break;
          } // end switch( td->getType() )
        } // end of loop over materials 
      } // end of loop over variables for printing values
    } // end of if ts % freq	

    //increments to next timestep
    ts++;
    std::cout << " completed." << std::endl;
  } // end of loop over time

} // end asci()
