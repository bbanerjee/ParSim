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


#include <StandAlone/tools/puda/rtdata.h>

#include <StandAlone/tools/puda/util.h>

#include <Core/DataArchive/DataArchive.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ShareAssignParticleVariable.h>

#include <Core/OS/Dir.h>

#include <iomanip>
#include <cstdio>

using namespace Uintah;
struct MaterialData {
  std::vector<ShareAssignParticleVariable<double> > pv_double_list;
  std::vector<ShareAssignParticleVariable<float> > pv_float_list;
  std::vector<ShareAssignParticleVariable<Point> > pv_point_list;
  std::vector<ShareAssignParticleVariable<Vector> > pv_vector_list;
  std::vector<ShareAssignParticleVariable<Matrix3> > pv_matrix3_list;
  ShareAssignParticleVariable<Point> p_x;

  MaterialData() {}
  MaterialData(const MaterialData& other)
        : pv_double_list(other.pv_double_list),
          pv_float_list(other.pv_float_list),
          pv_point_list(other.pv_point_list),
          pv_vector_list(other.pv_vector_list),
          pv_matrix3_list(other.pv_matrix3_list),
          p_x(other.p_x) // Use the copy constructor of ShareAssignParticleVariable
  {
  }

  MaterialData& operator=(const MaterialData& other) 
  {
    if (this != &other) {
      pv_double_list = other.pv_double_list;
      pv_float_list = other.pv_float_list;
      pv_point_list = other.pv_point_list;
      pv_vector_list = other.pv_vector_list;
      pv_matrix3_list = other.pv_matrix3_list;
      p_x = other.p_x;
    }
    return *this;
  }
};

/////
///  Helper Functions:
///

// takes a string and replaces all occurances of old with newch
string replaceChar( string s, char old, char newch );
string makeFileName( const string & raydatadir,
                     const string & variable_file,
                     const string & time_file, 
                     const string & patchID_file,
                     const string &  materialType_file );
// Use the setupOutFiles() function to open a pair of files for
// outputing data to the real-time raytracer.
//
// in: pointers to the pointer to the files data and header
//     the file names
// out: inialized files for writing
//      boolean reporting the success of the file creation
bool setupOutFiles(FILE** data, FILE** header, string name, string head);
///
/////

void
Uintah::rtdata( DataArchive * da, CommandLineFlags & clf )
{
  // Create a directory if it's not already there.
  // The exception occurs when the directory is already there
  // and the Dir.create fails.  This exception is ignored. 
  if( clf.raydatadir != "" ) {
    Dir rayDir;
    try {
      rayDir.create( clf.raydatadir );
    }
    catch (Exception& e) {
      std::cerr <<  "Caught exception: " << e.message() << std::endl;
    }
  }

  // Set up the file that contains a list of all the files
  FILE * filelist;
  const string filelistname = clf.raydatadir + string("/") + string("timelist");
  filelist = fopen(filelistname.c_str(),"w");
  if (!filelist) {
    std::cerr <<  "Can't open output file " << filelistname << std::endl;
    abort();
  }

  std::vector<std::string> vars;
  std::vector<int> num_matl;
  std::vector<const Uintah::TypeDescription*> types;
  da->queryVariables(vars, num_matl, types);
  ASSERTEQ(vars.size(), types.size());
  std::cout << "There are " << vars.size() << " variables:\n";
      
  std::vector<int> index;
  std::vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());
  std::cout << "There are " << index.size() << " timesteps:\n";

  std::string time_file;
  std::string variable_file;
  std::string patchID_file;
  std::string materialType_file;
      
  findTimestep_loopLimits( clf.tslow_set, clf.tsup_set, times, clf.time_step_lower, clf.time_step_upper );
      
  // for all timesteps
  for( unsigned long t = clf.time_step_lower; t <= clf.time_step_upper; t++ ) {
    double time = times[t];
     std::ostringstream tempstr_time;
    tempstr_time << std::setprecision(17) << time;
    time_file = replaceChar(string(tempstr_time.str()),'.','_');
    GridP grid = da->queryGrid(t);
    fprintf(filelist,"<TIMESTEP>\n");
    if(clf.do_verbose) {
      std::cout << "time = " << time << std::endl;
    }
    // Create a directory if it's not already there.
    // The exception occurs when the directory is already there
    // and the Dir.create fails.  This exception is ignored. 
    Dir rayDir;
    try {
      rayDir.create(clf.raydatadir + string("/TS_") + time_file);
    }
    catch (Exception& e) {
      std::cerr <<  "Caught directory making exception: " << e.message() << std::endl;
    }
    // for each level in the grid
    for(int l=0;l<grid->numLevels();l++){
      LevelP level = grid->getLevel(l);
          
      // for each patch in the level
      for(Level::const_patch_iterator iter = level->patchesBegin();
          iter != level->patchesEnd(); iter++){
        const Patch* patch = *iter;
         std::ostringstream tempstr_patch;
        tempstr_patch << patch->getID();
        patchID_file = tempstr_patch.str();
        fprintf(filelist,"<PATCH>\n");

        std::vector<MaterialData> material_data_list; 
                    
        // for all vars in one timestep in one patch
        for(int v=0;v<(int)vars.size();v++){
          std::string var = vars[v];
          //cerr << "***********Found variable " << var << "*********\n";
          variable_file = replaceChar(var,'.','_');
          const Uintah::TypeDescription* td = types[v];
          const Uintah::TypeDescription* subtype = td->getSubType();

          ConsecutiveRangeSet matls = da->queryMaterials(var, patch, t);
          // loop over materials
          for(ConsecutiveRangeSet::iterator matlIter = matls.begin();
              matlIter != matls.end(); matlIter++){
            int matl = *matlIter;
             std::ostringstream tempstr_matl;
            tempstr_matl << matl;
            materialType_file = tempstr_matl.str();

            MaterialData material_data;

            if (matl <(int) material_data_list.size())
              material_data = material_data_list[matl];
                
            switch(td->getType()){
            case Uintah::TypeDescription::Type::ParticleVariable:
              if (clf.do_PTvar) {
                switch(subtype->getType()){
                case Uintah::TypeDescription::Type::double_type:
                  {
                    ParticleVariable<double> value;
                    da->query(value, var, matl, patch, t);
                    material_data.pv_double_list.push_back(value);
                  }
                break;
                case Uintah::TypeDescription::Type::float_type:
                  {
                    ParticleVariable<float> value;
                    da->query(value, var, matl, patch, t);
                    material_data.pv_float_list.push_back(value);
                  }
                break;
                case Uintah::TypeDescription::Type::Point:
                  {
                    ParticleVariable<Point> value;
                    da->query(value, var, matl, patch, t);

                    //cout << __LINE__ << ":var = ("<<var<<") for material index "<<matl<<"\n";
                    if (var == "p.x") {
                      //cout << "Found p.x!\n";
                      material_data.p_x.copyPointer(value);
                    } else {
                      material_data.pv_point_list.push_back(value);
                    }
                  }
                break;
                case Uintah::TypeDescription::Type::Vector:
                  {
                    ParticleVariable<Vector> value;
                    da->query(value, var, matl, patch, t);
                    material_data.pv_vector_list.push_back(value);
                  }
                break;
                case Uintah::TypeDescription::Type::Matrix3:
                  {
                    ParticleVariable<Matrix3> value;
                    da->query(value, var, matl, patch, t);
                    material_data.pv_matrix3_list.push_back(value);
                  }
                break;
                default:
                  std::cerr <<  __LINE__ << ":Particle Variable of unknown type: " << subtype->getName() << std::endl;
                  break;
                }
                break;
              }
              break;
            case Uintah::TypeDescription::Type::NCVariable:
              switch(subtype->getType()){
              case Uintah::TypeDescription::Type::double_type:
                {
                  if (clf.do_NCvar_double) {
                    // setup output files
                    const string raydatafile = makeFileName( clf.raydatadir, variable_file, time_file, patchID_file ,materialType_file );
                    FILE* datafile;
                    FILE* headerfile;
                    if (!setupOutFiles(&datafile,&headerfile,raydatafile,string("hdr")))
                      abort();

                    // addfile to filelist
                    fprintf(filelist,"%s\n",raydatafile.c_str());
                    // get the data and write it out
                    double min = 0.0, max = 0.0;
                    NCVariable<double> value;
                    da->query(value, var, matl, patch, t);
                    IntVector dim(value.getHighIndex()-value.getLowIndex());
                    if(dim.x() && dim.y() && dim.z()){
                      NodeIterator iter = patch->getNodeIterator();
                      min=max=value[*iter];
                      for(;!iter.done(); iter++){
                        min=Min(min, value[*iter]);
                        max=Max(max, value[*iter]);
                        float temp_value = (float)value[*iter];
                        fwrite(&temp_value, sizeof(float), 1, datafile);
                      }   
                    }
                        
                    Point b_min = patch->getExtraBox().lower();
                    Point b_max = patch->getExtraBox().upper();
                        
                    // write the header file
                    fprintf(headerfile, "%d %d %d\n",dim.x(), dim.y(), dim.z());
                    fprintf(headerfile, "%f %f %f\n",(float)b_min.x(),(float)b_min.y(),(float)b_min.z());
                    fprintf(headerfile, "%f %f %f\n",(float)b_max.x(),(float)b_max.y(),(float)b_max.z());
                    fprintf(headerfile, "%f %f\n",(float)min,(float)max);

                    fclose(datafile);
                    fclose(headerfile);
                  }
                }
              break;
              case Uintah::TypeDescription::Type::float_type:
                {
                  if (clf.do_NCvar_float) {
                    // setup output files
                    const string raydatafile = makeFileName( clf.raydatadir, variable_file, time_file, patchID_file, materialType_file );
                    FILE* datafile;
                    FILE* headerfile;
                    if (!setupOutFiles(&datafile,&headerfile,raydatafile,string("hdr")))
                      abort();

                    // addfile to filelist
                    fprintf(filelist,"%s\n",raydatafile.c_str());
                    // get the data and write it out
                    float min = 0.0, max = 0.0;
                    NCVariable<float> value;
                    da->query(value, var, matl, patch, t);
                    IntVector dim(value.getHighIndex()-value.getLowIndex());
                    if(dim.x() && dim.y() && dim.z()){
                      NodeIterator iter = patch->getNodeIterator();
                      min=max=value[*iter];
                      for(;!iter.done(); iter++){
                        min=Min(min, value[*iter]);
                        max=Max(max, value[*iter]);
                        float temp_value = value[*iter];
                        fwrite(&temp_value, sizeof(float), 1, datafile);
                      }   
                    }
                        
                    Point b_min = patch->getExtraBox().lower();
                    Point b_max = patch->getExtraBox().upper();
                        
                    // write the header file
                    fprintf(headerfile, "%d %d %d\n",dim.x(), dim.y(), dim.z());
                    fprintf(headerfile, "%f %f %f\n",(float)b_min.x(),(float)b_min.y(),(float)b_min.z());
                    fprintf(headerfile, "%f %f %f\n",(float)b_max.x(),(float)b_max.y(),(float)b_max.z());
                    fprintf(headerfile, "%f %f\n",(float)min,(float)max);

                    fclose(datafile);
                    fclose(headerfile);
                  }
                }
              break;
              case Uintah::TypeDescription::Type::Point:
                {
                  if (clf.do_NCvar_point) {
                    // not implemented at this time
                  }
                }
              break;
              case Uintah::TypeDescription::Type::Vector:
                {
                  if (clf.do_NCvar_vector) {
                    // not implemented at this time
                  }
                }
              break;
              case Uintah::TypeDescription::Type::Matrix3:
                {
                  if (clf.do_NCvar_matrix3) {
                    // not implemented at this time
                  }
                }
              break;
              default:
                std::cerr <<  "NC variable of unknown type: " << subtype->getName() << std::endl;
                break;
              }
              break;
            case Uintah::TypeDescription::Type::CCVariable:
              switch(subtype->getType()){
              case Uintah::TypeDescription::Type::double_type:
                {
                  if (clf.do_CCvar_double) {
                    // setup output files
                    const string raydatafile = makeFileName( clf.raydatadir, variable_file, time_file, patchID_file, materialType_file );
                    FILE* datafile;
                    FILE* headerfile;
                    if (!setupOutFiles(&datafile,&headerfile,raydatafile,string("hdr")))
                      abort();

                    // addfile to filelist
                    fprintf(filelist,"%s\n",raydatafile.c_str());
                    // get the data and write it out
                    double min = 0.0, max = 0.0;
                    CCVariable<double> value;
                    da->query(value, var, matl, patch, t);
                    IntVector dim(value.getHighIndex()-value.getLowIndex());
                    if(dim.x() && dim.y() && dim.z()){
                      NodeIterator iter = patch->getNodeIterator();
                      min=max=value[*iter];
                      for(;!iter.done(); iter++){
                        min=Min(min, value[*iter]);
                        max=Max(max, value[*iter]);
                        float temp_value = (float)value[*iter];
                        fwrite(&temp_value, sizeof(float), 1, datafile);
                      }   
                    }
                        
                    Point b_min = patch->getExtraBox().lower();
                    Point b_max = patch->getExtraBox().upper();
                        
                    // write the header file
                    fprintf(headerfile, "%d %d %d\n",dim.x(), dim.y(), dim.z());
                    fprintf(headerfile, "%f %f %f\n",(float)b_min.x(),(float)b_min.y(),(float)b_min.z());
                    fprintf(headerfile, "%f %f %f\n",(float)b_max.x(),(float)b_max.y(),(float)b_max.z());
                    fprintf(headerfile, "%f %f\n",(float)min,(float)max);

                    fclose(datafile);
                    fclose(headerfile);
                  }
                }
              break;
              case Uintah::TypeDescription::Type::float_type:
                {
                  if (clf.do_CCvar_float) {
                    // setup output files
                    const string raydatafile = makeFileName( clf.raydatadir, variable_file, time_file, patchID_file, materialType_file );
                    FILE* datafile;
                    FILE* headerfile;
                    if (!setupOutFiles(&datafile,&headerfile,raydatafile,string("hdr")))
                      abort();

                    // addfile to filelist
                    fprintf(filelist,"%s\n",raydatafile.c_str());
                    // get the data and write it out
                    float min = 0.0, max = 0.0;
                    CCVariable<float> value;
                    da->query(value, var, matl, patch, t);
                    IntVector dim(value.getHighIndex()-value.getLowIndex());
                    if(dim.x() && dim.y() && dim.z()){
                      NodeIterator iter = patch->getNodeIterator();
                      min=max=value[*iter];
                      for(;!iter.done(); iter++){
                        min=Min(min, value[*iter]);
                        max=Max(max, value[*iter]);
                        float temp_value = value[*iter];
                        fwrite(&temp_value, sizeof(float), 1, datafile);
                      }   
                    }
                        
                    Point b_min = patch->getExtraBox().lower();
                    Point b_max = patch->getExtraBox().upper();
                        
                    // write the header file
                    fprintf(headerfile, "%d %d %d\n",dim.x(), dim.y(), dim.z());
                    fprintf(headerfile, "%f %f %f\n",(float)b_min.x(),(float)b_min.y(),(float)b_min.z());
                    fprintf(headerfile, "%f %f %f\n",(float)b_max.x(),(float)b_max.y(),(float)b_max.z());
                    fprintf(headerfile, "%f %f\n",(float)min,(float)max);

                    fclose(datafile);
                    fclose(headerfile);
                  }
                }
              break;
              case Uintah::TypeDescription::Type::Point:
                {
                  if (clf.do_CCvar_point) {
                    // not implemented at this time
                  }
                }
              break;
              case Uintah::TypeDescription::Type::Vector:
                {
                  if (clf.do_CCvar_vector) {
                    // not implemented at this time
                  }
                }
              break;
              case Uintah::TypeDescription::Type::Matrix3:
                {
                  if (clf.do_CCvar_matrix3) {
                    // not implemented at this time
                  }
                }
              break;
              default:
                std::cerr <<  "CC variable of unknown type: " << subtype->getName() << std::endl;
                break;
              }
              break;
            default:
              std::cerr <<  "Variable of unknown type: " << td->getName() << std::endl;
              break;
            } // end switch(td->getType())
            if (matl < (int)material_data_list.size())
              material_data_list[matl] = material_data;
            else
              material_data_list.push_back(material_data);
          } // end matl
              
        } // end vars
        // after all the variable data has been collected write it out
        if (clf.do_PTvar) {
          FILE* datafile;
          FILE* headerfile;
          //--------------------------------------------------
          // set up the first min/max
          Point min, max;
          std::vector<double> d_min,d_max,f_min,f_max,v_min,v_max,m_min,m_max;
          bool data_found = false;
          int total_particles = 0;
              
          // loops until a non empty material_data set has been
          // found and inialized the mins and maxes
          for(int m = 0; m <(int) material_data_list.size(); m++) {
            // determine the min and max
            MaterialData md = material_data_list[m];
            //cerr << "First md = " << m << std::endl;
            ParticleSubset* pset = md.p_x.getParticleSubset();
            if (!pset) {
              std::cerr <<  __LINE__ << ":No particle location variable found for material index "<<m<<"\n";
              continue;
              abort();
            }
            int numParticles = pset->numParticles();
            if(numParticles > 0){
              ParticleSubset::iterator iter = pset->begin();

              // setup min/max for p.x
              min=max=md.p_x[*iter];
              // setup min/max for all others
              if (clf.do_PTvar_all) {
                for(int i = 0; i <(int) md.pv_double_list.size(); i++) {
                  d_min.push_back(md.pv_double_list[i][*iter]);
                  d_max.push_back(md.pv_double_list[i][*iter]);
                }
                for(int i = 0; i <(int) md.pv_float_list.size(); i++) {
                  f_min.push_back(md.pv_float_list[i][*iter]);
                  f_max.push_back(md.pv_float_list[i][*iter]);
                }
                for(int i = 0; i < (int)md.pv_vector_list.size(); i++) {
                  v_min.push_back(md.pv_vector_list[i][*iter].length());
                  v_max.push_back(md.pv_vector_list[i][*iter].length());
                }
                for(int i = 0; i < (int)md.pv_matrix3_list.size(); i++) {
                  m_min.push_back(md.pv_matrix3_list[i][*iter].Norm());
                  m_max.push_back(md.pv_matrix3_list[i][*iter].Norm());
                }
              }
              // initialized mins/maxes
              data_found = true;
              // setup output files
              const string raydatafile = makeFileName( clf.raydatadir, string(""), time_file, patchID_file,string("") );
              if (!setupOutFiles(&datafile,&headerfile,raydatafile,string("meta")))
                abort();
              // addfile to filelist
              fprintf(filelist,"%s\n",raydatafile.c_str());
                  
              break;
            }
                
          }

          //--------------------------------------------------
          // extract data and write it to a file MaterialData at a time

          if (clf.do_verbose)
            std::cerr <<  "---Extracting data and writing it out  ";
          for(int m = 0; m <(int) material_data_list.size(); m++) {
            MaterialData md = material_data_list[m];
            ParticleSubset* pset = md.p_x.getParticleSubset();
            // a little redundant, but may not have been cought
            // by the previous section
            if (!pset) {
              std::cerr <<  __LINE__ << ":No particle location variable found\n";
              continue;
              abort();
            }
                
            int numParticles = pset->numParticles();
            total_particles+= numParticles;
            if(numParticles > 0){
              ParticleSubset::iterator iter = pset->begin();
              for(;iter != pset->end(); iter++){
                // p_x
                min=Min(min, md.p_x[*iter]);
                max=Max(max, md.p_x[*iter]);
                float temp_value = (float)(md.p_x[*iter]).x();
                fwrite(&temp_value, sizeof(float), 1, datafile);
                temp_value = (float)(md.p_x[*iter]).y();
                fwrite(&temp_value, sizeof(float), 1, datafile);
                temp_value = (float)(md.p_x[*iter]).z();
                fwrite(&temp_value, sizeof(float), 1, datafile);
                if (clf.do_PTvar_all) {
                  // double data
                  for(int i = 0; i <(int) md.pv_double_list.size(); i++) {
                    double value = md.pv_double_list[i][*iter];
                    d_min[i]=Min(d_min[i],value);
                    d_max[i]=Max(d_max[i],value);
                    temp_value = (float)value;
                    fwrite(&temp_value, sizeof(float), 1, datafile);
                  }
                  // float data
                  for(int i = 0; i <(int) md.pv_float_list.size(); i++) {
                    float value = md.pv_float_list[i][*iter];
                    f_min[i]=Min(f_min[i],(double)value);
                    f_max[i]=Max(f_max[i],(double)value);
                    temp_value = value;
                    fwrite(&temp_value, sizeof(float), 1, datafile);
                  }
                  // vector data
                  for(int i = 0; i < (int)md.pv_vector_list.size(); i++) {
                    double value = md.pv_vector_list[i][*iter].length();
                    v_min[i]=Min(v_min[i],value);
                    v_max[i]=Max(v_max[i],value);
                    temp_value = (float)value;
                    fwrite(&temp_value, sizeof(float), 1, datafile);
                  }
                  // matrix3 data
                  for(int i = 0; i < (int)md.pv_matrix3_list.size(); i++) {
                    double value = md.pv_matrix3_list[i][*iter].Norm();
                    m_min[i]=Min(m_min[i],value);
                    m_max[i]=Max(m_max[i],value);
                    temp_value = (float)value;
                    fwrite(&temp_value, sizeof(float), 1, datafile);
                  }
                  if (clf.do_patch) {
                    temp_value = (float)patch->getID();
                    fwrite(&temp_value, sizeof(float), 1, datafile);
                  }
                  if (clf.do_material) {
                    temp_value = (float)m;
                    fwrite(&temp_value, sizeof(float), 1, datafile);
                  }
                }
              }
            }
          }
              
          //--------------------------------------------------
          // write the header file

          if( clf.do_verbose ) {
            std::cerr <<  "---Writing header file\n";
          }
          if (data_found) {
            fprintf(headerfile,"%d\n",total_particles);
            fprintf(headerfile,"%.17g\n",(max.x()-min.x())/total_particles);
            fprintf(headerfile,"%.17g %.17g\n",min.x(),max.x());
            fprintf(headerfile,"%.17g %.17g\n",min.y(),max.y());
            fprintf(headerfile,"%.17g %.17g\n",min.z(),max.z());
            if (clf.do_PTvar_all) {
              for(int i = 0; i < (int)d_min.size(); i++) {
                fprintf(headerfile,"%.17g %.17g\n",d_min[i],d_max[i]);
              }
              for(int i = 0; i < (int)f_min.size(); i++) {
                fprintf(headerfile,"%.17g %.17g\n",f_min[i],f_max[i]);
              }
              for(int i = 0; i < (int)v_min.size(); i++) {
                fprintf(headerfile,"%.17g %.17g\n",v_min[i],v_max[i]);
              }
              for(int i = 0; i < (int)m_min.size(); i++) {
                fprintf(headerfile,"%.17g %.17g\n",m_min[i],m_max[i]);
              }
              if (clf.do_patch) {
                fprintf(headerfile,"%.17g %.17g\n",(float)patch->getID(),(float)patch->getID());
              }
              if (clf.do_material) {
                fprintf(headerfile,"%.17g %.17g\n",0.0,(float)material_data_list.size());
              }
            }
          }
          fclose(datafile);
          fclose(headerfile);
        }
        fprintf(filelist,"</PATCH>\n");
      } // end patch
    } // end level
    fprintf(filelist,"</TIMESTEP>\n");
  } // end timestep
  fclose(filelist);

} // end rtdata()


string
replaceChar(string s, char old, char newch) {
  string result;
  for (int i = 0; i<(int)s.size(); i++)
    if (s[i] == old)
      result += newch;
    else
      result += s[i];
  return result;
}

// given the various parts of the name we piece together the full name
string 
makeFileName( const string & raydatadir,
              const string & variable_file,
              const string & time_file, 
              const string & patchID_file,
              const string &  materialType_file )
{
  string raydatafile;
  if (raydatadir != "")
    raydatafile+= raydatadir + string("/");
  raydatafile+= string("TS_") + time_file + string("/");
  if (variable_file != "")
    raydatafile+= string("VAR_") + variable_file + string(".");
  if (materialType_file != "")
    raydatafile+= string("MT_") + materialType_file + string(".");
  raydatafile+= string("PI_") + patchID_file;
  return raydatafile;
}

bool
setupOutFiles(FILE** data, FILE** header, string name, string head)
{
  FILE * datafile;
  FILE * headerfile;
  const string headername = name + string(".") + head;

  datafile = fopen(name.c_str(),"w");
  if (!datafile) {
    std::cerr <<  "Can't open output file " << name << std::endl;
    return false;
  }
  
  headerfile = fopen(headername.c_str(),"w");
  if (!headerfile) {
    std::cerr <<  "Can't open output file " << headername << std::endl;
    return false;
  }
  
  *data = datafile;
  *header = headerfile;
  return true;
}
