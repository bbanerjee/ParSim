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
//                            udaFileFormat.C                   //
// ************************************************************************* //

#include <VisIt/udaReaderMTMD/udaFileFormat.h>

#include <iostream>
#include <string>
#include <cmath>
#include <string>
#include <vector>
#include <limits> 

#ifdef PARALLEL
#  include <mpi.h>
#endif

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkUnstructuredGrid.h>

#include <DebugStream.h>
#include <InvalidVariableException.h>
#include <InvalidFilesException.h>

using std::max;
using std::min;
using std::vector;
using std::string;

//#ifdef PARALLEL
// we only get the option for  serialized reads if compiling the parallel version
//#  define SERIALIZED_READS
//#endif

const double NAN_REPLACE_VAL=1.0E9;

// ****************************************************************************
//  Method: udaFileFormat constructor
//
//  Programmer: sshankar 
//  Creation:   Tue May 13 19:02:26 PST 2008
//
// ****************************************************************************
udaFileFormat::udaFileFormat(const char *filename)
  useExtraCells(true),
  archive(nullptr),
  grid(nullptr)
{
  for (int i = 0; attrs != 0 && i<attrs->GetNumberOfOptions(); ++i) {
    if (attrs->GetName(i) == "Load extra cells") {
      useExtraCells = attrs->GetBool("Load extra cells");
    }
  }
    
  // Verify that it is a UDA index.xml file:
  // The 2nd line should look like this <Uintah_DataArchive>.
  FILE * fp = fopen( filename, "r" );
  if( fp == nullptr ) {
    std::ostringstream error;
    error << "Failed to open file: " << filename;
    throw InvalidFilesExpection(error.str(), __FILE__, __LINE__);
  }

  char line[1024];
  char * result = fgets( line, 1024, fp );
  if( result ) { 
    result = fgets( line, 1024, fp );
  }

  string lineStr = line;
  if( !result || lineStr.find( "<Uintah_DataArchive>" ) == string::npos ) {
    std::ostringstream error;
    error << std::string(filename) << " does not appear to be a <Uintah_DataArchive>.";
    throw InvalidFilesExpection(error.str(), __FILE__, __LINE__);
  }
  fclose( fp );

  // use the folder name, not the index.xml file name to open the archive
  string folder(filename);
  size_t found = folder.find_last_of("/");
  folder = folder.substr(0, found);
  archive = openDataArchive(folder);

  // timestep times
  cycleTimes = getCycleTimes(archive);

  // haven't loaded any timestep data yet
  stepInfo = nullptr;  
  currTimeStep = -1;
}


// Destructor
udaFileFormat::~udaFileFormat()
{
  if (grid) {
    releaseGrid(grid);
  }

  if (archive) {
    closeDataArchive(archive);
  }

  if (stepInfo) {
    delete stepInfo;
  }
}


// ****************************************************************************
//  Method: udaFileFormat::GetNTimesteps
//
//  Purpose:
//      Tells the rest of the code how many timesteps there are in this file.
//
//  Programmer: sshankar 
//  Creation:   Tue May 13 19:02:26 PST 2008
//
// ****************************************************************************

int
udaFileFormat::GetNTimesteps(void)
{
  return cycleTimes.size();
}


// ****************************************************************************
// Method: udaFileFormat::GetTime
//
// Purpose: 
//   Get the time.
//
// Programmer: sshankar 
// Creation:   Fri Feb 6 15:31 MST 2009 
//
// ****************************************************************************

double 
udaFileFormat::GetTime(int ts)
{
  return cycleTimes[ts];
}


// ****************************************************************************
//  Method: getBounds
//
//  Purpose:
//   Returns the bounds for the given patch of the specified mesh 
//   based on periodicity and type.
//
//  Node centered data uses the same mesh as cell centered, 
//  but face centered meshes need an extra value for one axis,
//  unless they are periodic on that axis.
//
//  use patch_id=-1 to query all patches.
//
// ****************************************************************************
void getBounds(int low[3], int high[3], const string meshName, 
               const LevelInfo &levelInfo,int patch_id=-1)
{
  levelInfo.getBounds(low,high,meshName,patch_id);
  
  //debug5<<"getBounds("<<meshName<<",id="<<patch_id<<")=["<<low[0]<<","<<low[1]<<","<<low[2]<<"] to ["<<high[0]<<","<<high[1]<<","<<high[2]<<"]\n";
}

// ****************************************************************************
//  Method: udaFileFormat::ReadMetaData
//
//  Purpose:
//      Does the actual work for PopulateMetaData().
//
//     ADDS MESHES
// ****************************************************************************
void
udaFileFormat::ReadMetaData(avtDatabaseMetaData *md, int timeState)
{
  ActivateTimestep(timeState);

  int numLevels = stepInfo->levelInfo.size();

  int totalPatches = 0;
  for (int i = 0; i < numLevels; i++)
    totalPatches +=  stepInfo->levelInfo[i].patchInfo.size();
  //debug5 << "udaFileFormat::ReadMetaData: Levels: " << numLevels << " Patches: " << totalPatches << endl;

  vector<int> groupIds(totalPatches);
  vector<string> pieceNames(totalPatches);

  for (int i = 0; i < totalPatches; i++) {
    char tmpName[64];
    int level, local_patch;

    GetLevelAndLocalPatchNumber(i, level, local_patch);
    sprintf(tmpName,"level%d, patch%d", level, local_patch);

    groupIds[i] = level;
    pieceNames[i] = tmpName;
  }

  // compute the bounding box of the mesh from the grid indices of level 0
  LevelInfo &levelInfo = stepInfo->levelInfo[0];

  // don't add proc id unless CC_Mesh or NC_Mesh exists (some only have SFCk_MESH)
  bool addProcId=false;
  string mesh_for_procid("NC_Mesh");

  // grid meshes are shared between materials, and particle meshes are
  // shared between variables - keep track of what has been added so they're only added once
  std::set<string> meshes_added;

  // If a variable exists in multiple materials, we don't want to add it more than
  // once to the meta data - it can mess up visit's expressions variable lists.
  std::set<string> mesh_vars_added;

  //get CC bounds
  int low[3],high[3];
  getBounds(low,high,"CC_Mesh",levelInfo);

  //this can be done once for everything because the spatial range is the same for all meshes
  double box_min[3] = { levelInfo.anchor[0] + low[0] * levelInfo.spacing[0],
                        levelInfo.anchor[1] + low[1] * levelInfo.spacing[1],
                        levelInfo.anchor[2] + low[2] * levelInfo.spacing[2] };
  double box_max[3] = { levelInfo.anchor[0] + high[0] * levelInfo.spacing[0],
                        levelInfo.anchor[1] + high[1] * levelInfo.spacing[1],
                        levelInfo.anchor[2] + high[2] * levelInfo.spacing[2] };
  //debug5<<"box_min/max=["<<box_min[0]<<","<<box_min[1]<<","<<box_min[2]<<"] to ["<<box_max[0]<<","<<box_max[1]<<","<<box_max[2]<<"]\n";

  int logical[3];
  for (int i=0; i<3; i++)
    logical[i] = high[i]-low[i];
  //debug5 <<"logical: "<< logical[0] << ", "<< logical[1] << ", "<< logical[2]<<endl;

  for (int i=0; i<(int)stepInfo->varInfo.size(); i++) {
    if (stepInfo->varInfo[i].type.find("ParticleVariable") == string::npos) {
      string varname = stepInfo->varInfo[i].name;
      string vartype = stepInfo->varInfo[i].type;

      string mesh_for_this_var;
      avtCentering cent=AVT_ZONECENT;

      if (vartype.find("NC") != string::npos) {
        cent = AVT_NODECENT;
        mesh_for_this_var.assign("NC_Mesh"); 
        addProcId=true;
      }  
      else if (vartype.find("CC") != string::npos) {  
        cent = AVT_ZONECENT;
        mesh_for_this_var.assign("CC_Mesh");
        addProcId=true;
        mesh_for_procid=mesh_for_this_var;
      }
      else if (vartype.find("SFC") != string::npos) { 
        cent = AVT_ZONECENT;

        if (vartype.find("SFCX") != string::npos)		
          mesh_for_this_var.assign("SFCX_Mesh");
        else if (vartype.find("SFCY") != string::npos)		
          mesh_for_this_var.assign("SFCY_Mesh");
        else if (vartype.find("SFCZ") != string::npos)		
          mesh_for_this_var.assign("SFCZ_Mesh");
      }  
      else
        debug5<<"unknown vartype: "<<vartype<<endl;

      if (meshes_added.find(mesh_for_this_var)==meshes_added.end()) {
        avtMeshMetaData *mesh = new avtMeshMetaData;

        mesh->name = mesh_for_this_var;
        mesh->meshType = AVT_AMR_MESH;
        mesh->topologicalDimension = 3;
        mesh->spatialDimension = 3;

        mesh->numBlocks = totalPatches;
        mesh->blockTitle = "patches";
        mesh->blockPieceName = "patch";
        mesh->numGroups = numLevels;
        mesh->groupTitle = "levels";
        mesh->groupPieceName = "level";
        mesh->blockNames = pieceNames;
        mesh->containsExteriorBoundaryGhosts = false;

        mesh->hasSpatialExtents = true; 
        mesh->minSpatialExtents[0] = box_min[0];
        mesh->maxSpatialExtents[0] = box_max[0];
        mesh->minSpatialExtents[1] = box_min[1];
        mesh->maxSpatialExtents[1] = box_max[1];
        mesh->minSpatialExtents[2] = box_min[2];
        mesh->maxSpatialExtents[2] = box_max[2];

        mesh->hasLogicalBounds = true;
        mesh->logicalBounds[0] = logical[0];
        mesh->logicalBounds[1] = logical[1];
        mesh->logicalBounds[2] = logical[2];

        md->Add(mesh);
        meshes_added.insert(mesh_for_this_var);
      }

      //Add meshvars
      for (int j=0; j<(int)stepInfo->varInfo[i].materials.size(); j++) {
        char buffer[128];
        string newVarname = varname;
        sprintf(buffer, "%d", stepInfo->varInfo[i].materials[j]);
        newVarname.append("/");
        newVarname.append(buffer);

        if (mesh_vars_added.find(mesh_for_this_var+newVarname)==mesh_vars_added.end()) {
          mesh_vars_added.insert(mesh_for_this_var+newVarname);

          if (vartype.find("Vector") != string::npos)
            AddVectorVarToMetaData(md, newVarname, mesh_for_this_var, cent, 3); // 3 -> vector dimension
          else if (vartype.find("Matrix3") != string::npos)
            AddTensorVarToMetaData(md, newVarname, mesh_for_this_var, cent, 9); // 9 -> tensor 
          else 
            AddScalarVarToMetaData(md, newVarname, mesh_for_this_var, cent);
        }
      }
    }   
  }

  // add a proc id enum variable
  if (addProcId)
  {
    avtScalarMetaData *scalar = new avtScalarMetaData();

    scalar->name = "proc_id";
    scalar->meshName = mesh_for_procid;
    scalar->centering = AVT_ZONECENT;
    scalar->hasDataExtents = false;
    scalar->treatAsASCII = false;
    md->Add(scalar);
  }
  

  // Nothing needs to be modifed for particle data, as they exist only on a single level
  for (int i=0; i<(int)stepInfo->varInfo.size(); i++) {
    if (stepInfo->varInfo[i].type.find("ParticleVariable") != string::npos) {
      string varname = stepInfo->varInfo[i].name;
      string vartype = stepInfo->varInfo[i].type;

      // j=-1 -> all materials (*)
      for (int j=-1; j<(int)stepInfo->varInfo[i].materials.size(); j++) {
        string mesh_for_this_var = string("Particle_Mesh/");
        string newVarname = varname+"/";

        if (j >= 0) {
          char buffer[128];
          sprintf(buffer, "%d", stepInfo->varInfo[i].materials[j]);
          mesh_for_this_var.append(buffer);
          newVarname.append(buffer);
        }
        else {
          mesh_for_this_var.append("*");
          newVarname.append("*");
        }

        if (meshes_added.find(mesh_for_this_var)==meshes_added.end()) {

          avtMeshMetaData *mesh = new avtMeshMetaData;

          mesh->name = mesh_for_this_var;
          mesh->meshType = AVT_POINT_MESH;
          mesh->topologicalDimension = 0;
          mesh->spatialDimension = 3;

          mesh->numBlocks = totalPatches;
          mesh->blockTitle = "patches";
          mesh->blockPieceName = "patch";
          mesh->numGroups = numLevels;
          mesh->groupTitle = "levels";
          mesh->groupPieceName = "level";
          mesh->blockNames = pieceNames;

          mesh->hasSpatialExtents = true; 
          mesh->minSpatialExtents[0] = box_min[0];
          mesh->maxSpatialExtents[0] = box_max[0];
          mesh->minSpatialExtents[1] = box_min[1];
          mesh->maxSpatialExtents[1] = box_max[1];
          mesh->minSpatialExtents[2] = box_min[2];
          mesh->maxSpatialExtents[2] = box_max[2];

          mesh->hasLogicalBounds = true;
          mesh->logicalBounds[0] = logical[0];
          mesh->logicalBounds[1] = logical[1];
          mesh->logicalBounds[2] = logical[2];

          md->Add(mesh); 
          meshes_added.insert(mesh_for_this_var);
        }

        if (mesh_vars_added.find(mesh_for_this_var+newVarname)==mesh_vars_added.end()) {
          mesh_vars_added.insert(mesh_for_this_var+newVarname);

          avtCentering cent = AVT_NODECENT;
          if ((vartype.find("Vector") != string::npos) ||
              (vartype.find("Point") != string::npos))
            AddVectorVarToMetaData(md, newVarname, mesh_for_this_var, cent, 3); // 3 -> vector dimension
          else if (vartype.find("Matrix3") != string::npos)
            AddTensorVarToMetaData(md, newVarname, mesh_for_this_var, cent, 9); // 9 -> tensor 
          else
            AddScalarVarToMetaData(md, newVarname, mesh_for_this_var, cent);
        }
      }
    }   
  }
  
  md->AddGroupInformation(numLevels, totalPatches, groupIds);
  md->AddDefaultSILRestrictionDescription(std::string("!TurnOnAll"));


  AddExpressionsToMetadata(md);
}


// ****************************************************************************
//  Method: udaFileFormat::PopulateDatabaseMetaData
//
//  Purpose:
//      This database meta-data object is like a table of contents for the
//      file.  By populating it, you are telling the rest of VisIt what
//      information it can request from you.
//
//  Programmer: sshankar -- generated by xml2avt
//  Creation:   Tue May 13 19:02:26 PST 2008
//
// ****************************************************************************
void
udaFileFormat::PopulateDatabaseMetaData(avtDatabaseMetaData *md, int timeState)
{
#ifdef SERIALIZED_READS

  int numProcs, rank;
  int msg = 128, tag = 256;
  MPI_Status status;

  MPI_Comm_size(VISIT_MPI_COMM, &numProcs);
  MPI_Comm_rank(VISIT_MPI_COMM, &rank);
  //debug5 << "Proc: " << rank << " sent to mdserver" << endl;  

  if (rank == 0) {
    ReadMetaData(md, timeState);
    MPI_Send(&msg, 1, MPI_INT, 1, tag, VISIT_MPI_COMM);
  }
  else {
    MPI_Recv(&msg, 1, MPI_INT, rank - 1, tag, VISIT_MPI_COMM, &status);
    if (msg == 128 && tag == 256) {
      ReadMetaData(md, timeState);
      if (rank < (numProcs - 1))
        MPI_Send(&msg, 1, MPI_INT, rank + 1, tag, VISIT_MPI_COMM);
    }
  }
#else      
  ReadMetaData(md, timeState);
#endif
}



// ****************************************************************************
//  Method: udaFileFormat::GetGlobalDomainNumber
//
//  Purpose:
//      Translates the level and local patch number into a global patch id.
//  
// ****************************************************************************
int udaFileFormat::GetGlobalDomainNumber(int level, int local_patch) {
  int g=0;
  for (int l=0; l<level; l++)
    g += stepInfo->levelInfo[l].patchInfo.size();
  g += local_patch;

  return g;
}


