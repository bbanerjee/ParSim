#include <Common/DataArchive.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <CCA/Components/ProblemSpecification/ProblemSpecReader.h>
#include <CCA/Ports/InputContext.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/LoadBalancer.h>
#include <Core/Domain/Domain.h>
#include <Core/Domain/UnknownVariable.h>
#include <Core/Domain/Variables/VarLabel.h>
#include <Core/Domain/Level.h>
#include <Core/Math/MiscMath.h>

#include <Core/Exceptions/InternalError.h>
#include <Core/Util/Assert.h>
#include <Core/Thread/Time.h>
#include <Core/Util/DebugStream.h>
#include <Core/Containers/OffsetArray1.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <fcntl.h>

#ifdef _WIN32
#  include <io.h>
#else
#  include <sys/param.h>
#  include <unistd.h>
#endif


using namespace std;
using namespace Matiti;
using namespace SCIRun;

DebugStream DataArchive::dbg("DataArchive", false);

DataArchive::DataArchive(const std::string& filebase,
                         bool verbose /* = true */ ) :
  timestep_cache_size(10), default_cache_size(10), 
  d_filebase(filebase), 
{
  if( d_filebase == "" ) {
    throw InternalError("DataArchive::DataArchive 'filebase' cannot be empty (\"\").", __FILE__, __LINE__);
  }

  while( d_filebase[ d_filebase.length() - 1 ] == '/' ) {
    // Remove '/' from the end of the filebase (if there is one).
    d_filebase = d_filebase.substr( 0, filebase.length() - 1 );
  }

  string index = d_filebase + "/index.xml";
  if( verbose && processor == 0) {
    // cerr << "Parsing " << index << "\n";
  }

  d_indexDoc = ProblemSpecReader().readInputFile( index );

  d_globalEndianness = "";
  d_globalNumBits = -1;
  queryEndiannessAndBits(d_indexDoc, d_globalEndianness, d_globalNumBits);
}


DataArchive::~DataArchive()
{
  //d_indexDoc->releaseDocument();
}

// static, so can be called from either DataArchive or TimeData
void
DataArchive::queryEndiannessAndBits(ProblemSpecP doc, string& endianness, int& numBits)
{
  ProblemSpecP meta = doc->findBlock("Meta");

  if (meta == 0) {
    return;
  }

  ProblemSpecP endian_node = meta->findBlock("endianness");
  if (endian_node) {
    endianness = endian_node->getNodeValue();
  }

  ProblemSpecP nBits_node = meta->findBlock("nBits");
  if( nBits_node) {
    numBits = atoi(nBits_node->getNodeValue().c_str());
  }
}

void
DataArchive::queryTimesteps( std::vector<int>& index,
                             std::vector<double>& times )
{
  double start = Time::currentSeconds();
  if(d_timeData.size() == 0){
    if(d_timeData.size() == 0){
      ProblemSpecP ts = d_indexDoc->findBlock("timesteps");
      if(ts == 0) {
        throw InternalError("DataArchive::queryTimesteps 'timesteps' node not found in index.xml",
                            __FILE__, __LINE__);
      }
      for(ProblemSpecP t = ts->getFirstChild(); t != 0; t = t->getNextSibling()){
        if(t->getNodeType() == ProblemSpec::ELEMENT_NODE){
          map<string,string> attributes;
          t->getAttributes(attributes);
          string tsfile = attributes["href"];
          if(tsfile == "")
            throw InternalError("DataArchive::queryTimesteps:timestep href not found",
                                  __FILE__, __LINE__);
          

          int timestepNumber;
          double currentTime;
          string ts = d_filebase + "/" + tsfile;
          ProblemSpecP timestepDoc = 0;

          if(attributes["time"] == "") {
            // This block if for earlier versions of the index.xml file that do not
            // contain time information as attributes of the timestep field.

            timestepDoc = ProblemSpecReader().readInputFile( ts );
            
            ProblemSpecP time = timestepDoc->findBlock("Time");
            if(time == 0)
              throw InternalError("DataArchive::queryTimesteps:Cannot find Time block",
                                  __FILE__, __LINE__);
            
            if(!time->get("timestepNumber", timestepNumber))
              throw InternalError("DataArchive::queryTimesteps:Cannot find timestepNumber",
                                  __FILE__, __LINE__);
            
            if(!time->get("currentTime", currentTime))
              throw InternalError("DataArchive::queryTimesteps:Cannot find currentTime",
                                  __FILE__, __LINE__);
          } else {
            // This block will read delt and time info from the index.xml file instead of
            // opening every single timestep.xml file to get this information
            istringstream timeVal(attributes["time"]);
            istringstream timestepVal(t->getNodeValue());

            timeVal >> currentTime;
            timestepVal >> timestepNumber;

            if( !timeVal || !timestepVal ) {
              printf( "WARNING: DataArchive.cc: stringstream failed...\n" );
            }

          }

          d_tsindex.push_back(timestepNumber);
          d_tstimes.push_back(currentTime);
          d_timeData.push_back(TimeData(this, timestepDoc, ts));
        }
      }
    }
  }
  index=d_tsindex;
  times=d_tstimes;
  dbg << "DataArchive::queryTimesteps completed in " << Time::currentSeconds()-start << " seconds\n";
}

DataArchive::TimeData& 
DataArchive::getTimeData(int index)
{
  //ASSERTRANGE(index, 0, (int)d_timeData.size());
  TimeData& td = d_timeData[index];
  if (!td.d_initialized)
    td.init();

  list<int>::iterator is_cached = std::find(d_lastNtimesteps.begin(), d_lastNtimesteps.end(), index);
  if (is_cached != d_lastNtimesteps.end()) {
    // It's in the list, so yank it in preperation for putting it at
    // the top of the list.
    dbg << "Already cached, putting at top of list.\n";
    d_lastNtimesteps.erase(is_cached);
  } else {
    dbg << "Not in list.\n";
    // Not in the list.  If the list is maxed out, purge the cache
    // of the last item by removing it from the list.  If
    // timestep_cache_size is <= 0, there is an unlimited size to
    // the cache, so don't purge.
    dbg << "timestep_cache_size = "<<timestep_cache_size<<", d_lastNtimesteps.size() = "<<d_lastNtimesteps.size()<<"\n";
    if (timestep_cache_size > 0 && (int)(d_lastNtimesteps.size()) >= timestep_cache_size) {
      int cacheTimestep = d_lastNtimesteps.back();
      d_lastNtimesteps.pop_back();
      dbg << "Making room.  Purging index "<< cacheTimestep <<"\n";
      d_timeData[cacheTimestep].purgeCache();
    }
  }
  // Finally insert our new candidate at the top of the list.
  d_lastNtimesteps.push_front(index);
  return td;
}

DomainP
DataArchive::queryDomain( int index, const ProblemSpec* ups)
{

  double start = Time::currentSeconds();
  TimeData& timedata = getTimeData(index);

  if (timedata.d_domain != 0) {
    return timedata.d_domain;
  } 

  timedata.d_meshInfo.clear();
  timedata.d_matlInfo.clear();

  const ProblemSpecP top = timedata.d_tstop;
  if (top == 0)
    throw InternalError("DataArchive::queryDomain:Cannot find Domain in timestep",
                        __FILE__, __LINE__);
  ProblemSpecP domainnode = top->findBlock("Domain");
  if(domainnode == 0)
    throw InternalError("DataArchive::queryDomain:Cannot find Domain in timestep",
                        __FILE__, __LINE__);
  DomainP domain = scinew Domain;
  for(ProblemSpecP n = domainnode->getFirstChild(); n != 0; n=n->getNextSibling()){

      timedata.d_meshInfo.push_back(vector<MeshData>());
      timedata.d_matlInfo.push_back(vector<bool>());

      int numMeshes = -1234;
      for(ProblemSpecP r = n->getFirstChild(); r != 0; r=r->getNextSibling()){
        if(r->getNodeName() == "numMeshes") {
          if(!r->get(numMeshes))
            throw InternalError("DataArchive::queryDomain:Error parsing numRegions",
                                __FILE__, __LINE__);
        } else if(r->getNodeName() == "Mesh") {
          int id;
          if(!r->get("id", id))
            throw InternalError("DataArchive::queryDomain:Error parsing mesh id",
                                __FILE__, __LINE__);
          IntVector lowIndex;
          if(!r->get("lowIndex", lowIndex))
            throw InternalError("DataArchive::queryDomain:Error parsing mesh lowIndex",
                                __FILE__, __LINE__);
          IntVector highIndex;
          if(!r->get("highIndex", highIndex))
            throw InternalError("DataArchive::queryDomain:Error parsing mesh highIndex",
                                __FILE__, __LINE__);
          IntVector inLowIndex = lowIndex;
          IntVector inHighIndex = highIndex;
          r->get("interiorLowIndex", inLowIndex);
          r->get("interiorHighIndex", inHighIndex);
          domain->addMesh(lowIndex, highIndex,inLowIndex, inHighIndex,domain.get_rep(),id);
          MeshData pi;
          timedata.d_meshInfo.push_back(pi);
      }
	//ASSERTEQ(level->numMeshes(), numMeshes);
      
      if (ups) {
        // this is not necessary on non-restarts.
        ProblemSpecP domain_ps = ups->findBlock("Domain");
        domain->assignBCS(domain_ps,0);
       }
  }
  
  domain->performConsistencyCheck();

  timedata.d_domain = domain;

  //ASSERTEQ(domain->numLevels(), numLevels);
  dbg << "DataArchive::queryDomain completed in " << Time::currentSeconds()-start << " seconds\n";
  return domain;
}

void
DataArchive::queryVariables( vector<string>& names,
                             vector<const Matiti::TypeDescription*>& types)
{
  double start = Time::currentSeconds();
  ProblemSpecP vars = d_indexDoc->findBlock("variables");
  if(vars == 0)
    throw InternalError("DataArchive::queryVariables:variables section not found\n",
                        __FILE__, __LINE__);
  queryVariables(vars, names, types);

  dbg << "DataArchive::queryVariables completed in " << Time::currentSeconds()-start << " seconds\n";
}

void
DataArchive::queryGlobals( vector<string>& names,
                           vector<const Matiti::TypeDescription*>& types)
{
  double start = Time::currentSeconds();
  ProblemSpecP vars = d_indexDoc->findBlock("globals");
  if(vars == 0)
    return;
  queryVariables(vars, names, types);


  dbg << "DataArchive::queryGlobals completed in " << Time::currentSeconds()-start << " seconds\n";   
}

void
DataArchive::queryVariables(ProblemSpecP vars, vector<string>& names,
                            vector<const Matiti::TypeDescription*>& types)
{
  for(ProblemSpecP n = vars->getFirstChild(); n != 0; n = n->getNextSibling()){
    if(n->getNodeName() == "variable") {
      map<string,string> attributes;
      n->getAttributes(attributes);

      string type = attributes["type"];
      if(type == "")
        throw InternalError("DataArchive::queryVariables:Variable type not found",
                            __FILE__, __LINE__);
      const TypeDescription* td = TypeDescription::lookupType(type);
      if(!td){
        static TypeDescription* unknown_type = 0;
        if(!unknown_type)
          unknown_type = scinew TypeDescription(TypeDescription::Unknown,
                                                "-- unknown type --",
                                                false, MPI_Datatype(-1));
        td = unknown_type;
      }
      types.push_back(td);
      string name = attributes["name"];
      if(name == "")
        throw InternalError("DataArchive::queryVariables:Variable name not found",
                            __FILE__, __LINE__);
      names.push_back(name);
    } else if(n->getNodeType() != ProblemSpec::TEXT_NODE){
      cerr << "DataArchive::queryVariables:WARNING: Unknown variable data: " << n->getNodeName() << '\n';
    }
  }
}

void
DataArchive::query( Variable& var, const std::string& name, int matlIndex, 
                    const Mesh* mesh, int index, DataFileInfo* dfi /* = 0 */)
{
  double tstart = Time::currentSeconds();
  string url;

#if !defined( _WIN32 ) && !defined( DISABLE_SCI_MALLOC )
  const char* tag = AllocatorSetDefaultTag("QUERY");
#endif

  TimeData& timedata = getTimeData(index);

  //ASSERT(timedata.d_initialized);
  // make sure info for this mesh gets parsed from p*****.xml.
  timedata.parseMesh(mesh);

  VarData& varinfo = timedata.d_varInfo[name];
  string dataurl;
  int meshid;
  if (mesh) {
    // we need to use the real_mesh (in case of periodic boundaries) to get the data, but we need the
    // passed in mesh to allocate the mesh to the proper virtual region... (see var.allocate below)
    const Mesh* real_mesh = mesh->getRealMesh();
    MeshData& meshinfo = timedata.d_meshInfo[real_mesh->getLevel()->getIndex()][real_mesh->getLevelIndex()];
    //ASSERT(meshinfo.parsed);
    meshid = real_mesh->getID();

    ostringstream ostr;
    // append l#/datafilename to the directory
    ostr << timedata.d_tsurldir << "l" << mesh->getLevel()->getIndex() << "/" << meshinfo.datafilename;
    dataurl = ostr.str();
  }
  else {
    // reference reduction file 'global.data' will a null mesh
    meshid = -1;
    dataurl = timedata.d_tsurldir + timedata.d_globaldata;
  }

  // on a call from restartInitialize, we already have the information from the dfi,
  // otherwise get it from the hash table info
  DataFileInfo datafileinfo;
  if (!dfi) {
    // if this is a virtual mesh, grab the real mesh, but only do that here - in the next query, we want
    // the data to be returned in the virtual coordinate space
    if (!timedata.d_datafileInfo.lookup(VarnameMatlMesh(name, matlIndex, meshid), datafileinfo)) {
      cerr << "VARIABLE NOT FOUND: " << name << ", material index " << matlIndex << ", mesh " << mesh->getID() << ", time index " << index << "\nPlease make sure the correct material index is specified\n";
      throw InternalError("DataArchive::query:Variable not found",
                          __FILE__, __LINE__);
    }
    dfi = &datafileinfo;
  }
  const TypeDescription* td = var.virtualGetTypeDescription();
  //ASSERT(td->getName() == varinfo.type);
  
  if (td->getType() == TypeDescription::MeshNodeVariable) {
    if(dfi->numMeshNodes == -1)
      throw InternalError("DataArchive::query:Cannot get numMeshNodes",
                          __FILE__, __LINE__);

    if (mesh->isVirtual())
      throw InternalError("DataArchive::query: MeshNode query on virtual meshes "
                          "not finished.  We need to adjust the particle positions to virtual space...", __FILE__, __LINE__);
    psetDBType::key_type key(matlIndex, mesh);
    MeshNodeSubset* psubset = 0;
    psetDBType::iterator psetIter = d_psetDB.find(key);
    if(psetIter != d_psetDB.end()) {
      psubset = (*psetIter).second.get_rep();
    }
    if (psubset == 0 || psubset->numMeshNodes() != dfi->numMeshNodes)
    {
      psubset = scinew MeshNodeSubset(dfi->numMeshNodes, matlIndex, mesh);
      //      cout << "numMeshNodes: " << dfi->numMeshNodes << "\n";
      //      cout << "d_pset size: " << d_psetDB.size() << "\n";
      //      cout << "1. key is: " << key.first << "\n";
      //      cout << "2. key is: " << key.second << "\n";
      d_psetDB[key] = psubset;
    }
    (static_cast<MeshNodeVariableBase*>(&var))->allocate(psubset);
//      (dynamic_cast<MeshNodeVariableBase*>(&var))->allocate(psubset);
  }
  else if (td->getType() != TypeDescription::ReductionVariable) {
    var.allocate(mesh, varinfo.boundaryLayer);
  }
  
#ifdef _WIN32
  int fd = open(dataurl.c_str(), O_RDONLY|O_BINARY);
#else
  int fd = open(dataurl.c_str(), O_RDONLY);
#endif
  if(fd == -1) {
    cerr << "Error opening file: " << dataurl.c_str() << ", errno=" << errno << '\n';
    throw ErrnoException("DataArchive::query (open call)", errno, __FILE__, __LINE__);
  }
  off_t ls = lseek(fd, dfi->start, SEEK_SET);

  if(ls == -1) {
    cerr << "Error lseek - file: " << dataurl.c_str() << ", errno=" << errno << '\n';
    throw ErrnoException("DataArchive::query (lseek call)", errno, __FILE__, __LINE__);
  }
  InputContext ic(fd, dataurl.c_str(), dfi->start);
  double starttime = Time::currentSeconds();
  var.read(ic, dfi->end, timedata.d_swapBytes, timedata.d_nBytes, varinfo.compression);

  dbg << "DataArchive::query: time to read raw data: "<<Time::currentSeconds() - starttime<<endl;
  ASSERTEQ(dfi->end, ic.cur);
  int s = close(fd);
  if(s == -1) {
    cerr << "Error closing file: " << dataurl.c_str() << ", errno=" << errno << '\n';
    throw ErrnoException("DataArchive::query (close call)", errno, __FILE__, __LINE__);
  }

#if !defined( _WIN32 ) && !defined( DISABLE_SCI_MALLOC )
  AllocatorSetDefaultTag(tag);
#endif
  dbg << "DataArchive::query() completed in "
      << Time::currentSeconds()-tstart << " seconds\n";
}

void DataArchive::query( Variable& var, const string& name, int matlIndex, 
                         const Mesh* mesh, int timeIndex,
                         Ghost::GhostType gt, int ngc)
{
  if (ngc == 0)
    query(var, name, matlIndex, mesh, timeIndex, 0);
  else {
    TimeData& td = getTimeData(timeIndex);
    td.parseMesh(mesh); // make sure vars is actually populated
    if (td.d_varInfo.find(name) != td.d_varInfo.end()) {
      VarData& varinfo = td.d_varInfo[name];
      const TypeDescription* type = TypeDescription::lookupType(varinfo.type);
      IntVector low, high;
      mesh->computeVariableExtents(type->getType(), varinfo.boundaryLayer, gt, ngc, low, high);
      queryRegion(var, name, matlIndex, mesh->getLevel(), timeIndex, low, high);
    }
    else {
      cerr << "VARIABLE NOT FOUND: " << name << ", material index " << matlIndex << ", mesh " << mesh->getID() << ", time index " 
           << timeIndex << "\nPlease make sure the correct material index is specified\n";
      throw InternalError("DataArchive::query:Variable not found",
                          __FILE__, __LINE__);
    }
  }
}

void DataArchive::queryRegion(Variable& var, const string& name, int matlIndex, 
                              const Level* level, int timeIndex, IntVector low, IntVector high)
{
  // NOTE - this is not going to do error checking like making sure the entire volume is filled.  
  //   We'll assume that if there were bad regions, they would have been caught in the simulation.
  DomainVariableBase* domainvar = dynamic_cast<DomainVariableBase*>(&var);
  ASSERT(domainvar);
  domainvar->allocate(low, high);

  TimeData& td = getTimeData(timeIndex);
  const TypeDescription* type = 0;
  Mesh::VariableBasis basis = Mesh::NodeBased; // not sure if this is a reasonable default...
  Mesh::selectType meshes;
  
  level->selectMeshes(low, high, meshes);
  for(int i=0;i<meshes.size();i++){
    const Mesh* mesh = meshes[i];
    
    if (type == 0) {
      td.parseMesh(mesh); // make sure varInfo is loaded
      VarData& varinfo = td.d_varInfo[name];
      type = TypeDescription::lookupType(varinfo.type);
      basis = Mesh::translateTypeToBasis(type->getType(), false);
    }
    IntVector l, h;

    l = Max(mesh->getExtraLowIndex(basis, IntVector(0, 0, 0)), low);
    h = Min(mesh->getExtraHighIndex(basis, IntVector(0, 0, 0)), high);
    if (l.x() >= h.x() || l.y() >= h.y() || l.z() >= h.z())
      continue;
    DomainVariableBase* tmpVar = domainvar->cloneType();
    query(*tmpVar, name, matlIndex, mesh, timeIndex);

    if (mesh->isVirtual()) {
      // if mesh is virtual, it is probable a boundary layer/extra cell that has been requested (from AMR)
      // let Bryan know if this doesn't work.  We need to adjust the source but not the dest by the virtual offset
      tmpVar->offset(mesh->getVirtualOffset());
    }
    try {
      domainvar->copyMesh(tmpVar, l, h);
    } catch (InternalError& e) {
      cout << " Bad range: " << low << " " << high << ", mesh intersection: " << l << " " << h 
           << " actual mesh " << mesh->getLowIndex(basis) << " " << mesh->getHighIndex(basis) 
           << " var range: "  << tmpVar->getLow() << " " << tmpVar->getHigh() << endl;
      throw e;
    }
    delete tmpVar;
  }
}



void 
DataArchive::findMeshAndIndex(DomainP domain, Mesh*& mesh, particleIndex& idx,
                               long64 particleID, int matlIndex, int levelIndex,
                               int index)
{
  Mesh *local = mesh;
  if( mesh != NULL ){
    MeshNodeVariable<long64> var;
    query(var, "p.particleID", matlIndex, mesh, index);
    //  cerr<<"var["<<idx<<"] = "<<var[idx]<<endl;
    if( idx < var.getMeshNodeSubset()->numMeshNodes() && var[idx] == particleID )
      return;
    else {
      MeshNodeSubset* subset = var.getMeshNodeSubset();
      for(MeshNodeSubset::iterator p_iter = subset->begin();
          p_iter != subset->end(); p_iter++){
        if( var[*p_iter] == particleID){
          idx = *p_iter;
          return;
        }
      }
    }
  }
  mesh = NULL;
//   for (int level_nr = 0;
//        (level_nr < domain->numLevels()) && (mesh == NULL); level_nr++) {
    
//     const LevelP level = domain->getLevel(level_nr);
    const LevelP level = domain->getLevel(levelIndex);
    
    for (Level::const_meshIterator iter = level->meshesBegin();
         (iter != level->meshesEnd()) && (mesh == NULL); iter++) {
      if( *iter == local ) continue;
      MeshNodeVariable<long64> var;
      query(var, "p.particleID", matlIndex, *iter, index);
      MeshNodeSubset* subset = var.getMeshNodeSubset();
      for(MeshNodeSubset::iterator p_iter = subset->begin();
          p_iter != subset->end(); p_iter++){
        if( var[*p_iter] == particleID){
          mesh = *iter;
          idx = *p_iter;
          //      cerr<<"var["<<*p_iter<<"] = "<<var[*p_iter]<<endl;
          break;
        }
      }
      
      if( mesh != NULL )
        break;
    }
//  }
}

void
DataArchive::restartInitialize(int index, const DomainP& domain, DataWarehouse* dw,
                               LoadBalancer* lb, double* pTime)
{
  vector<int> indices;
  vector<double> times;
  queryTimesteps(indices, times);

  vector<string> names;
  vector< const TypeDescription *> typeDescriptions;
  queryVariables(names, typeDescriptions);
  queryGlobals(names, typeDescriptions);  
  
  map<string, VarLabel*> varMap;
  for (unsigned i = 0; i < names.size(); i++) {
    VarLabel * vl = VarLabel::find(names[i]);
    if( vl == NULL ) {
//      proc0cout << "Warning, VarLabel for " << names[i] << " was not found... attempting to create.\n"
//          << "However, it is possible that this may cause problems down the road...\n";
      //***** THIS ASSUMES A SINGLE GHOST CELL ***** BE CAREFUL ********
      // check if we have extracells specified. This affects Wasatch only and should have no impact on other components.
      const bool hasExtraCells = (domain->getMeshByID(0,0)->getExtraCells() != SCIRun::IntVector(0,0,0));
      // if extracells are specified, then create varlabels that are consistent with Wasatch varlabels.
      vl = VarLabel::create( names[i], typeDescriptions[i], hasExtraCells? IntVector(0,0,0) : IntVector(1,1,1) );
    }
    varMap[names[i]] = vl;
  }

  TimeData& timedata = getTimeData(index);

  *pTime = times[index];

  if (lb)
    lb->restartInitialize(this, index, timedata.d_tstop, timedata.d_tsurl, domain);

  // set here instead of the SimCont because we need the DW ID to be set 
  // before saving particle subsets
  dw->setID( indices[index]);
  
  // make sure to load all the data so we can iterate through it 
  for (int l = 0; l < domain->numLevels(); l++) {
    LevelP level = domain->getLevel(l);
    for (int p = 0; p < level->numMeshes(); p++) {
      const Mesh* mesh = level->getMesh(p);
      if (lb->getMeshwiseProcessorAssignment(mesh) == d_processor)
        timedata.parseMesh(mesh);
    }
  }

  // iterate through all entries in the VarData hash table, and loading the 
  // variables if that data belongs on this processor
  VarHashMapIterator iter(&timedata.d_datafileInfo);
  iter.first();
  for (; iter.ok(); ++iter) {
    VarnameMatlMesh& key = iter.get_key();
    DataFileInfo& data = iter.get_data();

    // get the Mesh from the Mesh ID (ID of -1 = NULL - for reduction vars)
    const Mesh* mesh = key.meshid_ == -1 ? NULL : domain->getMeshByID(key.meshid_, 0);
    int matl = key.matlIndex_;

    VarLabel* label = varMap[key.name_];
    if (label == 0) {
      throw UnknownVariable(key.name_, dw->getID(), mesh, matl,
                            "on DataArchive::scheduleRestartInitialize",
                            __FILE__, __LINE__);
    }

    if (!mesh || !lb || lb->getMeshwiseProcessorAssignment(mesh) == d_processor) {
      Variable* var = label->typeDescription()->createInstance();
      query(*var, key.name_, matl, mesh, index, &data);

      MeshNodeVariableBase* particles;
      if ((particles = dynamic_cast<MeshNodeVariableBase*>(var))) {
        if (!dw->haveMeshNodeSubset(matl, mesh)) {
          dw->saveMeshNodeSubset(particles->getMeshNodeSubset(), matl, mesh);
        }
        else {
          ASSERTEQ(dw->getMeshNodeSubset(matl, mesh), particles->getMeshNodeSubset());
        }
      }
      dw->put(var, label, matl, mesh); 
      delete var; // should have been cloned when it was put
    }
  }
}

bool
DataArchive::queryRestartTimestep(int& timestep)
{
  ProblemSpecP restartNode = d_indexDoc->findBlock("restart");
  if (restartNode == 0) {
    ProblemSpecP restartsNode = d_indexDoc->findBlock("restarts");
    if (restartsNode == 0)
      return false;
    
    restartNode = restartsNode->findBlock("restart");
    if (restartNode == 0)
      return false;

    // get the last restart tag in the restarts list
    while (restartNode->findNextBlock("restart") != 0)
      restartNode = restartNode->findNextBlock("restart");
  }
  
  map<string,string> attributes;
  restartNode->getAttributes(attributes);
  string ts = attributes["timestep"];
  if (ts == "")
    return false;
  timestep = atoi(ts.c_str());
  return true;
}

// We want to cache at least a single timestep, so that we don't have
// to reread the timestep for every mesh queried.  This sets the
// cache size to one, so that this condition is held.
void
DataArchive::turnOffXMLCaching() {
  setTimestepCacheSize(1);
}

// Sets the number of timesteps to cache back to the default_cache_size
void
DataArchive::turnOnXMLCaching() {
  setTimestepCacheSize(default_cache_size);
}

// Sets the timestep cache size to whatever you want.  This is useful
// if you want to override the default cache size determined by
// TimeHashMaps.
void
DataArchive::setTimestepCacheSize(int new_size) {
  // Now we need to reduce the size
  int current_size = (int)d_lastNtimesteps.size();
  dbg << "current_size = "<<current_size<<"\n";
  if (timestep_cache_size >= current_size) {
    // everything's fine
    return;
  }

  int kill_count = current_size - timestep_cache_size;
  dbg << "kill_count = "<<kill_count<<"\n";
  for(int i = 0; i < kill_count; i++) {
    int cacheTimestep = d_lastNtimesteps.back();
    dbg << "Making room.  Purging time index "<< cacheTimestep <<"\n";

    d_lastNtimesteps.pop_back();
    d_timeData[cacheTimestep].purgeCache();
  }
}

DataArchive::TimeData::TimeData(DataArchive* da, ProblemSpecP timestepDoc, string timestepURL) :
  d_initialized(false), d_tstop(timestepDoc), d_tsurl(timestepURL), da(da)
{
  d_tsurldir = timestepURL.substr(0, timestepURL.find_last_of('/')+1);
}

DataArchive::TimeData::~TimeData()
{
  purgeCache();
}

void
DataArchive::TimeData::init()
{
  d_initialized=true;
  //  cerr << "MeshHashMaps["<<time<<"]::init\n";
  // grab the data xml files from the timestep xml file
  if (d_tstop == 0) {
    d_tstop = ProblemSpecReader().readInputFile( d_tsurl );
  }

  // Handle endianness and number of bits
  string endianness = da->d_globalEndianness;
  int    numbits    = da->d_globalNumBits;
  DataArchive::queryEndiannessAndBits(d_tstop, endianness, numbits);

  static bool endian_warned = false;
  static bool bits_warned = false;

  if (endianness == "") {
    endianness = string(SCIRun::endianness());
    if (!endian_warned) {
      endian_warned = true;
      cout<<"\nXML Warning: endianness node not found.\n"<<
        "Assuming data was created on a " << SCIRun::endianness() << " machine.\n"<<
        "To eliminate this message and express the correct\n"<<
        "endianess, please add either\n"<<
        "\t<endianness>little_endian</endianness>\n"<<
        "or\n\t<endianness>big_endian</endianness>\n"<<
        "to the <Meta> section of the index.xml file.\n\n";
    }
  }
  if (numbits == -1) {
    numbits = sizeof(unsigned long) * 8;
    if (!bits_warned) {
      cout<<"\nXML Warning: nBits node not found.\n"<<
        "Assuming data was created using " << sizeof(unsigned long) * 8 << " bits.\n"
        "To eliminate this message and express the correct\n"<<
        "number of bits, please add either\n"<<
        "\t<nBits>32</nBits>\n"<<
        "or\n\t<nBits>64</nBits>\n"<<
        "to the <Meta> section of the index.xml file.\n\n";
    }
  }

  d_swapBytes = endianness != string(SCIRun::endianness());
  d_nBytes = numbits / 8;

  ProblemSpecP datanode = d_tstop->findBlock("Data");
  if(datanode == 0)
    throw InternalError("Cannot find Data in timestep", __FILE__, __LINE__);
  for(ProblemSpecP n = datanode->getFirstChild(); n != 0; n=n->getNextSibling()){
    if(n->getNodeName() == "Datafile") {
      map<string,string> attributes;
      n->getAttributes(attributes);
      string proc = attributes["proc"];
      /* - Remove this check for restarts.  We need to accurately
         determine which mesh goes on which proc, and for the moment
         we need to be able to parse all pxxxx.xml files.  --BJW  
      if (proc != "") {
        int procnum = atoi(proc.c_str());
        if ((procnum % numProcessors) != processor)
          continue;
      }
      */
      string datafile = attributes["href"];
      if(datafile == "")
        throw InternalError("timestep href not found", __FILE__, __LINE__);

      if (datafile == "global.xml") {
        parseFile(d_tsurldir + datafile, -1, -1);
      }
      else {

        // get level info out of the xml file: should be lX/pxxxxx.xml
        unsigned level = 0;
        string::size_type start = datafile.find_first_of("l",0, datafile.length()-3);
        string::size_type end = datafile.find_first_of("/");
        if (start != string::npos && end != string::npos && end > start && end-start <= 2)
          level = atoi(datafile.substr(start+1, end-start).c_str());

        if (level >= d_xmlUrls.size()) {
          d_xmlUrls.resize(level+1);
          d_xmlParsed.resize(level+1);
        }

        string url = d_tsurldir + datafile;
        d_xmlUrls[level].push_back(url);
        d_xmlParsed[level].push_back(false);
      }
    }
    else if(n->getNodeType() != ProblemSpec::TEXT_NODE){
      cerr << "WARNING: Unknown element in Data section: " << n->getNodeName() << '\n';
    }
  }
}

void
DataArchive::TimeData::purgeCache()
{
  d_domain = 0;
  d_tstop = 0;

  d_datafileInfo.remove_all();
  d_meshInfo.clear(); 
  d_varInfo.clear();
  d_xmlUrls.clear();
  d_xmlParsed.clear();
  d_initialized = false;
}

// This is the function that parses the p*****.xml file for a single processor.
void
DataArchive::TimeData::parseFile(string urlIt, int levelNum, int baseMesh)
{
  // parse the file
  ProblemSpecP top = ProblemSpecReader().readInputFile( urlIt );
  
  // materials are the same for all meshes on a level - don't parse them for more than one file
  bool addMaterials = levelNum >= 0 && d_matlInfo[levelNum].size() == 0;

  for(ProblemSpecP vnode = top->getFirstChild(); vnode != 0; vnode=vnode->getNextSibling()){
    if(vnode->getNodeName() == "Variable") {
      string varname;
      if(!vnode->get("variable", varname))
        throw InternalError("Cannot get variable name", __FILE__, __LINE__);
      
      int meshid;
      if(!vnode->get("mesh", meshid) && !vnode->get("region", meshid))
        throw InternalError("Cannot get mesh id", __FILE__, __LINE__);
      
      int index;
      if(!vnode->get("index", index))
        throw InternalError("Cannot get index", __FILE__, __LINE__);
      
      if (addMaterials) {
        // set the material to existing.  index+1 to use matl -1
        if (index+1 >= (int)d_matlInfo[levelNum].size())
          d_matlInfo[levelNum].resize(index+2);
        d_matlInfo[levelNum][index] = true;
      }

      map<string,string> attributes;
      vnode->getAttributes(attributes);

      string type = attributes["type"];
      if(type == "")
        throw InternalError("DataArchive::query:Variable doesn't have a type",
                            __FILE__, __LINE__);
      long start;
      if(!vnode->get("start", start))
        throw InternalError("DataArchive::query:Cannot get start", __FILE__, __LINE__);
      long end;
      if(!vnode->get("end", end))
        throw InternalError("DataArchive::query:Cannot get end",
                            __FILE__, __LINE__);
      string filename;  
      if(!vnode->get("filename", filename))
        throw InternalError("DataArchive::query:Cannot get filename",
                            __FILE__, __LINE__);

      // not required
      string compressionMode = "";  
      IntVector boundary(0,0,0);
      int numMeshNodes = -1;

      vnode->get("compression", compressionMode);      
      vnode->get("boundaryLayer", boundary);
      vnode->get("numMeshNodes", numMeshNodes);

      if (d_varInfo.find(varname) == d_varInfo.end()) {
        VarData& varinfo = d_varInfo[varname];
        varinfo.type = type;
        varinfo.compression = compressionMode;
        varinfo.boundaryLayer = boundary;
      }
      else if (compressionMode != "") {
        // For particles variables of size 0, the uda doesn't say it
        // has a compressionMode...  (FYI, why is this?  Because it is
        // ambiguous... if there is no data, is it compressed?)
        //
        // To the best of my understanding, we only look at the variables stats
        // the first time we encounter it... even if there are multiple materials.
        // So we run into a problem is the variable has 0 data the first time it
        // is looked at... The problem there is that it doesn't mark it as being
        // compressed, and therefore the next time we see that variable (eg, in
        // another material) we (used to) assume it was not compressed... the 
        // following lines compenstate for this problem:
        VarData& varinfo = d_varInfo[varname];
        varinfo.compression = compressionMode;
      }

      if (levelNum == -1) { // global file (reduction vars)
        d_globaldata = filename;
      }
      else {
        ASSERTRANGE(meshid-baseMesh, 0, (int)d_meshInfo[levelNum].size());

        MeshData& meshinfo = d_meshInfo[levelNum][meshid-baseMesh];
        if (!meshinfo.parsed) {
          meshinfo.parsed = true;
          meshinfo.datafilename = filename;
        }
      }
      VarnameMatlMesh vmp(varname, index, meshid);
      DataFileInfo dummy;

      if (d_datafileInfo.lookup(vmp, dummy) == 1) {
        //cerr << "Duplicate variable name: " << name << endl;
      }
      else {
        DataFileInfo dfi(start, end, numMeshNodes);
        d_datafileInfo.insert(vmp, dfi);
      }
    } else if(vnode->getNodeType() != ProblemSpec::TEXT_NODE){
      cerr << "WARNING: Unknown element in Variables section: " << vnode->getNodeName() << '\n';
    }
  }
  //top->releaseDocument();
}


void
DataArchive::TimeData::parseMesh(const Mesh* mesh)
{
  //ASSERT(d_domain != 0);
  if (!mesh) return;

  const Mesh* real_mesh = mesh->getMesh();
  
  // make sure the data for this mesh has been processed.
  // Return straightaway if we have parsed this mesh
  int meshIndex = real_mesh->getIndex();

  MeshData& meshinfo = d_meshInfo[meshIndex];
  if (meshinfo.parsed)
    return;

  //If this is a newer uda, the mesh info in the domain will store the processor where the data is
  if (meshinfo.proc != -1) {
    ostringstream file;
    file << d_tsurldir << "l" << (int) real_mesh->getIndex() << "/p" << setw(5) << setfill('0') << (int) meshinfo.proc << ".xml";
    parseFile(file.str());
  }
  // Try making a guess as to the processor.  First go is to try
  // the processor of the same index as the mesh.  Many datasets
  // have only one mesh per processor, so this is a reasonable
  // first attempt.  Future attemps could perhaps be smarter.
  if (!meshinfo.parsed && meshIndex < (int)d_xmlParsed[levelIndex].size() && !d_xmlParsed[levelIndex][meshIndex]) {
    parseFile(d_xmlUrls[levelIndex][meshIndex], levelIndex, levelBaseMeshID);
    d_xmlParsed[levelIndex][meshIndex] = true;
  }

  // failed the guess - parse the entire dataset for this level
  if (!meshinfo.parsed) {
    for (unsigned proc = 0; proc < d_xmlUrls[levelIndex].size(); proc++) {
      parseFile(d_xmlUrls[levelIndex][proc], levelIndex, levelBaseMeshID);
      d_xmlParsed[levelIndex][proc] = true;
    }
  }

}


ConsecutiveRangeSet
DataArchive::queryMaterials( const string& varname,
                             const Mesh* mesh,
                             int index )
{
  double start = Time::currentSeconds();

  TimeData& timedata = getTimeData(index);
  timedata.parseMesh(mesh);

  ConsecutiveRangeSet matls;

  for (unsigned i = 0; i < timedata.d_matlInfo[mesh->getIndex()].size(); i++) {
    // i-1, since the matlInfo is adjusted to allow -1 as entries
    VarnameMatlMesh vmp(varname, i-1, mesh->getMesh()->getID());
    DataFileInfo dummy;

    if (timedata.d_datafileInfo.lookup(vmp, dummy) == 1)
      matls.addInOrder(i-1);

  }

  dbg << "DataArchive::queryMaterials completed in " << Time::currentSeconds()-start << " seconds\n";

  return matls;
}

int
DataArchive::queryNumMaterials(const Mesh* mesh, int index)
{
  double start = Time::currentSeconds();

  TimeData& timedata = getTimeData(index);

  timedata.parseMesh(mesh);

  int numMatls = -1;

  for (unsigned i = 0; i < timedata.d_matlInfo[mesh->getIndex()].size(); i++) {
    if (timedata.d_matlInfo[mesh->getIndex()][i]) {
      numMatls++;
    }
  }

  dbg << "DataArchive::queryNumMaterials completed in " << Time::currentSeconds()-start << " seconds\n";

  return numMatls;
}

