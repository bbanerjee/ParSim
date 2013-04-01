#ifndef MATITI_DATAARCHIVE_H
#define MATITI_DATAARCHIVE_H

#include <Mesh/MeshNodeVariable.h>
#include <Mesh/MeshElementVariable.h>
#include <Core/Domain/Domain.h>
#include <Core/Domain/DomainP.h>
#include <Core/Domain/Variables/VarnameMatlMesh.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Common/Handle.h>
#include <Core/Exceptions/VariableNotFoundInDomain.h>
#include <Core/Thread/Mutex.h>
#include <Core/Thread/Time.h>
#include <Core/Util/DebugStream.h>
#include <Core/Containers/ConsecutiveRangeSet.h>
#include <Core/Containers/HashTable.h>

#include   <string>
#include   <vector>
#include   <list>

#include <fcntl.h>

#ifndef _WIN32
#  include <unistd.h>
#endif


namespace Matiti {

  using namespace SCIRun;

  class VarLabel;
  class DataWarehouse;


  //! Container to hold data when read in from disk.
  class DataArchive {
    private:

      // what we need to store on a per-variable basis
      // everything else can be retrieved from a higher level
      struct DataFileInfo {
        DataFileInfo(long s, long e, long np) : start(s), end(e), numMeshNodes(np) {}
        DataFileInfo() {}
        long start;
        long end;
        int numMeshNodes;
      };

      // store these in separate arrays so we don't have to store nearly as many of them
      struct VarData {
        std::string type;
        std::string compression;
        IntVector boundaryLayer;
      };

      struct MeshData {
        MeshData() : parsed(false) {}
        bool parsed;
        std::string datafilename;
      };

      typedef HashTable<VarnameMatlMesh, DataFileInfo> VarHashMap;
      typedef HashTableIter<VarnameMatlMesh, DataFileInfo> VarHashMapIterator;

      //! Top of DataArchive structure for storing hash maps of variable data
      //! - containing data for each time step.
      class TimeData {
        public:    
          TimeData(DataArchive* da, ProblemSpecP timestepDoc, std::string timestepURL);
          ~TimeData();
          VarData& findVariableInfo(const std::string& name, const Mesh* mesh, int matl);

          // reads timestep.xml and prepares the data xml files to be read
          void init();
          void purgeCache(); // purge the cached data

          // makes sure that (if the mesh data exists) then it is parsed.  Try logic to pick
          // the right file first, and if you can't, parse everything
          void parseMesh(const Mesh* mesh);

          // parse an individual data file and load appropriate storage
          void parseFile(std::string file, int baseMesh);

          // This would be private data, except we want DataArchive to have access,
          // so we would mark DataArchive as 'friend', but we're already a private
          // nested class of DataArchive...

          // info in the data file about the mesh-matl-var
          VarHashMap d_datafileInfo;

          // Mesh info (separate by levels) - proc, whether parsed, datafile, etc.
          // Gets expanded and proc is set during queryDomain.  Other fields are set
          // when parsed
          // Organized in a contiguous array, by mesh-level-index
          std::vector<std::vector<MeshData> > d_meshInfo; 

          // Wheter a material is active per level
          std::vector<std::vector<bool> > d_matlInfo;

          // var info - type, compression, and boundary layer
          std::map<std::string, VarData> d_varInfo; 

          // xml urls referred to in timestep.xml
          std::vector<std::vector<std::string> > d_xmlUrls;
          std::vector<std::vector<bool> > d_xmlParsed;

          std::string d_globaldata;

          ConsecutiveRangeSet matls;  // materials available this timestep

          DomainP d_domain;               // store the domain...
          bool d_initialized;         // Flagged once this mesh's init is called
          ProblemSpecP d_tstop;       // ProblemSpecP of timestep.xml
          std::string d_tsurl;        // path to timestep.xml
          std::string d_tsurldir;     // dir that contains timestep.xml
          bool d_swapBytes;
          int d_nBytes;
          DataArchive* da;            // pointer for parent DA.  Need for backward-compatibility with endianness, etc.
      };

    public:
      DataArchive(const std::string& filebase,
          bool verbose = true ); // If you want error messages printed to the screen.

      // GROUP: Destructors
      //////////
      // Destructor
      virtual ~DataArchive();

      TimeData& getTimeData(int index);


      std::string name(){ return d_filebase;}

      //! Set up data arachive for restarting a Matiti simulation   
      void restartInitialize(int timestep, const DomainP& domain, DataWarehouse* dw,
          double* pTime /* passed back */);

      inline ProblemSpecP getTimestepDoc(int index) { return getTimeData(index).d_tstop; }

      static void queryEndiannessAndBits(ProblemSpecP, std::string& endianness, int& numBits);  

      // GROUP:  Information Access
      //////////
      // However, we need a means of determining the names of existing
      // variables. We also need to determine the type of each variable.
      // Get a list of scalar or vector variable names and  
      // a list of corresponding data types
      void queryVariables( std::vector< std::string>& names,
          std::vector< const TypeDescription *>&  );
      void queryGlobals( std::vector< std::string>& names,
          std::vector< const TypeDescription *>&  );
      void queryTimesteps( std::vector<int>& index,
          std::vector<double>& times );

      //! the ups is for the assignBCS that needs to happen
      //! if we are reading the simulation domain from the uda,
      //! and thus is only necessary on a true restart.
      DomainP queryDomain(int index, const ProblemSpec* ups = 0);




#if 0
      //////////
      // Does a variable exist in a particular mesh?
      bool exists(const std::string&, const Mesh*, int) {
        return true;
      }
#endif

      ConsecutiveRangeSet queryMaterials(const std::string& varname,
          const Mesh* mesh, int index);

      int queryNumMaterials(const Mesh* mesh, int index);

      // Queries a variable for a material, mesh, and index in time.
      // Optionally pass in DataFileInfo if you're iterating over
      // entries in the hash table (like restartInitialize does)
      void query( Variable& var, const std::string& name, int matlIndex, 
          const Mesh* mesh, int timeIndex, DataFileInfo* dfi = 0);

      //////////
      // query the variable value for a particular meshNode  overtime;
      template<class T>
        void query(std::vector<T>& values, const std::string& name,
            int matlIndex, long64 meshNodeID, 
            double startTime, double endTime) ;
      //////////
      // similarly, we want to be able to track variable values in a particular
      // mesh element over time.
      template<class T>
        void query(std::vector<T>& values, const std::string& name, int matlIndex,
            IntVector loc, double startTime, double endTime);

      //////////
      // Pass back the timestep number specified in the "restart" tag of the
      // index file, or return false if such a tag does not exist.
      bool queryRestartTimestep(int& timestep);
#if 0
      //////////
      // In other cases we will have noticed something interesting and we
      // will want to access some small portion of a mesh.  We will need
      // to request some range of data in index space.
      template<class T> void get(T& data, const std::string& name,
          const Mesh* mesh, elementIndex min, elementIndex max);
#endif

      // Only cache a single timestep
      void turnOnXMLCaching();
      // Cache the default number of timesteps
      void turnOffXMLCaching();
      // Cache new_size number of timesteps.  Calls the
      // TimeHashMaps::updateCacheSize function with new_size.  See
      // corresponding documentation.
      void setTimestepCacheSize(int new_size);

      // This is a list of the last n timesteps accessed.  Data from
      // only the last timestep_cache_size timesteps is stored, unless
      // timestep_cache_size is less than or equal to zero then the size
      // is unbounded.
      std::list<int> d_lastNtimesteps;

      // Tells you the number of timesteps to cache. Less than or equal to
      // zero means to cache all of them.
      int timestep_cache_size;

      // This will be the default number of timesteps cached, determined
      // by the number of processors.
      int default_cache_size;

    protected:
      DataArchive();

    private:
      friend class DataArchive::TimeData;
      DataArchive(const DataArchive&);
      DataArchive& operator=(const DataArchive&);

      void queryVariables( const ProblemSpecP vars, std::vector<std::string>& names,
          std::vector<const TypeDescription*>& types);

      std::string d_filebase;  
      ProblemSpecP d_indexDoc;
      ProblemSpecP d_restartTimestepDoc;
      std::string d_restartTimestepURL;

      bool d_simRestart;

      std::vector<TimeData> d_timeData;
      std::vector<int> d_tsindex;
      std::vector<double> d_tstimes;

      // global bits and endianness - read from index.xml ONLY if not in timestep.xml
      std::string d_globalEndianness;
      int d_globalNumBits;

      typedef std::map<std::pair<int, const Mesh*>, Handle<MeshNodeSubset> > psetDBType;
      psetDBType d_psetDB;

      void findMeshAndIndex(DomainP domain, Mesh*& mesh, meshNodeIndex& idx,
          long64 meshNodeID, int matIndex, int index);

      static DebugStream dbg;
  };


  template<class T>
    void DataArchive::query(std::vector<T>& values, const std::string& name,
        int matlIndex, long64 meshNodeID,
        double startTime, double endTime)
    {
      double call_start = SCIRun::Time::currentSeconds();

      std::vector<int> index;
      std::vector<double> times;
      queryTimesteps(index, times); // build timesteps if not already done

      // figure out what kind of variable we're looking for
      std::vector<std::string> type_names;
      std::vector<const TypeDescription*> type_descriptions;
      queryVariables(type_names, type_descriptions);
      const TypeDescription* type = NULL;
      std::vector<std::string>::iterator name_iter = type_names.begin();
      std::vector<const TypeDescription*>::iterator type_iter = type_descriptions.begin();
      for ( ; name_iter != type_names.end() && type == NULL;
          name_iter++, type_iter++) {
        if (*name_iter == name)
          type = *type_iter;
      }
      if (type == NULL)
        throw InternalError("Unable to determine variable type", __FILE__, __LINE__);
      if (type->getType() != TypeDescription::MeshNodeVariable)    
        throw InternalError("Variable type is not MeshNodeVariable", __FILE__, __LINE__);
      // find the first timestep
      int ts = 0;
      while ((ts < (int)d_tstimes.size()) && (startTime > d_tstimes[ts]))
        ts++;

      // idx needs to be initialized before it is used in findMeshAndIndex.
      meshNodeIndex idx = 0;
      for ( ; (ts < (int)d_tstimes.size()) && (d_tstimes[ts] <= endTime); ts++) {
        // figure out what mesh contains the element. As far as I can tell,
        // nothing prevents this from changing between timesteps, so we have to
        // do this every time -- if that can't actually happen we might be able
        // to speed this up.
        Mesh* mesh = NULL;
        DomainP domain = queryDomain( ts);
        findMeshAndIndex(domain, mesh, idx, meshNodeID, matlIndex, ts);
        //    std::cerr <<" Mesh = 0x"<<hex<<mesh<<dec<<", index = "<<idx;
        if (mesh == NULL)
          throw VariableNotFoundInDomain(name,meshNodeID,matlIndex,
              "DataArchive::query", __FILE__, __LINE__);

        MeshNodeVariable<T> var;
        query(var, name, matlIndex, mesh, ts);
        //now find the index that corresponds to the meshNodeID
        //std::cerr <<" time = "<<t<<",  value = "<<var[idx]<<std::endl;
        values.push_back(var[idx]);

      }
      dbg << "DataArchive::query(values) completed in "
        << (SCIRun::Time::currentSeconds() - call_start) << " seconds\n";
    }  

  template<class T>
    void DataArchive::query(std::vector<T>& values, const std::string& name,
        int matlIndex, IntVector loc,
        double startTime, double endTime)
    {
      double call_start = SCIRun::Time::currentSeconds();

      std::vector<int> index;
      std::vector<double> times;
      queryTimesteps(index, times); // build timesteps if not already done

      // figure out what kind of variable we're looking for
      std::vector<std::string> type_names;
      std::vector<const TypeDescription*> type_descriptions;
      queryVariables(type_names, type_descriptions);
      const TypeDescription* type = NULL;
      std::vector<std::string>::iterator name_iter = type_names.begin();
      std::vector<const TypeDescription*>::iterator type_iter = type_descriptions.begin();
      for ( ; name_iter != type_names.end() && type == NULL;
          name_iter++, type_iter++) {
        if (*name_iter == name)
          type = *type_iter;
      }
      if (type == NULL)
        throw InternalError("Unable to determine variable type", __FILE__, __LINE__);

      // find the first timestep
      int ts = 0;
      while ((ts < (int)d_tstimes.size()) && (startTime > d_tstimes[ts]))
        ts++;

      for ( ; (ts < (int)d_tstimes.size()) && (d_tstimes[ts] <= endTime); ts++) {
        // figure out what mesh contains the element. As far as I can tell,
        // nothing prevents this from changing between timesteps, so we have to
        // do this every time -- if that can't actually happen we might be able
        // to speed this up.
        Mesh* mesh = NULL;
        DomainP domain = queryDomain(ts);

        switch (type->getType()) {
          case TypeDescription::MeshElementVariable:
            for (Level::const_meshIterator iter = domain->meshesBegin();
                  (iter != domain->meshesEnd()) && (mesh == NULL); iter++) {
              if ((*iter)->containsMeshElement(loc)) {
                mesh = *iter;
                // We found our mesh, quit looking.
                break;
              }
            }
            break;

          case TypeDescription::MeshNodeVariable:
            for (Level::const_meshIterator iter = domain->meshesBegin();
                  (iter != domain->meshesEnd()) && (mesh == NULL); iter++) {
                if ((*iter)->containsMeshNode(loc)) {
                  mesh = *iter;
                  break;
                }
              }
              break;

            default:
              std::cerr << "Variable of unsupported type for this element-based query: " << type->getType() << '\n';
              break;
          }
        }
        if (mesh == NULL) {
          throw VariableNotFoundInDomain(name,loc,matlIndex,"DataArchive::query", __FILE__, __LINE__);
        }

        switch (type->getType()) {
          case TypeDescription::MeshElementVariable: {
                                              MeshElementVariable<T> var;
                                              query(var, name, matlIndex, mesh, ts);
                                              values.push_back(var[loc]);
                                            } break;

          case TypeDescription::MeshNodeVariable: {
                                              NCVariable<T> var;
                                              query(var, name, matlIndex, mesh, ts);
                                              values.push_back(var[loc]);
                                            } break;

          default:
                                              // Dd: Is this correct?  Error here?
                                              break;
        }
        //std::cerr << "DataArchive::query:data extracted" << std::endl;
      }

      dbg << "DataArchive::query(values) completed in "
        << (SCIRun::Time::currentSeconds() - call_start) << " seconds\n";
    }
  
} // end namespace Matiti

#endif

