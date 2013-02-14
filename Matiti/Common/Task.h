#ifndef MATITI_TASK_H
#define MATITI_TASK_H

#include <Common/VarLabel.h>
#include <Common/ComputeSet.h>
#include <Common/DataWarehouseP.h>
#include <Common/constHandle.h>

#include <map>
#include <set>
#include <vector>
#include <string>
#include <iostream>

namespace Matiti {

  class DataWarehouse;

  class Task {

  protected:

    // base CPU Action class
    class ActionBase {
    public:
      virtual ~ActionBase();
      virtual void doit(const MaterialSubset* matls,
                        DataWarehouse* fromDW,
                        DataWarehouse* toDW) = 0;
    };


  private:

    // begin CPU Action constructors
    template<class T>
    class Action : public ActionBase {
      
      T* ptr;
      void (T::*pmf)(const MaterialSubset* matls,
                     DataWarehouse* fromDW,
                     DataWarehouse* toDW);
    public: // class Action
      Action(T* ptr,
             void (T::*pmf)(const MaterialSubset* matls,
                            DataWarehouse* fromDW,
                            DataWarehouse*toDW) )
             : ptr(ptr), pmf(pmf) {}
      virtual ~Action() {}
      
      virtual void doit(const MaterialSubset* matls,
                        DataWarehouse* fromDW,
                        DataWarehouse* toDW) {
        (ptr->*pmf)(matls, fromDW, toDW);
      }
    }; // end class Action

    template<class T, class Arg1>
    class Action1 : public ActionBase {

      T* ptr;
      void (T::*pmf)(const MaterialSubset* matls,
                     DataWarehouse* fromDW,
                     DataWarehouse* toDW,
                     Arg1 arg1);
      Arg1 arg1;
    public: // class Action1
      Action1( T* ptr,
               void (T::*pmf)(const MaterialSubset* matls,
                              DataWarehouse* fromDW,
                              DataWarehouse* toDW,
                              Arg1 arg1),
               Arg1 arg1 )
        : ptr(ptr), pmf(pmf), arg1(arg1) {}
      virtual ~Action1() {}

      //////////
      // Insert Documentation Here:
      virtual void doit(const MaterialSubset* matls,
                        DataWarehouse* fromDW,
                        DataWarehouse* toDW) {
        (ptr->*pmf)(matls, fromDW, toDW, arg1);
      }
    }; // end class Action1

    template<class T, class Arg1, class Arg2>
    class Action2 : public ActionBase {

      T* ptr;
      void (T::*pmf)(const MaterialSubset* matls,
                     DataWarehouse* fromDW,
                     DataWarehouse* toDW,
                     Arg1 arg1, Arg2 arg2);
      Arg1 arg1;
      Arg2 arg2;
    public: // class Action2
      Action2( T* ptr,
               void (T::*pmf)(const MaterialSubset* matls,
                              DataWarehouse*,
                              DataWarehouse*,
                              Arg1, Arg2),
               Arg1 arg1, Arg2 arg2 )
        : ptr(ptr), pmf(pmf), arg1(arg1), arg2(arg2) {}
      virtual ~Action2() {}

      //////////
      // Insert Documentation Here:
      virtual void doit(const MaterialSubset* matls,
                        DataWarehouse* fromDW,
                        DataWarehouse* toDW) {
        (ptr->*pmf)(matls, fromDW, toDW, arg1, arg2);
      }
    }; // end class Action2

    template<class T, class Arg1, class Arg2, class Arg3>
    class Action3 : public ActionBase {

      T* ptr;
      void (T::*pmf)(const MaterialSubset* matls,
                     DataWarehouse* fromDW,
                     DataWarehouse* toDW,
                     Arg1 arg1, Arg2 arg2, Arg3 arg3);
      Arg1 arg1;
      Arg2 arg2;
      Arg3 arg3;
    public: // class Action3
      Action3( T* ptr,
               void (T::*pmf)(const MaterialSubset* matls,
                              DataWarehouse* fromDW,
                              DataWarehouse* toDW,
                              Arg1, Arg2, Arg3),
               Arg1 arg1, Arg2 arg2, Arg3 arg3 )
        : ptr(ptr), pmf(pmf), arg1(arg1), arg2(arg2), arg3(arg3) {}
      virtual ~Action3() {}

      //////////
      // Insert Documentation Here:
      virtual void doit(const MaterialSubset* matls,
                        DataWarehouse* fromDW,
                        DataWarehouse* toDW) {
        (ptr->*pmf)(matls, fromDW, toDW, arg1, arg2, arg3);
      }
    }; // end Action3

    template<class T, class Arg1, class Arg2, class Arg3, class Arg4>
    class Action4 : public ActionBase {

      T* ptr;
      void (T::*pmf)(const MaterialSubset* matls,
                     DataWarehouse* fromDW,
                     DataWarehouse* toDW,
                     Arg1 arg1, Arg2 arg2, Arg3 arg3, Arg4 arg4);
      Arg1 arg1;
      Arg2 arg2;
      Arg3 arg3;
      Arg4 arg4;
    public: // class Action4
      Action4( T* ptr,
               void (T::*pmf)(const MaterialSubset* matls,
                              DataWarehouse* fromDW,
                              DataWarehouse* toDW,
                              Arg1, Arg2, Arg3, Arg4),
               Arg1 arg1, Arg2 arg2, Arg3 arg3, Arg4 arg4)
        : ptr(ptr), pmf(pmf), arg1(arg1), arg2(arg2), arg3(arg3), arg4(arg4) {}
      virtual ~Action4() {}

      //////////
      // Insert Documentation Here:
      virtual void doit(const MaterialSubset* matls,
                        DataWarehouse* fromDW,
                        DataWarehouse* toDW) {
        (ptr->*pmf)(matls, fromDW, toDW, arg1, arg2, arg3, arg4);
      }
    }; // end Action4

    template<class T, class Arg1, class Arg2, class Arg3, class Arg4, class Arg5>
    class Action5 : public ActionBase {

      T* ptr;
      void (T::*pmf)(const MaterialSubset* matls,
                     DataWarehouse* fromDW,
                     DataWarehouse* toDW,
                     Arg1 arg1, Arg2 arg2, Arg3 arg3, Arg4 arg4, Arg5 arg5);
      Arg1 arg1;
      Arg2 arg2;
      Arg3 arg3;
      Arg4 arg4;
      Arg5 arg5;
    public: // class Action5
      Action5( T* ptr,
               void (T::*pmf)(const MaterialSubset* matls,
                              DataWarehouse* fromDW,
                              DataWarehouse* toDW,
                              Arg1, Arg2, Arg3, Arg4, Arg5),
               Arg1 arg1, Arg2 arg2, Arg3 arg3, Arg4 arg4, Arg5 arg5)
        : ptr(ptr), pmf(pmf), arg1(arg1), arg2(arg2), arg3(arg3), arg4(arg4), arg5(arg5) {}
      virtual ~Action5() {}

      //////////
      // Insert Documentation Here:
      virtual void doit(const MaterialSubset* matls,
                        DataWarehouse* fromDW,
                        DataWarehouse* toDW) {
        (ptr->*pmf)(matls, fromDW, toDW, arg1, arg2, arg3, arg4, arg5);
      }
    }; // end Action5
    // end CPU Action constructors

  public: // class Task
    
    enum WhichDW {
      OldDW=0,
      NewDW=1,
      ParentOldDW=2,
      ParentNewDW=3,
      TotalDWs=4
    };

    enum {
      NoDW = -1,
      InvalidDW = -2
    };
    
    enum TaskType {
      Normal,
      Reduction,
      InitialSend,
      Output
    };
    
    Task(const std::string& taskName, TaskType type)
      :  d_taskName(taskName),
         d_action(0),
    {
      d_tasktype = type;
      initialize();
    }
    
    // begin CPU Task constructors
    template<class T>
    Task(const std::string& taskName,
         T* ptr,
         void (T::*pmf)(const MaterialSubset* matls,
                        DataWarehouse* fromDW,
                        DataWarehouse* toDW))
      : d_taskName(taskName),
        d_action(new Action<T>(ptr, pmf)),
    {
      d_tasktype = Normal;
      initialize();
    }

    template<class T, class Arg1>
    Task(const std::string& taskName,
         T* ptr,
         void (T::*pmf)(const MaterialSubset* matls,
                        DataWarehouse* fromDW,
                        DataWarehouse* toDW,
                        Arg1 arg1),
         Arg1 arg1)
      : d_taskName(taskName),
        d_action(new Action1<T, Arg1>(ptr, pmf, arg1)),
    {
      d_tasktype = Normal;
      initialize();
    }

    template<class T, class Arg1, class Arg2>
    Task(const std::string& taskName,
         T* ptr,
         void (T::*pmf)(const MaterialSubset* matls,
                        DataWarehouse* fromDW,
                        DataWarehouse* toDW,
                        Arg1, Arg2),
         Arg1 arg1, Arg2 arg2)
      : d_taskName(taskName),
        d_action(new Action2<T, Arg1, Arg2>(ptr, pmf, arg1, arg2)),
    {
      d_tasktype = Normal;
      initialize();
    }

    template<class T, class Arg1, class Arg2, class Arg3>
    Task(const std::string& taskName,
         T* ptr,
         void (T::*pmf)(const MaterialSubset* matls,
                        DataWarehouse* fromDW,
                        DataWarehouse* toDW,
                        Arg1, Arg2, Arg3),
         Arg1 arg1, Arg2 arg2, Arg3 arg3)
      : d_taskName(taskName),
        d_action(new Action3<T, Arg1, Arg2, Arg3>(ptr, pmf, arg1, arg2, arg3)),
    {
      d_tasktype = Normal;
      initialize();
    }

    template<class T, class Arg1, class Arg2, class Arg3, class Arg4>
    Task(const std::string& taskName,
         T* ptr,
         void (T::*pmf)(const MaterialSubset* matls,
                        DataWarehouse* fromDW,
                        DataWarehouse* toDW,
                        Arg1, Arg2, Arg3, Arg4),
         Arg1 arg1, Arg2 arg2, Arg3 arg3, Arg4 arg4)
      : d_taskName(taskName),
        d_action(new Action4<T, Arg1, Arg2, Arg3, Arg4>(ptr, pmf, arg1, arg2, arg3, arg4)),
    {
      d_tasktype = Normal;
      initialize();
    }

    template<class T, class Arg1, class Arg2, class Arg3, class Arg4, class Arg5>
    Task(const std::string& taskName,
         T* ptr,
         void (T::*pmf)(const MaterialSubset* matls,
                        DataWarehouse* fromDW,
                        DataWarehouse* toDW,
                        Arg1, Arg2, Arg3, Arg4, Arg5),
         Arg1 arg1, Arg2 arg2, Arg3 arg3, Arg4 arg4, Arg5 arg5)
      : d_taskName(taskName),
        d_action(new Action5<T, Arg1, Arg2, Arg3, Arg4, Arg5>(ptr, pmf, arg1, arg2, arg3, arg4, arg5)),
    {
      d_tasktype = Normal;
      initialize();
    }
    // end CPU Task constructors

    void initialize();
    
    virtual ~Task();
    
    void hasSubScheduler(bool state = true);
    inline bool getHasSubScheduler() const { return d_hasSubScheduler; }
    
    enum MaterialDomainSpec {
      NormalDomain,  // <- Normal/default setting
      OutOfDomain,   // <- Require things from all material 
    };

    //////////
    // Most general case
    void requires(WhichDW, const VarLabel*, const MaterialSubset* matls, 
                  MaterialDomainSpec matls_dom)

    void requires(WhichDW, const VarLabel*, const MaterialSubset* matls)
    
    void requires(WhichDW, const VarLabel*) 
    
    //////////
    // Most general case
    void computes(const VarLabel*, const MaterialSubset* matls, 
                  MaterialDomainSpec matls_domain);
    
    void computes(const VarLabel*, const MaterialSubset* matls);

    void computes(const VarLabel*);
    
     //////////
    // Most general case
    void modifies(const VarLabel*, const MaterialSubset* matls, 
                  MaterialDomainSpec matls_domain)
    
    void modifies(const VarLabel*, const MaterialSubset* matls);
    
    void modifies(const VarLabel*); 

    //////////
    // Tells the task to actually execute the function assigned to it.
    virtual void doit(const MaterialSubset*, vector<DataWarehouseP>& dws);

    inline const MaterialSet* getMaterialSet() const {
      return matl_set;
    }
    
    struct Edge;
      
    int d_phase;  //synchronized phase id, for dynamic task scheduling
    std::set<Task*> childTasks;
    std::set<Task*> allChildTasks;
    
    enum DepType {
      Modifies, Computes, Requires
    };
    
    struct  Dependency {
      Dependency* next;
      DepType deptype;
      Task* task;
      const VarLabel*  var;
      const MaterialSubset* matls;
      Edge* req_head;   // Used in compiling the task graph.
      Edge* req_tail;
      Edge* comp_head;
      Edge* comp_tail;
      MaterialDomainSpec matls_dom;
      WhichDW whichdw;  // Used only by Requires
      
      // in the multi-TG construct, this will signify that the required
      // var will be constructed by the old TG
      int numGhostCells;
      int mapDataWarehouse() const {
        return task->mapDataWarehouse(whichdw);
      }
      
      Dependency(DepType deptype, Task* task, WhichDW dw, const VarLabel* var,
                 bool oldtg, 
                 const MaterialSubset* matls,
                 MaterialDomainSpec matls_dom = NormalDomain);
      ~Dependency();
      inline void addComp(Edge* edge);
      inline void addReq(Edge* edge);

      constHandle<MaterialSubset>
      getMaterialsUnderDomain(const MaterialSubset* domainMaterials) const;

    private:
      Dependency();
      Dependency& operator=(const Dependency& copy);
      Dependency(const Dependency&);
    }; // end struct Dependency
    
    
    struct  Edge {
      const Dependency* comp;
      Edge* compNext;
      const Dependency* req;
      Edge* reqNext;
      inline Edge(const Dependency* comp, const Dependency * req)
        : comp(comp), compNext(0), req(req), reqNext(0)
      {
      }
    };

    typedef std::multimap<const VarLabel*, Dependency*, VarLabel::Compare> DepMap;
    
    const Dependency* getComputes() const {
      return comp_head;
    }
    const Dependency* getRequires() const {
      return req_head;
    }
    const Dependency* getModifies() const {
      return mod_head;
    }
    
    Dependency* getComputes() {
      return comp_head;
    }
    Dependency* getRequires() {
      return req_head;
    }
    Dependency* getModifies() {
      return mod_head;
    }

    // finds if it computes or modifies var
    bool hasComputes(const VarLabel* var, int matlIndex) const;

    // finds if it requires or modifies var
    bool hasRequires(const VarLabel* var, int matlIndex, 
                     WhichDW dw) const;

    // finds if it modifies var
    bool hasModifies(const VarLabel* var, int matlIndex) const;

    bool isReductionTask() const {
      return d_tasktype == Reduction;
    }
    
    void setType(TaskType tasktype) {
      d_tasktype = tasktype;
    }
    TaskType getType() const {
      return d_tasktype;
    }
    
    //////////
    // Prints out information about the task...
    void display( std::ostream & out ) const;
    
    //////////
    // Prints out all information about the task, including dependencies
    void displayAll( std::ostream & out ) const;
    
    int mapDataWarehouse(WhichDW dw) const;
    DataWarehouse* mapDataWarehouse(WhichDW dw, vector<DataWarehouseP>& dws) const;

    int getSortedOrder() const {
      return sortedOrder;
    }

    void setSortedOrder(int order) {
      sortedOrder = order;
    }

    void setMapping(int dwmap[TotalDWs]);

    void setSets(const MaterialSet* matls);

    
  private: // class Task
    Dependency* isInDepMap(const DepMap& depMap, const VarLabel* var,
                           int matlIndex) const;
    
    //////////
    // Insert Documentation Here:
    std::string d_taskName;

protected:
    ActionBase*    d_action;

private:
    Dependency* comp_head;
    Dependency* comp_tail;
    Dependency* req_head;
    Dependency* req_tail;
    Dependency* mod_head;
    Dependency* mod_tail;

    DepMap d_requiresOldDW;
    DepMap d_computes; // also contains modifies
    DepMap d_requires; // also contains modifies
    DepMap d_modifies;
    
    const MaterialSet* matl_set;
    
    bool     d_hasSubScheduler;
    TaskType d_tasktype;
    
    Task(const Task&);
    Task& operator=(const Task&);

    static const MaterialSubset* getGlobalMatlSubset();
    static MaterialSubset* globalMatlSubset;

    int dwmap[TotalDWs];
    int sortedOrder;

     friend std::ostream & operator << ( std::ostream & out, const Matiti::Task & task );
     friend std::ostream & operator << ( std::ostream & out, const Matiti::Task::TaskType & tt );
     friend std::ostream & operator << ( std::ostream & out, const Matiti::Task::Dependency & dep );

  }; // end class Task
  
  inline void Task::Dependency::addComp(Edge* edge)
    {
      if(comp_tail)
        comp_tail->compNext=edge;
      else
        comp_head=edge;
      comp_tail=edge;
    }
  inline void Task::Dependency::addReq(Edge* edge)
    {
      if(req_tail)
        req_tail->reqNext=edge;
      else
        req_head=edge;
      req_tail=edge;
    }

} // End namespace Matiti

// This mus tbe at the bottom
#include <CCA/Ports/DataWarehouse.h>

#endif
