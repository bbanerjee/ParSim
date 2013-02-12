#ifndef MATITI_VARIABLE_H
#define MATITI_VARIABLE_H

#include <string>

namespace Matiti {

  class RefCounted;

  class Variable {

  public:
    virtual ~Variable();
  
    //marks a variable as invalid (for example, it is in the process of receiving mpi)
    void setValid() { d_valid=true;} 
    void setInvalid() { d_valid=false;} 
    //returns if a variable is marked valid or invalid
    bool isValid() const {return d_valid;}

    virtual void emitNormal(std::ostream& out) = 0;
    virtual void readNormal(std::istream& in) = 0;

    virtual void allocate() = 0;

    virtual void getSizeInfo(std::string& elems, unsigned long& totsize, void*& ptr) const = 0;
 
    virtual void copyPointer(Variable&) = 0;

    virtual RefCounted* getRefCounted() = 0;

  protected:
    Variable();

  private:    
    Variable(const Variable&);
    Variable& operator=(const Variable&);

    //signals of the variable is valid, an mpi variable is not valid until mpi has been recieved
    bool d_valid;
 };

} // End namespace Matiti

#endif
