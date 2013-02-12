#ifndef MATITI_CONSTVARIABLE_H
#define MATITI_CONSTVARIABLE_H

#include <Mesh/constVariableBase.h>

namespace Matiti {

  template<class VariableBase, class Variable, class T, class Index> 
  class constVariable : public constVariableBase<VariableBase> {
  public:
    typedef T value_type;

    constVariable()
      : rep_() {}

    constVariable(const Variable& copy)
      : rep_(copy) {}

    constVariable<VariableBase, Variable, T, Index>&
    operator=(const constVariable<VariableBase, Variable, T, Index>& v)
    { copyPointer(v.rep_); return *this; }

    constVariable<VariableBase, Variable, T, Index>& operator=(const Variable& v)
    { copyPointer(v); return *this; }
    
    constVariableBase<VariableBase>&
    operator=(const constVariableBase<VariableBase>& v)
    {
      const constVariable<VariableBase, Variable, T, Index>* cvp =
        dynamic_cast<const constVariable<VariableBase, Variable, T, Index>*>(&v);
      //ASSERT(cvp != 0);
      copyPointer(cvp->rep_);
      return *this;
    }

    constVariableBase<VariableBase>&
    operator=(const VariableBase& v)
    { copyPointer(v); return *this; }

    virtual ~constVariable() {}

    operator const Variable&() const
    { return this->rep_; }
    virtual const VariableBase& getBaseRep() const
    { return this->rep_; }

    // It's ok for a constVariable to copyPointer of a const variable
    // (even though a non-const variable can't).
    virtual void copyPointer(const VariableBase& copy)
    { this->rep_.copyPointer(const_cast<VariableBase&>(copy)); }

    virtual const VariableBase* clone() const
      // need to cast it if it is a GridVariable
    { return dynamic_cast<const VariableBase*>(this->rep_.clone()); }

    virtual VariableBase* cloneType() const
      // need to cast it if it is a GridVariable
    { return dynamic_cast<VariableBase*>(this->rep_.cloneType()); }

    inline const T& operator[](Index idx) const
    { return this->rep_[idx]; }
     
  protected:
    Variable rep_;
  };

} // end namespace Matiti


#endif

