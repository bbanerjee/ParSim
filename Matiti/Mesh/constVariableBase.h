#ifndef MATITI_CONSTVARIABLEBASE_H
#define MATITI_CONSTVARIABLEBASE_H


namespace Matiti {

  template<class VariableBase> 
  class constVariableBase {
  public:   
    virtual ~constVariableBase() {}

    virtual constVariableBase& operator=(const VariableBase&) = 0;
    virtual constVariableBase& operator=(const constVariableBase&) = 0;

    operator const VariableBase&() const
    { return getBaseRep(); }
   
    virtual const VariableBase& getBaseRep() const = 0;
    
    virtual void copyPointer(const VariableBase& copy) = 0;

    virtual const VariableBase* clone() const = 0;
    virtual VariableBase* cloneType() const = 0;

    virtual const TypeDescription* virtualGetTypeDescription() const = 0;

  protected:
    constVariableBase() {}
    constVariableBase(const constVariableBase&) {}
  private:
  };

} // end namespace Matiti


#endif

