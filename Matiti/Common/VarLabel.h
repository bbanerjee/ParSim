#ifndef MATITI_VARLABEL_H
#define MATITI_VARLABEL_H

#include <string>
#include <iosfwd>
#include <Common/RefCounted.h>

namespace Matiti {

  using std::string;

  class TypeDescription;

  class VarLabel : public RefCounted {
  public:
    enum VarType {
      Normal,
      PositionVariable
    };

    // Ensure the uniqueness of VarLabel names (same name, same object).
    static VarLabel* create(const string&, const TypeDescription*,
                            VarType vartype = Normal);

    static bool destroy(const VarLabel* label);

    inline const string& getName() const {
      return d_name;
    }
    string getFullName(int matlIndex) const;
    bool isPositionVariable() const {
      return d_vartype == PositionVariable;
    }

    const TypeDescription* typeDescription() const {
      return d_td;
    }

    static VarLabel* find(string name);

    class Compare {
    public:
      inline bool operator()(const VarLabel* v1, const VarLabel* v2) const {
        // because of uniqueness, we can use pointer comparisons
        //return v1 < v2;
        // No we cannot, because we need the order to be the same on different processes
        if(v1 == v2)
          return false;
        return v1->getName() < v2->getName();
      }
    private:
    };
    
    bool equals(const VarLabel* v2) const {
      // because of uniqueness, we can use pointer comparisons
      return this == v2;
      /* old way
         if(this == v2)
         return true;
         return getName() == v2->getName();
      */
    }

    static void printAll(); // for debugging
     
    string                 d_name;

     friend std::ostream & operator<<( std::ostream & out, const Matiti::VarLabel & vl );

  private:
    // You must use VarLabel::create.
    VarLabel(const string&, const TypeDescription*, VarType vartype);
    // You must use destroy.
    ~VarLabel();   
     
    const TypeDescription* d_td;
    VarType                d_vartype;
    
    VarLabel(const VarLabel&);
    VarLabel& operator=(const VarLabel&);
  };
} // End namespace Matiti

#endif
