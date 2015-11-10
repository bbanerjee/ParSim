/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-     Parresia Research Limited, New Zealand
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

#ifndef VAANGO_CCA_COMPONENTS_SCHEDULERS_KeyDatabase_H
#define VAANGO_CCA_COMPONENTS_SCHEDULERS_KeyDatabase_H

#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarLabelMatl.h>
#include <CCA/Components/Schedulers/HashVarLabelMatl.h>

#include <map>
#include <sci_hash_map.h>


namespace Uintah {

  template<class DomainType> class DWDatabase;

  /**************************************
     
     CLASS
       KeyDatabase
      
     GENERAL INFORMATION
      
       KeyDatabase.h
      
     KEYWORDS
      
     DESCRIPTION
      
     WARNING
      
  ****************************************/

  template<class DomainType>
  class KeyDatabase {

    template<class T> friend class DWDatabase;

  public:
    KeyDatabase();

    ~KeyDatabase();

    void clear();

    void insert(const VarLabel* label,
                int matlIndex,
                const DomainType* dom);

    int lookup(const VarLabel* label,
               int matlIndex,
               const DomainType* dom);

    void merge(const KeyDatabase<DomainType>& newDB);

  private:

    typedef hashmap<VarLabelMatl<DomainType>, int> keyDBtype;
    keyDBtype keys;
    int keycount;
  };

  template<class DomainType>
  KeyDatabase<DomainType>::KeyDatabase():keycount(0)
  {
  }

  template<class DomainType>
  KeyDatabase<DomainType>::~KeyDatabase()
  {
  }

  //______________________________________________________________________
  //
  template<class DomainType>
  int KeyDatabase<DomainType>::lookup(const VarLabel* label, int matlIndex, const DomainType* dom)
  {
    VarLabelMatl<DomainType> v(label, matlIndex, getRealDomain(dom));
    typename keyDBtype::const_iterator iter = keys.find(v);
    if (iter == keys.end()) {
      return -1;
    }
    else {
      return iter->second;
    }
  }

  //______________________________________________________________________
  //
  template<class DomainType>
  void KeyDatabase<DomainType>::merge(const KeyDatabase<DomainType>& newDB){
    for (auto keyiter = newDB.keys.begin(); keyiter != newDB.keys.end(); keyiter++) {
      typename keyDBtype::const_iterator iter = keys.find(keyiter->first);
      if (iter == keys.end()) {
        keys.insert(std::pair<VarLabelMatl<DomainType>, int>(keyiter->first, keycount++));
      }
    }
  }

  //______________________________________________________________________
  //
  template<class DomainType>
  void KeyDatabase<DomainType>::insert(const VarLabel* label, int matlIndex, const DomainType* dom)
  {
    VarLabelMatl<DomainType> v(label, matlIndex, getRealDomain(dom));
    typename keyDBtype::const_iterator iter = keys.find(v);
    if (iter == keys.end())
      keys.insert(std::pair<VarLabelMatl<DomainType>, int>(v, keycount++));
  }

  //______________________________________________________________________
  //
  template<class DomainType>
  void KeyDatabase<DomainType>::clear()
  {
    keys.clear();
    keycount = 0;
  }

}
#endif
