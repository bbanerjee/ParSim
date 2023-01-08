/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#ifndef VAANGO_CCA_COMPONENTS_SCHEDULERS_HashVarLableMatl_H
#define VAANGO_CCA_COMPONENTS_SCHEDULERS_HashVarLableMatl_H

//
// Hash function for VarLabelMatl
//
#ifdef HAVE_GNU_HASHMAP

  namespace __gnu_cxx
  {
//    using Uintah::KeyDatabase;
 //   using Uintah::DWDatabase;
    using Uintah::VarLabelMatl;
    template <class DomainType>
    struct hash<VarLabelMatl<DomainType> > : public std::unary_function<VarLabelMatl<DomainType>, size_t>
    {
      size_t operator()(const VarLabelMatl<DomainType>& v) const
      {
        size_t h=0;
        char *str =const_cast<char*> (v.label_->getName().data());
        while (int c = *str++) h = h*7+c;
        return ( ( ((size_t)v.label_) << (sizeof(size_t)/2) ^ ((size_t)v.label_) >> (sizeof(size_t)/2) )
                 ^ (size_t)v.domain_ ^ (size_t)v.matlIndex_ );
      }
    };
  }

#elif HAVE_TR1_HASHMAP || HAVE_C11_HASHMAP 

  namespace std {
#if HAVE_TR1_HASHMAP 
    namespace tr1 {
#endif 
//      using Uintah::KeyDatabase;
//      using Uintah::DWDatabase;
      using Uintah::VarLabelMatl;
      template <class DomainType>
      struct hash<VarLabelMatl<DomainType> > : public unary_function<VarLabelMatl<DomainType>, size_t>
      {
        size_t operator()(const VarLabelMatl<DomainType>& v) const
        {
          size_t h=0;
          char *str =const_cast<char*> (v.label_->getName().data());
          while (int c = *str++) h = h*7+c;
          return ( ( ((size_t)v.label_) << (sizeof(size_t)/2) ^ ((size_t)v.label_) >> (sizeof(size_t)/2) )
                   ^ (size_t)v.domain_ ^ (size_t)v.matlIndex_ );
        }
      };
#if HAVE_TR1_HASHMAP 
    } // end namespace tr1
#endif 
  } // end namespace std

#endif // #ifdef HAVE_GNU_HASHMAP

#endif //VAANGO_CCA_COMPONENTS_SCHEDULERS_HashVarLableMatl_H
