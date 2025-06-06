// *****************************************************************************
// * This file is part of the FreeFileSync project. It is distributed under    *
// * GNU General Public License: http://www.gnu.org/licenses/gpl-3.0           *
// * Copyright (C) Zenju (zenju AT freefilesync DOT org) - All Rights Reserved *
// *****************************************************************************

#ifndef CVRT_STRUC_H_018727409908342709743
#define CVRT_STRUC_H_018727409908342709743

#include "dom.h"

namespace zen {
/**
\file
\brief Handle conversion of arbitrary types to and from XML elements.
See comments in cvrt_text.h
*/

/// Convert XML element to structured user data
/**
  \param input The input XML element.
  \param value Conversion target value.
  \return "true" if value was read successfully.
*/
template <class T>
bool readStruc(const XmlElement& input, T& value);
/// Convert structured user data into an XML element
/**
  \param value The value to be converted.
  \param output The output XML element.
*/
template <class T>
void writeStruc(const T& value, XmlElement& output);

//------------------------------ implementation
//-------------------------------------
namespace impl_2384343 {
ZEN_INIT_DETECT_MEMBER_TYPE(value_type);
ZEN_INIT_DETECT_MEMBER_TYPE(iterator);
ZEN_INIT_DETECT_MEMBER_TYPE(const_iterator);

ZEN_INIT_DETECT_MEMBER(begin) //
ZEN_INIT_DETECT_MEMBER(end) // we don't know the exact declaration of the member
                            // attribute: may be in a base class!
ZEN_INIT_DETECT_MEMBER(insert) //
}

template <typename T>
struct IsStlContainer
  : StaticBool<impl_2384343::HasMemberType_value_type<T>::value &&
               impl_2384343::HasMemberType_iterator<T>::value &&
               impl_2384343::HasMemberType_const_iterator<T>::value &&
               impl_2384343::HasMember_begin<T>::value &&
               impl_2384343::HasMember_end<T>::value &&
               impl_2384343::HasMember_insert<T>::value>
{
};

namespace impl_2384343 {
ZEN_INIT_DETECT_MEMBER_TYPE(first_type);
ZEN_INIT_DETECT_MEMBER_TYPE(second_type);

ZEN_INIT_DETECT_MEMBER(first)  // we don't know the exact declaration of the
                               // member attribute: may be in a base class!
ZEN_INIT_DETECT_MEMBER(second) //
}

template <typename T>
struct IsStlPair
  : StaticBool<impl_2384343::HasMemberType_first_type<T>::value &&
               impl_2384343::HasMemberType_second_type<T>::value &&
               impl_2384343::HasMember_first<T>::value &&
               impl_2384343::HasMember_second<T>::value>
{
};

//######################################################################################

// Conversion from arbitrary types to an XML element
enum ValueType
{
  VALUE_TYPE_STL_CONTAINER,
  VALUE_TYPE_STL_PAIR,
  VALUE_TYPE_OTHER,
};

template <class T>
struct GetValueType
  : StaticEnum<
      ValueType,
      GetTextType<T>::value != TEXT_TYPE_OTHER
        ? VALUE_TYPE_OTHER
        : // some string classes are also STL containers, so check this first
        IsStlContainer<T>::value
          ? VALUE_TYPE_STL_CONTAINER
          : IsStlPair<T>::value ? VALUE_TYPE_STL_PAIR : VALUE_TYPE_OTHER>
{
};

template <class T, ValueType type>
struct ConvertElement;
/* -> expected interface
{
    void writeStruc(const T& value, XmlElement& output) const;
    bool readStruc(const XmlElement& input, T& value) const;
};
*/

// partial specialization: handle conversion for all STL-container types!
template <class T>
struct ConvertElement<T, VALUE_TYPE_STL_CONTAINER>
{
  void writeStruc(const T& value, XmlElement& output) const
  {
    for (const typename T::value_type& childVal : value) {
      XmlElement& newChild = output.addChild("Item");
      zen::writeStruc(childVal, newChild);
    }
  }
  bool readStruc(const XmlElement& input, T& value) const
  {
    bool success = true;
    value.clear();

    auto itPair = input.getChildren("Item");
    for (auto it = itPair.first; it != itPair.second; ++it) {
      typename T::value_type
        childVal; // MSVC 2010 bug: cannot put this into a lambda body
      if (zen::readStruc(*it, childVal))
        value.insert(value.end(), childVal);
      else
        success = false;
    }
    return success;
  }
};

// partial specialization: handle conversion for std::pair
template <class T>
struct ConvertElement<T, VALUE_TYPE_STL_PAIR>
{
  void writeStruc(const T& value, XmlElement& output) const
  {
    XmlElement& child1 = output.addChild("one"); // don't use "1st/2nd", this
                                                 // will confuse a few pedantic
                                                 // XML parsers
    zen::writeStruc(value.first, child1);

    XmlElement& child2 = output.addChild("two");
    zen::writeStruc(value.second, child2);
  }
  bool readStruc(const XmlElement& input, T& value) const
  {
    bool success = true;
    const XmlElement* child1 = input.getChild("one");
    if (!child1 || !zen::readStruc(*child1, value.first))
      success = false;

    const XmlElement* child2 = input.getChild("two");
    if (!child2 || !zen::readStruc(*child2, value.second))
      success = false;

    return success;
  }
};

// partial specialization: not a pure structured type, try text conversion
// (thereby respect user specializations of writeText()/readText())
template <class T>
struct ConvertElement<T, VALUE_TYPE_OTHER>
{
  void writeStruc(const T& value, XmlElement& output) const
  {
    std::string tmp;
    writeText(value, tmp);
    output.setValue(tmp);
  }
  bool readStruc(const XmlElement& input, T& value) const
  {
    std::string rawStr;
    input.getValue(rawStr);
    return readText(rawStr, value);
  }
};

template <class T>
inline void
writeStruc(const T& value, XmlElement& output)
{
  ConvertElement<T, GetValueType<T>::value>().writeStruc(value, output);
}

template <class T>
inline bool
readStruc(const XmlElement& input, T& value)
{
  return ConvertElement<T, GetValueType<T>::value>().readStruc(input, value);
}
}

#endif // CVRT_STRUC_H_018727409908342709743
