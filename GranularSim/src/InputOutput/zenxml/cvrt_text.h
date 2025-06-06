// *****************************************************************************
// * This file is part of the FreeFileSync project. It is distributed under    *
// * GNU General Public License: http://www.gnu.org/licenses/gpl-3.0           *
// * Copyright (C) Zenju (zenju AT freefilesync DOT org) - All Rights Reserved *
// *****************************************************************************

#ifndef CVRT_TEXT_H_018727339083427097434
#define CVRT_TEXT_H_018727339083427097434

#include "string_tools.h"
#include "utf.h"

namespace zen {
/**
\file
\brief Handle conversion of string-convertible types to and from std::string.

It is \b not required to call these functions directly. They are implicitly used
by zen::XmlElement::getValue(),
zen::XmlElement::setValue(), zen::XmlElement::getAttribute() and
zen::XmlElement::setAttribute().
\n\n
Conversions for the following user types are supported by default:
    - strings - std::string, std::wstring, char*, wchar_t*, char, wchar_t,
ect..., all STL-compatible-string-classes
    - numbers - int, double, float, bool, long, ect..., all built-in numbers
    - STL containers - std::map, std::set, std::vector, std::list, ect..., all
STL-compatible-containers
    - std::pair

You can add support for additional types via template specialization. \n\n
Specialize zen::readStruc() and zen::writeStruc() to enable conversion from
structured user types to XML elements.
Specialize zen::readText() and zen::writeText() to enable conversion from
string-convertible user types to std::string.
Prefer latter if possible since it does not only enable conversions from XML
elements to user data, but also from and to XML attributes.
\n\n
<b> Example: </b> type "bool"
\code
namespace zen
{
template <> inline
void writeText(const bool& value, std::string& output)
{
    output = value ? "true" : "false";
}

template <> inline
bool readText(const std::string& input, bool& value)
{
    const std::string tmp = trimCpy(input);
    if (tmp == "true")
        value = true;
    else if (tmp == "false")
        value = false;
    else
        return false;
    return true;
}
}
\endcode
*/

/// Convert text to user data - used by XML elements and attributes
/**
  \param input Input text.
  \param value Conversion target value.
  \return "true" if value was read successfully.
*/
template <class T>
bool readText(const std::string& input, T& value);
/// Convert user data into text - used by XML elements and attributes
/**
  \param value The value to be converted.
  \param output Output text.
*/
template <class T>
void writeText(const T& value, std::string& output);

/* Different classes of data types:

-----------------------------
| structured                |  readStruc/writeStruc - e.g. string-convertible
types, STL containers, std::pair, structured user types
| ------------------------- |
| | to-string-convertible | |  readText/writeText   - e.g. string-like types,
all built-in arithmetic numbers, bool
| | ---------------       | |
| | | string-like |       | |  utfCvrtTo            - e.g. std::string,
wchar_t*, char[], wchar_t, wxString, MyStringClass, ...
| | ---------------       | |
| ------------------------- |
-----------------------------
*/

//------------------------------ implementation
//-------------------------------------

// Conversion from arbitrary types to text (for use with XML elements and
// attributes)
enum TextType
{
  TEXT_TYPE_BOOL,
  TEXT_TYPE_NUMBER,
  TEXT_TYPE_STRING,
  TEXT_TYPE_OTHER,
};

template <class T>
struct GetTextType
  : StaticEnum<TextType, IsSameType<T, bool>::value
                           ? TEXT_TYPE_BOOL
                           : IsStringLike<T>::value
                               ? TEXT_TYPE_STRING
                               : // string before number to correctly handle
                               // char/wchar_t -> this was an issue with Loki
                               // only!
                               IsArithmetic<T>::value ? TEXT_TYPE_NUMBER : //
                                 TEXT_TYPE_OTHER>
{
};

//######################################################################################

template <class T, TextType type>
struct ConvertText;
/* -> expected interface
{
    void writeText(const T& value, std::string& output) const;
    bool readText(const std::string& input, T& value) const;
};
*/

// partial specialization: type bool
template <class T>
struct ConvertText<T, TEXT_TYPE_BOOL>
{
  void writeText(bool value, std::string& output) const
  {
    output = value ? "true" : "false";
  }
  bool readText(const std::string& input, bool& value) const
  {
    const std::string tmp = trimCpy(input);
    if (tmp == "true")
      value = true;
    else if (tmp == "false")
      value = false;
    else
      return false;
    return true;
  }
};

// partial specialization: handle conversion for all built-in arithmetic types!
template <class T>
struct ConvertText<T, TEXT_TYPE_NUMBER>
{
  void writeText(const T& value, std::string& output) const
  {
    output = numberTo<std::string>(value);
  }
  bool readText(const std::string& input, T& value) const
  {
    value = stringTo<T>(input);
    return true;
  }
};

// partial specialization: handle conversion for all string-like types!
template <class T>
struct ConvertText<T, TEXT_TYPE_STRING>
{
  void writeText(const T& value, std::string& output) const
  {
    output = utfCvrtTo<std::string>(value);
  }
  bool readText(const std::string& input, T& value) const
  {
    value = utfCvrtTo<T>(input);
    return true;
  }
};

// partial specialization: unknown type
template <class T>
struct ConvertText<T, TEXT_TYPE_OTHER>
{
  //###########################################################################################################################################
  static_assert(sizeof(T) == -1, "");
  /*
      ATTENTION: The data type T is yet unknown to the zen::Xml framework!

      Please provide a specialization for T of the following two functions in
     order to handle conversions to XML elements and attributes

      template <> void zen::writeText(const T& value, std::string& output)
      template <> bool zen::readText(const std::string& input, T& value)

      If T is structured and cannot be converted to a text representation
     specialize these two functions to allow at least for conversions to XML
     elements:

      template <> void zen::writeStruc(const T& value, XmlElement& output)
      template <> bool zen::readStruc(const XmlElement& input, T& value)
  */
  //###########################################################################################################################################
};

template <class T>
inline void
writeText(const T& value, std::string& output)
{
  ConvertText<T, GetTextType<T>::value>().writeText(value, output);
}

template <class T>
inline bool
readText(const std::string& input, T& value)
{
  return ConvertText<T, GetTextType<T>::value>().readText(input, value);
}
}

#endif // CVRT_TEXT_H_018727339083427097434
