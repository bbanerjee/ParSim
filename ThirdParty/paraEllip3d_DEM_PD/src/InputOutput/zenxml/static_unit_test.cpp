// *****************************************************************************
// * This file is part of the FreeFileSync project. It is distributed under    *
// * GNU General Public License: http://www.gnu.org/licenses/gpl-3.0           *
// * Copyright (C) Zenju (zenju AT freefilesync DOT org) - All Rights Reserved *
// *****************************************************************************

#include <set>
#include <map>
#include <deque>
#include <vector>
#include <list>
#include <utility>
#include "xml.h"

using namespace zen;

namespace
{
void unit_test()
{
    class Dummy {};

    //compile time checks only

    static_assert(!IsStlContainer<wchar_t>                     ::value, "");
    static_assert(!IsStlContainer<char>                        ::value, "");
    static_assert(!IsStlContainer<Dummy>                       ::value, "");
    static_assert(!IsStlContainer<NullType>                    ::value, "");
    static_assert( IsStlContainer<std::set     <int>>          ::value, "");
    static_assert( IsStlContainer<std::deque   <float>>        ::value, "");
    static_assert( IsStlContainer<std::list    <size_t>>       ::value, "");
    static_assert( IsStlContainer<std::map     <double, char>> ::value, "");
    static_assert( IsStlContainer<std::multimap<short, double>>::value, "");
    static_assert( IsStlContainer<std::vector  <wchar_t>>      ::value, "");
    static_assert( IsStlPair     <std::pair<int, double>>      ::value, "");
    static_assert(!IsStlPair<Dummy>                            ::value, "");

    static_assert(!IsStringLike<Dummy>               ::value, "");
    static_assert(!IsStringLike<int>                 ::value, "");
    static_assert(!IsStringLike<double>              ::value, "");
    static_assert(!IsStringLike<short>               ::value, "");
    static_assert( IsStringLike<char>                ::value, "");
    static_assert( IsStringLike<const wchar_t>       ::value, "");
    static_assert( IsStringLike<const char>          ::value, "");
    static_assert( IsStringLike<wchar_t>             ::value, "");
    static_assert( IsStringLike<char*>               ::value, "");
    static_assert( IsStringLike<wchar_t*>            ::value, "");
    static_assert( IsStringLike<char* const>         ::value, "");
    static_assert( IsStringLike<wchar_t* const>      ::value, "");
    static_assert( IsStringLike<const char*>         ::value, "");
    static_assert( IsStringLike<const char* const>   ::value, "");
    static_assert( IsStringLike<const wchar_t*>      ::value, "");
    static_assert( IsStringLike<const wchar_t* const>::value, "");
    static_assert( IsStringLike<const char[4]>       ::value, "");
    static_assert( IsStringLike<const wchar_t[4]>    ::value, "");
    static_assert( IsStringLike<char[4]>             ::value, "");
    static_assert( IsStringLike<wchar_t[4]>          ::value, "");
    static_assert( IsStringLike<std::string>         ::value, "");
    static_assert( IsStringLike<std::wstring>        ::value, "");
    static_assert( IsStringLike<const std::string>   ::value, "");
    static_assert( IsStringLike<const std::wstring>  ::value, "");
    static_assert( IsStringLike<const std::string&>  ::value, "");
    static_assert( IsStringLike<const std::wstring&> ::value, "");
    static_assert( IsStringLike<std::string&>        ::value, "");
    static_assert( IsStringLike<std::wstring&>       ::value, "");

    static_assert(!IsSameType<GetCharType<int>                 ::Type, char>   ::value, "");
    static_assert(!IsSameType<GetCharType<double>              ::Type, char>   ::value, "");
    static_assert(!IsSameType<GetCharType<short>               ::Type, char>   ::value, "");
    static_assert( IsSameType<GetCharType<char>                ::Type, char>   ::value, "");
    static_assert( IsSameType<GetCharType<wchar_t>             ::Type, wchar_t>::value, "");
    static_assert( IsSameType<GetCharType<const char>          ::Type, char>   ::value, "");
    static_assert( IsSameType<GetCharType<const wchar_t>       ::Type, wchar_t>::value, "");
    static_assert( IsSameType<GetCharType<char*>               ::Type, char>   ::value, "");
    static_assert( IsSameType<GetCharType<wchar_t*>            ::Type, wchar_t>::value, "");
    static_assert( IsSameType<GetCharType<char* const>         ::Type, char>   ::value, "");
    static_assert( IsSameType<GetCharType<wchar_t* const>      ::Type, wchar_t>::value, "");
    static_assert( IsSameType<GetCharType<const char*>         ::Type, char>   ::value, "");
    static_assert( IsSameType<GetCharType<const char* const>   ::Type, char>   ::value, "");
    static_assert( IsSameType<GetCharType<const wchar_t*>      ::Type, wchar_t>::value, "");
    static_assert( IsSameType<GetCharType<const wchar_t* const>::Type, wchar_t>::value, "");
    static_assert( IsSameType<GetCharType<const char   [4]>    ::Type, char>   ::value, "");
    static_assert( IsSameType<GetCharType<const wchar_t[4]>    ::Type, wchar_t>::value, "");
    static_assert( IsSameType<GetCharType<char   [4]>          ::Type, char>   ::value, "");
    static_assert( IsSameType<GetCharType<wchar_t[4]>          ::Type, wchar_t>::value, "");
    static_assert( IsSameType<GetCharType<std::string>         ::Type, char>   ::value, "");
    static_assert( IsSameType<GetCharType<std::wstring>        ::Type, wchar_t>::value, "");
    static_assert( IsSameType<GetCharType<const std::string>   ::Type, char>   ::value, "");
    static_assert( IsSameType<GetCharType<const std::wstring>  ::Type, wchar_t>::value, "");
    static_assert( IsSameType<GetCharType<const std::string&>  ::Type, char>   ::value, "");
    static_assert( IsSameType<GetCharType<const std::wstring&> ::Type, wchar_t>::value, "");
    static_assert( IsSameType<GetCharType<std::string&>        ::Type, char>   ::value, "");
    static_assert( IsSameType<GetCharType<std::wstring&>       ::Type, wchar_t>::value, "");
}
}
