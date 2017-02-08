// *****************************************************************************
// * This file is part of the FreeFileSync project. It is distributed under    *
// * GNU General Public License: http://www.gnu.org/licenses/gpl-3.0           *
// * Copyright (C) Zenju (zenju AT freefilesync DOT org) - All Rights Reserved *
// *****************************************************************************

#ifndef UTF_H_01832479146991573473545
#define UTF_H_01832479146991573473545

#include <cstdint>
#include <iterator>
#include "string_tools.h" //copyStringTo

namespace zen
{
//convert all(!) char- and wchar_t-based "string-like" objects applying a UTF8 conversions (but only if necessary!)
template <class TargetString, class SourceString>
TargetString utfCvrtTo(const SourceString& str);

const char BYTE_ORDER_MARK_UTF8[] = "\xEF\xBB\xBF";

template <class CharString>
bool isValidUtf8(const CharString& str); //check for UTF-8 encoding errors

//---- explicit conversion: wide <-> utf8 ----
template <class CharString, class WideString>
CharString wideToUtf8(const WideString& str); //example: std::string tmp = wideToUtf8<std::string>(L"abc");

template <class WideString, class CharString>
WideString utf8ToWide(const CharString& str); //std::wstring tmp = utf8ToWide<std::wstring>("abc");

//access unicode characters in UTF-encoded string (char- or wchar_t-based)
template <class UtfString>
size_t unicodeLength(const UtfString& str); //return number of code points for UTF-encoded string

template <class UtfString>
size_t findUnicodePos(const UtfString& str, size_t unicodePos); //return position of unicode char in UTF-encoded string


















//----------------------- implementation ----------------------------------
namespace implementation
{
using CodePoint = std::uint32_t;
using Char16    = std::uint16_t;
using Char8     = unsigned char;

const CodePoint LEAD_SURROGATE      = 0xd800;
const CodePoint TRAIL_SURROGATE     = 0xdc00; //== LEAD_SURROGATE_MAX + 1
const CodePoint TRAIL_SURROGATE_MAX = 0xdfff;

const CodePoint REPLACEMENT_CHAR    = 0xfffd;
const CodePoint CODE_POINT_MAX      = 0x10ffff;


template <class Function> inline
void codePointToUtf16(CodePoint cp, Function writeOutput) //"writeOutput" is a unary function taking a Char16
{
    //http://en.wikipedia.org/wiki/UTF-16

    if (cp < LEAD_SURROGATE)
        writeOutput(static_cast<Char16>(cp));
    else if (cp <= TRAIL_SURROGATE_MAX) //invalid code point
        codePointToUtf16(REPLACEMENT_CHAR, writeOutput); //resolves to 1-character utf16
    else if (cp < 0x10000)
        writeOutput(static_cast<Char16>(cp));
    else if (cp <= CODE_POINT_MAX)
    {
        cp -= 0x10000;
        writeOutput(LEAD_SURROGATE  + static_cast<Char16>(cp >> 10));
        writeOutput(TRAIL_SURROGATE + static_cast<Char16>(cp & 0x3ff));
    }
    else //invalid code point
        codePointToUtf16(REPLACEMENT_CHAR, writeOutput); //resolves to 1-character utf16
}


inline
size_t getUtf16Len(Char16 ch) //ch must be first code unit! returns 0 on error!
{
    if (ch < LEAD_SURROGATE)
        return 1;
    else if (ch < TRAIL_SURROGATE)
        return 2;
    else if (ch <= TRAIL_SURROGATE_MAX)
        return 0; //unexpected trail surrogate!
    else
        return 1;
}


template <class CharIterator, class Function> inline
void utf16ToCodePoint(CharIterator first, CharIterator last, Function writeOutput) //"writeOutput" is a unary function taking a CodePoint
{
    static_assert(sizeof(typename std::iterator_traits<CharIterator>::value_type) == 2, "");

    for ( ; first != last; ++first)
    {
        CodePoint cp = static_cast<Char16>(*first);
        switch (getUtf16Len(static_cast<Char16>(cp)))
        {
            case 0: //invalid utf16 character
                cp = REPLACEMENT_CHAR;
                break;
            case 1:
                break;
            case 2:
                if (++first != last) //trail surrogate expected!
                {
                    const Char16 ch = static_cast<Char16>(*first);
                    if (TRAIL_SURROGATE <= ch && ch <= TRAIL_SURROGATE_MAX) //trail surrogate expected!
                    {
                        cp = ((cp - LEAD_SURROGATE) << 10) + (ch - TRAIL_SURROGATE) + 0x10000;
                        break;
                    }
                }
                --first;
                cp = REPLACEMENT_CHAR;
                break;
        }
        writeOutput(cp);
    }
}


template <class Function> inline
void codePointToUtf8(CodePoint cp, Function writeOutput) //"writeOutput" is a unary function taking a Char8
{
    //http://en.wikipedia.org/wiki/UTF-8
    //assert(cp < LEAD_SURROGATE || TRAIL_SURROGATE_MAX < cp); //code points [0xd800, 0xdfff] are reserved for UTF-16 and *should* not be encoded in UTF-8

    if (cp < 0x80)
        writeOutput(static_cast<Char8>(cp));
    else if (cp < 0x800)
    {
        writeOutput(static_cast<Char8>((cp >> 6  ) | 0xc0));
        writeOutput(static_cast<Char8>((cp & 0x3f) | 0x80));
    }
    else if (cp < 0x10000)
    {
        writeOutput(static_cast<Char8>( (cp >> 12       ) | 0xe0));
        writeOutput(static_cast<Char8>(((cp >> 6) & 0x3f) | 0x80));
        writeOutput(static_cast<Char8>( (cp & 0x3f      ) | 0x80));
    }
    else if (cp <= CODE_POINT_MAX)
    {
        writeOutput(static_cast<Char8>( (cp >> 18        ) | 0xf0));
        writeOutput(static_cast<Char8>(((cp >> 12) & 0x3f) | 0x80));
        writeOutput(static_cast<Char8>(((cp >> 6)  & 0x3f) | 0x80));
        writeOutput(static_cast<Char8>( (cp & 0x3f       ) | 0x80));
    }
    else //invalid code point
        codePointToUtf8(REPLACEMENT_CHAR, writeOutput); //resolves to 3-byte utf8
}


inline
size_t getUtf8Len(unsigned char ch) //ch must be first code unit! returns 0 on error!
{
    if (ch < 0x80)
        return 1;
    if (ch >> 5 == 0x6)
        return 2;
    if (ch >> 4 == 0xe)
        return 3;
    if (ch >> 3 == 0x1e)
        return 4;
    return 0; //innvalid begin of UTF8 encoding
}


template <class CharIterator> inline
bool decodeTrail(CharIterator& first, CharIterator last, CodePoint& cp) //decode trailing surrogate byte
{
    if (++first != last) //trail surrogate expected!
    {
        const Char8 ch = static_cast<Char8>(*first);
        if (ch >> 6 == 0x2) //trail surrogate expected!
        {
            cp = (cp << 6) + (ch & 0x3f);
            return true;
        }
    }
    --first;
    cp = REPLACEMENT_CHAR;
    return false;
}

template <class CharIterator, class Function> inline
void utf8ToCodePoint(CharIterator first, CharIterator last, Function writeOutput) //"writeOutput" is a unary function taking a CodePoint
{
    static_assert(sizeof(typename std::iterator_traits<CharIterator>::value_type) == 1, "");

    for ( ; first != last; ++first)
    {
        CodePoint cp = static_cast<Char8>(*first);
        switch (getUtf8Len(static_cast<Char8>(cp)))
        {
            case 0: //invalid utf8 character
                cp = REPLACEMENT_CHAR;
                break;
            case 1:
                break;
            case 2:
                cp &= 0x1f;
                decodeTrail(first, last, cp);
                break;
            case 3:
                cp &= 0xf;
                if (decodeTrail(first, last, cp))
                    decodeTrail(first, last, cp);
                break;
            case 4:
                cp &= 0x7;
                if (decodeTrail(first, last, cp))
                    if (decodeTrail(first, last, cp))
                        decodeTrail(first, last, cp);
                if (cp > CODE_POINT_MAX) cp = REPLACEMENT_CHAR;
                break;
        }
        writeOutput(cp);
    }
}


template <class CharString> inline
size_t unicodeLength(const CharString& str, char) //utf8
{
    using CharType = typename GetCharType<CharString>::Type;

    const CharType*       strFirst  = strBegin(str);
    const CharType* const strLast   = strFirst + strLength(str);

    size_t len = 0;
    while (strFirst < strLast) //[!]
    {
        ++len;
        size_t utf8len = getUtf8Len(*strFirst);
        if (utf8len == 0) ++utf8len; //invalid utf8 character
        strFirst += utf8len;
    }
    return len;
}


template <class WideString> inline
size_t unicodeLengthWide(const WideString& str, Int2Type<2>) //windows: utf16-wchar_t
{
    using CharType = typename GetCharType<WideString>::Type;

    const CharType*       strFirst = strBegin(str);
    const CharType* const strLast  = strFirst + strLength(str);

    size_t len = 0;
    while (strFirst < strLast) //[!]
    {
        ++len;
        size_t utf16len = getUtf16Len(*strFirst);
        if (utf16len == 0) ++utf16len; //invalid utf16 character
        strFirst += utf16len;
    }
    return len;
}


template <class WideString> inline
size_t unicodeLengthWide(const WideString& str, Int2Type<4>) //other OS: utf32-wchar_t
{
    return strLength(str);
}


template <class WideString> inline
size_t unicodeLength(const WideString& str, wchar_t)
{
    return unicodeLengthWide(str, Int2Type<sizeof(wchar_t)>());
}
}


template <class UtfString> inline
size_t unicodeLength(const UtfString& str) //return number of code points
{
    return implementation::unicodeLength(str, typename GetCharType<UtfString>::Type());
}


namespace implementation
{
template <class CharString> inline
size_t findUnicodePos(const CharString& str, size_t unicodePos, char) //utf8-char
{
    using CharType = typename GetCharType<CharString>::Type;

    const CharType* strFirst = strBegin(str);
    const size_t strLen = strLength(str);

    size_t utfPos = 0;
    while (unicodePos-- > 0)
    {
        if (utfPos >= strLen)
            return strLen;

        size_t utf8len = getUtf8Len(strFirst[utfPos]);
        if (utf8len == 0) ++utf8len; //invalid utf8 character
        utfPos += utf8len;
    }
    if (utfPos >= strLen)
        return strLen;
    return utfPos;
}


template <class WideString> inline
size_t findUnicodePosWide(const WideString& str, size_t unicodePos, Int2Type<2>) //windows: utf16-wchar_t
{
    using CharType = typename GetCharType<WideString>::Type;

    const CharType* strFirst = strBegin(str);
    const size_t strLen = strLength(str);

    size_t utfPos = 0;
    while (unicodePos-- > 0)
    {
        if (utfPos >= strLen)
            return strLen;

        size_t utf16len = getUtf16Len(strFirst[utfPos]);
        if (utf16len == 0) ++utf16len; //invalid utf16 character
        utfPos += utf16len;
    }
    if (utfPos >= strLen)
        return strLen;
    return utfPos;
}


template <class WideString> inline
size_t findUnicodePosWide(const WideString& str, size_t unicodePos, Int2Type<4>) //other OS: utf32-wchar_t
{
    return std::min(strLength(str), unicodePos);
}


template <class UtfString> inline
size_t findUnicodePos(const UtfString& str, size_t unicodePos, wchar_t)
{
    return findUnicodePosWide(str, unicodePos, Int2Type<sizeof(wchar_t)>());
}
}


template <class UtfString> inline
size_t findUnicodePos(const UtfString& str, size_t unicodePos) //return position of unicode char in UTF-encoded string
{
    return implementation::findUnicodePos(str, unicodePos, typename GetCharType<UtfString>::Type());
}

//-------------------------------------------------------------------------------------------

namespace implementation
{
template <class WideString, class CharString> inline
WideString utf8ToWide(const CharString& str, Int2Type<2>) //windows: convert utf8 to utf16-wchar_t
{
    WideString output;
    utf8ToCodePoint(strBegin(str), strBegin(str) + strLength(str),
    [&](CodePoint cp) { codePointToUtf16(cp, [&](Char16 c) { output += static_cast<wchar_t>(c); }); });
    return output;
}


template <class WideString, class CharString> inline
WideString utf8ToWide(const CharString& str, Int2Type<4>) //other OS: convert utf8 to utf32-wchar_t
{
    WideString output;
    utf8ToCodePoint(strBegin(str), strBegin(str) + strLength(str),
    [&](CodePoint cp) { output += static_cast<wchar_t>(cp); });
    return output;
}


template <class CharString, class WideString> inline
CharString wideToUtf8(const WideString& str, Int2Type<2>) //windows: convert utf16-wchar_t to utf8
{
    CharString output;
    utf16ToCodePoint(strBegin(str), strBegin(str) + strLength(str),
    [&](CodePoint cp) { codePointToUtf8(cp, [&](char c) { output += c; }); });
    return output;
}


template <class CharString, class WideString> inline
CharString wideToUtf8(const WideString& str, Int2Type<4>) //other OS: convert utf32-wchar_t to utf8
{
    CharString output;
    std::for_each(strBegin(str), strBegin(str) + strLength(str),
    [&](CodePoint cp) { codePointToUtf8(cp, [&](char c) { output += c; }); });
    return output;
}
}


template <class CharString> inline
bool isValidUtf8(const CharString& str)
{
    using namespace implementation;
    bool valid = true;
    utf8ToCodePoint(strBegin(str), strBegin(str) + strLength(str),
                    [&](CodePoint cp)
    {
        if (cp == REPLACEMENT_CHAR)
            valid = false; //perf: should we use an (expensive) exception for iteration break?
    });
    return valid;
}


template <class WideString, class CharString> inline
WideString utf8ToWide(const CharString& str)
{
    static_assert(IsSameType<typename GetCharType<CharString>::Type, char   >::value, "");
    static_assert(IsSameType<typename GetCharType<WideString>::Type, wchar_t>::value, "");

    return implementation::utf8ToWide<WideString>(str, Int2Type<sizeof(wchar_t)>());
}


template <class CharString, class WideString> inline
CharString wideToUtf8(const WideString& str)
{
    static_assert(IsSameType<typename GetCharType<CharString>::Type, char   >::value, "");
    static_assert(IsSameType<typename GetCharType<WideString>::Type, wchar_t>::value, "");

    return implementation::wideToUtf8<CharString>(str, Int2Type<sizeof(wchar_t)>());
}

//-------------------------------------------------------------------------------------------

template <class TargetString, class SourceString> inline
TargetString utfCvrtTo(const SourceString& str, char, wchar_t) { return utf8ToWide<TargetString>(str); }

template <class TargetString, class SourceString> inline
TargetString utfCvrtTo(const SourceString& str, wchar_t, char) { return wideToUtf8<TargetString>(str); }

template <class TargetString, class SourceString> inline
TargetString utfCvrtTo(const SourceString& str, char, char) { return copyStringTo<TargetString>(str); }

template <class TargetString, class SourceString> inline
TargetString utfCvrtTo(const SourceString& str, wchar_t, wchar_t) { return copyStringTo<TargetString>(str); }

template <class TargetString, class SourceString> inline
TargetString utfCvrtTo(const SourceString& str)
{
    return utfCvrtTo<TargetString>(str,
                                   typename GetCharType<SourceString>::Type(),
                                   typename GetCharType<TargetString>::Type());
}
}

#endif //UTF_H_01832479146991573473545
