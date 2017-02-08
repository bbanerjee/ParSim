// *****************************************************************************
// * This file is part of the FreeFileSync project. It is distributed under    *
// * GNU General Public License: http://www.gnu.org/licenses/gpl-3.0           *
// * Copyright (C) Zenju (zenju AT freefilesync DOT org) - All Rights Reserved *
// *****************************************************************************

#ifndef IO_H_8917640501480763248343343
#define IO_H_8917640501480763248343343

#include <cstdio>
#include <cerrno>
#include "scope_guard.h"
#include "utf.h"
#include "error.h"

namespace zen
{
/**
\file
\brief Save and load byte streams from files
*/

#if !defined(ZEN_WIN) && !defined(ZEN_LINUX) && !defined(ZEN_MAC)
    #error Please specify your platform: #define ZEN_WIN, ZEN_LINUX or ZEN_MAC
#endif


///Exception thrown due to failed file I/O
struct XmlFileError : public XmlError
{
    using ErrorCode = int;

    explicit XmlFileError(ErrorCode ec) : lastError(ec) {}
    ///Native error code: errno
    ErrorCode lastError;
};


#ifdef ZEN_WIN
namespace implemenation //sad but true
{
template <class String> inline
FILE* fopen(const String& filename, const wchar_t* mode)
{
#ifdef _MSC_VER
    FILE* handle = nullptr;
    errno_t rv = ::_wfopen_s(&handle, utfCvrtTo<std::wstring>(filename).c_str(), mode); //really more secure now?
    (void)rv;
    return handle;
#else
    return ::_wfopen(utfCvrtTo<std::wstring>(filename).c_str(), mode);
#endif
}
}
#endif


///Save byte stream to a file
/**
\tparam String Arbitrary string-like type: e.g. std::string, wchar_t*, char[], wchar_t, wxString, MyStringClass, ...
\param stream Input byte stream
\param filename Output file name
\throw XmlFileError
*/
template <class String>
void saveStream(const std::string& stream, const String& filename) //throw XmlFileError
{
#ifdef ZEN_WIN
    FILE* handle = implemenation::fopen(utfCvrtTo<std::wstring>(filename).c_str(), L"wb");
#else
    FILE* handle = ::fopen(utfCvrtTo<std::string>(filename).c_str(), "w");
#endif
    if (!handle)
        throw XmlFileError(errno);
    ZEN_ON_SCOPE_EXIT(::fclose(handle));

    const size_t bytesWritten = ::fwrite(stream.c_str(), 1, stream.size(), handle);
    if (::ferror(handle) != 0)
        throw XmlFileError(errno);

    (void)bytesWritten;
    assert(bytesWritten == stream.size());
}


///Load byte stream from a file
/**
\tparam String Arbitrary string-like type: e.g. std::string, wchar_t*, char[], wchar_t, wxString, MyStringClass, ...
\param filename Input file name
\return Output byte stream
\throw XmlFileError
*/
template <class String>
std::string loadStream(const String& filename) //throw XmlFileError
{
#ifdef ZEN_WIN
    FILE* handle = implemenation::fopen(utfCvrtTo<std::wstring>(filename).c_str(), L"rb");
#else
    FILE* handle = ::fopen(utfCvrtTo<std::string>(filename).c_str(), "r");
#endif
    if (!handle)
        throw XmlFileError(errno);
    ZEN_ON_SCOPE_EXIT(::fclose(handle));

    std::string stream;
    const size_t blockSize = 64 * 1024;
    do
    {
        stream.resize(stream.size() + blockSize); //let's pray std::string implements exponential growth!

        const size_t bytesRead = ::fread(&*(stream.end() - blockSize), 1, blockSize, handle);
        if (::ferror(handle))
            throw XmlFileError(errno);
        if (bytesRead > blockSize)
            throw XmlFileError(0);
        if (bytesRead < blockSize)
            stream.resize(stream.size() - blockSize + bytesRead); //caveat: unsigned arithmetics
    }
    while (!::feof(handle));

    return stream;
}
}

#endif //IO_H_8917640501480763248343343
