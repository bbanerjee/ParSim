/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef _MATITI_HASH64_H
#define _MATITI_HASH64_H

#include <Types/Types.h>

namespace Matiti {

  // function object class for Hashing with lookup3
  struct Hash64 {

    std::size_t operator() (const long64& cellID) const {

      const u8* key = (const u8*) &cellID;
      u32 len = sizeof(cellID);
      u32 seed = 13;

      //return lookup3((const u8*) &key, sizeof(key), 13 );
    //}

    // lookup3 hash function
    //inline u32 lookup3( const u8 *key, u32 len, u32 seed ) {
      #if defined(_MSC_VER)
        #define rot(x,k) _rotl(x,k)
      #else
        #define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))
      #endif

      #define mix(a,b,c) \
      { \
        a -= c;  a ^= rot(c, 4);  c += b; \
        b -= a;  b ^= rot(a, 6);  a += c; \
        c -= b;  c ^= rot(b, 8);  b += a; \
        a -= c;  a ^= rot(c,16);  c += b; \
        b -= a;  b ^= rot(a,19);  a += c; \
        c -= b;  c ^= rot(b, 4);  b += a; \
      }

      #define final(a,b,c) \
      { \
        c ^= b; c -= rot(b,14); \
        a ^= c; a -= rot(c,11); \
        b ^= a; b -= rot(a,25); \
        c ^= b; c -= rot(b,16); \
        a ^= c; a -= rot(c,4);  \
        b ^= a; b -= rot(a,14); \
        c ^= b; c -= rot(b,24); \
      }

      u32 a, b, c;
      a = b = c = 0xdeadbeef + len + seed;

      const u32 *k = (const u32 *)key;

      while ( len > 12 ) {
        a += k[0];
        b += k[1];
        c += k[2];
        mix(a,b,c);
        len -= 12;
        k += 3;
      }

      switch( len ) {
        case 12: c+=k[2]; b+=k[1]; a+=k[0]; break;
        case 11: c+=k[2]&0xffffff; b+=k[1]; a+=k[0]; break;
        case 10: c+=k[2]&0xffff; b+=k[1]; a+=k[0]; break;
        case 9 : c+=k[2]&0xff; b+=k[1]; a+=k[0]; break;
        case 8 : b+=k[1]; a+=k[0]; break;
        case 7 : b+=k[1]&0xffffff; a+=k[0]; break;
        case 6 : b+=k[1]&0xffff; a+=k[0]; break;
        case 5 : b+=k[1]&0xff; a+=k[0]; break;
        case 4 : a+=k[0]; break;
        case 3 : a+=k[0]&0xffffff; break;
        case 2 : a+=k[0]&0xffff; break;
        case 1 : a+=k[0]&0xff; break;
        case 0 : return c;
      }

      final(a,b,c);
      return c;
    }

  };
}

#endif
