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

#include <InputOutput/ProblemSpecReader.h>
#include <Core/Exception.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>

#include <libxml/tree.h>
#include <libxml/parser.h>

using namespace Matiti;

ProblemSpecReader::ProblemSpecReader()
{
}

ProblemSpecReader::~ProblemSpecReader()
{
  d_xmlData = 0;

  for( unsigned int pos = 0; pos < d_upsFilename.size(); pos++ ) {
    delete d_upsFilename[ pos ];
  }
}

std::string
getPath( const std::string & filename )
{
  return filename.substr( 0, filename.rfind( "/" ) );
}

std::string
validateFilename( const std::string & filename, const xmlNode * parent )
{
  std::string fullFilename;
  std::string errorMsg;

  if( filename[0] != '/') { // If not absolute path, make it one...
          
    if( parent ) {
      fullFilename = getPath( *((std::string*)(parent->_private)) ) + "/" + filename;
      std::cout << "1) filename: " << fullFilename << "\n";
    }

    if( !parent ) { // Check to see if the file is relative to where the program was run... 

      char buffer[2000];
      char * str = getcwd( buffer, 2000 );
      if( str == NULL ) {
        throw Exception("**ERROR** Directory not returned by getcwd()", __FILE__, __LINE__);
      }
      else {
        fullFilename = std::string(buffer) + "/" + filename;
      }
    }
  }
  else {
    fullFilename = filename;
  }

  // Check file name open/close
  std::ifstream infile;
  infile.open(fullFilename.c_str());
  if (infile) {
    infile.close();
  } else {
    throw Exception("Invalid input file name", __FILE__, __LINE__);
  }

  return fullFilename;

} // end validateFilename()


Uintah::ProblemSpecP
ProblemSpecReader::readInputFile(const std::string & filename)
{
  if( d_xmlData != 0 ) {
    return d_xmlData;
  }

  static bool initialized = false;
  if (!initialized) {
    LIBXML_TEST_VERSION;
    initialized = true;
  }

  std::string full_filename = validateFilename( filename, NULL );

  xmlDocPtr doc = xmlReadFile( full_filename.c_str(), 0, XML_PARSE_PEDANTIC );
    
  // you must free doc when you are done.
  // Add the parser contents to the ProblemSpecP
  Uintah::ProblemSpecP prob_spec = new Uintah::ProblemSpec( xmlDocGetRootElement(doc), true );

  std::string * strPtr = new std::string( full_filename );

  d_upsFilename.push_back( strPtr );
  prob_spec->getNode()->_private = (void*)strPtr;

  d_xmlData = prob_spec;

  return prob_spec;

} // end readInputFile()


