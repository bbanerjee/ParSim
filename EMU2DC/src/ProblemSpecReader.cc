#include <ProblemSpecReader.h>
#include <Exception.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>

#include <libxml/tree.h>
#include <libxml/parser.h>

using namespace Emu2DC;

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
	std::cout << "WARNING: Directory not returned by getcwd()...\n";
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


