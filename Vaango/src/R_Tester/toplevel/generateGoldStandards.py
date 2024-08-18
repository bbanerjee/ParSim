#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2024 Biswajit Banerjee
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.
#

#! /usr/bin/env python
import os
import shutil
import platform
import socket
from optparse import OptionParser
from sys import argv, exit
from string import upper


# bulletproofing
if os.sys.version_info <= (2,5):
  print ""
  print "ERROR: Your python version [" + str( os.sys.version_info ) + "] is too old.\n" + \
        "       You must use version 2.5 or greater (but NOT version 3.x!). \n" + \
        "       If you're using either ember or updraft please add the following to your shell rc script: \n" + \
        "           csh/tcsh:  source  /uufs/chpc.utah.edu/sys/pkg/python/2.6.5/etc/python.csh \n" + \
        "           bash:      source /uufs/chpc.utah.edu/sys/pkg/python/2.6.5/etc/python.sh\n\n"
  print ""
  exit( 1 )

import subprocess

from helpers.runVaangoTests import nameoftest, testOS, input, num_processes, testOS, setGeneratingGoldStandards, userFlags

####################################################################################

vaango        = ""   # full path to vaango executable
inputs        = ""   # full path to src/Standalone/inputs/
OS            = platform.system()
debug_build   = ""
no_sci_malloc = ""
                      # HACK 
                      # 1 for GPU RT machine (albion), 0 otherwise.
                      #   need to make this generic, perhaps pycuda?
has_gpu       = 1 if socket.gethostname() == "albion" else 0

####################################################################################

#script_dir=os.sys.path[0]

usage = "%prog [options]\n\n" \
        "   generateGoldStandards creates a sub-directory for the component test_file.\n" \
        "   Note, multiple tests may be specified: -t ICE -t MPM etc."

parser = OptionParser( usage, add_help_option=False )

parser.set_defaults( verbose=False, 
                     parallelism=1 )

parser.add_option( "-b", dest="build_directory",         help="Vaango build directory [REQUIRED]",
                   action="store", type="string" )

parser.add_option( "-d", dest="is_debug",                help="Whether this is a debug build (use 'yes' or 'no')",
                   action="store", type="string" )

parser.add_option( "-h", "--help", action="help",        help="Show this help message" )

parser.add_option( "-j", type="int", dest="parallelism", help="Set make parallelism" )

parser.add_option( "-m", dest="sci_malloc_on",           help="Whether this is build has sci-malloc turned on (use 'yes' or 'no')",
                   action="store", type="string" )

parser.add_option( "-s", dest="src_directory",           help="Vaango src directory [defaults to .../bin/../src]",
                   action="store", type="string" )

parser.add_option( "-t", dest="component_test_file",               help="Name of specific test script (eg: ICE) [REQUIRED/Multiple allowed]",
                   action="append", type="string" )

parser.add_option( "-v", action="store_true", dest="verbose", help="Enable verbosity" )

parser.add_option( "-o", dest="output_directory",        help="Goldstandaard directory [defaults to build directory]",
                   action="store", type="string" )

####################################################################################

def error( error_msg ) :
    print ""
    print "ERROR: " + error_msg
    print ""
    parser.print_help()
    print ""
    exit( 1 )

####################################################################################

def validateArgs( options, args ) :
    global vaango, inputs, debug_build, no_sci_malloc

    if len( args ) > 0 :
        error( "Unknown command line args: " + str( args ) )

    if not options.build_directory :
        error( "Vaango build directory is required..." )
    elif options.build_directory[0] != "/" : 
        error( "Vaango build directory must be an absolute path (ie, it must start with '/')." )
    elif options.build_directory[-1] == "/" : 
        # Cut off the trailing '/'
        options.build_directory = options.build_directory[0:-1]

    if not options.component_test_file :
        error( "A component test file must be specified..." )

    if not os.path.isdir( options.build_directory ) :
        error( "Build directory '" + options.build_directory + "' does not exist." )

    vaango = options.build_directory + "/StandAlone/vaango"

    if not os.path.isfile( vaango ) :
        error( "'vaango' not here: '" + vaango + "'" )

    if not options.sci_malloc_on :
        error( "Whether this is build has sci-malloc turned on is not specified.  Please use <-m yes/no>." )
    else :
        if options.sci_malloc_on == "yes" :
            no_sci_malloc = False
        elif options.sci_malloc_on == "no" :
            no_sci_malloc = True
        else :
            error( "-d requires 'yes' or 'no'." )

    if not options.is_debug :
        error( "debug/optimized not specified.  Please use <-d yes/no>." )
    else :
        if options.is_debug == "yes" :
            debug_build = True
        elif options.is_debug == "no" :
            debug_build = False
        else :
            error( "-d requires 'yes' or 'no'." )

    if not options.src_directory :
        # Cut off the <bin> and replace it with 'src'
        last_slash = options.build_directory.rfind( "/" )
        options.src_directory = options.build_directory[0:last_slash] + "/src"

    if not os.path.isdir( options.src_directory ) :
        error( "Src directory '" + options.src_directory + "' does not exist." )

    inputs = options.src_directory + "/StandAlone/inputs"

    if not os.path.isdir( inputs ) :
        error( "'inputs' directory not found here: '" + inputs )

    if not options.output_directory :
       # Set output directory to build directory
       options.output_directory = os.path.join(options.build_directory, "gold_standards")

    setGeneratingGoldStandards( inputs )

####################################################################################

def generateGS() :

    global vaango, inputs, debug_build, no_sci_malloc, has_gpu
    try :
        (options, leftover_args ) = parser.parse_args()
    except :
        print "" # Print an extra newline at end of output for clarity
        exit( 1 )

    validateArgs( options, leftover_args )
    
    #__________________________________
    # Does mpirun command exist or has the environmental variable been set?
    try :
      MPIRUN = os.environ['MPIRUN']    # first try the environmental variable
    except :
      MPIRUN = "mpirun"
      rc = os.system("which mpirun > /dev/null 2>&1 ")

      if rc == 256:
        print "ERROR:generateGoldStandards.py "
        print "      mpirun command was not found and the environmental variable MPIRUN was not set."
        print "      You must either add mpirun to your path, or set the 'MPIRUN' environment variable."
        exit (1)
    print "Using mpirun: %s " % MPIRUN
        
    if options.verbose :
        #print "Building Gold Standards in " + os.getcwd()
        print "Building Gold Standards in " + options.output_directory

    ##############################################################
    # Determine if the code has been modified (svn stat)
    ##############################################################
    #process = subprocess.Popen( "git log HEAD..origin/master --oneline", shell=True, 
    #                            stdout=subprocess.PIPE, stderr=subprocess.PIPE )
    #( stdout, sterr ) = process.communicate()
    #result = process.returncode

    #if result != 0 :
    #    answer = ""
    #    while answer != "n" and answer != "y" :
    #        print ""
    #        print "WARNING:  'git log' failed to run correctly, so generateGoldStandards.py cannot tell"
    #        print "          if your tree is 'up to date'.  Are you sure you want to continue generating"
    #        print "          new gold standards at this time? [y/n]"
    #        print ""
 
    #         answer = os.sys.stdin.readline()[:-1]
    #        if answer == "n" :
    #            print ""
    #            print "Goodbye."
    #            print ""
    #            exit( 0 )
    #
    # !!!FIXME!!!: if svn fails to run, or returns differences, then
    # ask the user what they want to do.  This is not complete yet!!!
    ##############################################################

    ##############################################################
    # !!!FIXME!!!: Determine if configVars has changed... if not, let the user
    # know this...
    #
    #configVars = options.build_directory + "/configVars.mk"
    #shutil.copy( configVars, "." )
    ##############################################################

    ##############################################################
    # !!!FIXME!!!:
    # - Determine if/where mpirun is...
    # - Determine if vaango was built with MPI...
    # - Determine (ask the user?) if the (vaango) binary is up to date.
    ##############################################################

    components = options.component_test_file
    
    # Now change directory to output directory
    os.chdir( options.output_directory )

    print "Current directory: " + os.getcwd()

    # Exit if the component hasn't been compiled.  Note, not all components
    # are listed in the configVars.mk file 
    # configVars = options.build_directory + "/configVars.mk"
    # for component in components :
    #   
    #   searchString = "BUILD_%s=no" % upper(component)  # search for BUILD_<COMPONENT>=no
    #   for line in open(configVars):
    #     if searchString in line:
    #       print "\n ERROR: the component (%s) was not compiled.  You must compile it before you can generate the gold standards\n" % component
    #       exit( 1 ) 

    # Warn user if directories already exist
    some_dirs_already_exist = False

    for component in components :
        if os.path.isdir( component ) :
            if not some_dirs_already_exist :
                some_dirs_already_exist = True
                print ""
                print "Note, the following gold standards already exist: ",
            else :
                print ", ",
            os.sys.stdout.write( component )

    if some_dirs_already_exist :
        answer = ""
        while answer != "n" and answer != "y" :
            print ""
            print "Delete existing gold standards?  (If 'no', script will exit.) [y/n]"
            answer = os.sys.stdin.readline()[:-1]
            if answer == "n" :
                print ""
                print "Goodbye."
                print ""
                exit( 0 )

        for component in components :
            if os.path.isdir( component ) :
                print "Deleting " + component
                shutil.rmtree( component )

    #print components
    #print some_dirs_already_exist
    #print "Press enter to continue."
    #wait = os.sys.stdin.readline()[:-1]
    
    # Loop through component test files
    for component in components :

        # Pull the list of tests from the the 'component's python module's 'TESTS' variable:
        # (Need to 'import' the module first.)
        if options.verbose :
            print "Python importing " + component + ".py"

        try :
          THE_COMPONENT = __import__( component )
        except :
          print ""
          print "Error: loading the component '%s'." % component
          print "       Either that python file does not exist or there is a syntax error in the tests that have been defined.  Goodbye."
          print ""
          exit( -1 )

        os.mkdir( component )
        os.chdir( component )

        # Create a symbolic link to the 'inputs' directory so some .ups files will be able
        # to find what they need...
        if not os.path.islink( "inputs" ) :
            os.symlink( inputs, "inputs" )


        # determine which tests (local/nightly) to run default is local
        whichTests = os.getenv( 'WHICH_TESTS', "local" )
        
        print "Which_tests: %s " % whichTests
        if whichTests == "local" :
            tests = THE_COMPONENT.getLocalTests()
        elif whichTests == "nightly" :
            tests = THE_COMPONENT.getNightlyTests()
        else :
            print "\nThe environmental variable WHICH_TESTS:(%s) is not valid" % whichTests
            print "the valid options are local or nightly. \n"
            exit (-1)
          
                  
        if options.verbose :
            print "About to run tests for component: " + component

        for test in tests :
            if testOS( test ) != upper( OS ) and testOS( test ) != "ALL":
                continue
             
            #  Defaults
            vaango_options     = ""
            do_gpu          = 0    # run test if gpu is supported 
            
            #__________________________________
            # parse user flags for the gpu and vaango_options
            # override defaults if the flags have been specified
            if len(test) == 5:
              flags = userFlags(test)
              print "User Flags:"

              #  parse the user flags
              for i in range(len(flags)):
                if flags[i] == "gpu":
                  do_gpu = 1 

                tmp = flags[i].rsplit('=')
                if tmp[0] == "vaango_options":
                  vaango_options = tmp[1]
                  print "\n vaango_option: %s \n"%(vaango_options)
    
            if do_gpu == 1 and has_gpu == 0:
              print "\nWARNING: skipping this test.  This machine is not configured to run gpu tests\n"
              continue
              
            # FIXME: NOT SURE IF THIS IS RIGHT, BUT IT APPEARS TO MATCH WHAT THE RUN TESTS SCRIPT NEEDS:
            print "About to run test: " + nameoftest( test )
            os.mkdir( nameoftest( test ) )
            os.chdir( nameoftest( test ) )

            # Create (yet) another symbolic link to the 'inputs' directory so some .ups files will be able
            # to find what they need...  (Needed for, at least, methane8patch (ARCHES) test.)
            if not os.path.islink( "inputs" ) :
                os.symlink( inputs, "inputs" )

            np = float( num_processes( test ) )
            mpirun = ""

            
            MALLOC_FLAG = ""

            if debug_build :
                if no_sci_malloc :
                    print ""
                    print "WARNING!!! The build was not built with SCI Malloc on...  Memory tests will not be run."
                    print "WARNING!!! If you wish to perform memory checks, you must re-configure your debug build"
                    print "WARNING!!! with '--enable-sci-malloc', run 'make cleanreally', and re-compile everything."
                    print ""
                else :
                    os.environ['MALLOC_STRICT'] = "set"
                    os.environ['MALLOC_STATS'] = "malloc_stats"
                    MALLOC_FLAG = " -x MALLOC_STATS "

            #SVN_FLAGS = " -svnStat -svnDiff "
            SVN_FLAGS = "" # When debugging, if you don't want to spend time waiting for SVN, uncomment this line.

            if np > 1.0 :
                np = int( np )
                mpirun = "%s -np %s  " % (MPIRUN,np)

                command = mpirun + MALLOC_FLAG + vaango + SVN_FLAGS + " " + vaango_options + " " + inputs + "/" + component + "/" + input( test )  #+ " >> vaango_log.txt " 
            else :
                command = vaango + SVN_FLAGS + " " + vaango_options + " " + inputs + "/" + component + "/" + input( test )  #+ " >> vaango_log.txt " 

            print "Running command: " + command

            rc = os.system( command )
            
            # catch if vaango doesn't run to completion
            if rc != 0:
              print "\nERROR: %s: Test (%s) failed to complete\n" % (component,test)
              exit(-1)
            
            os.chdir( ".." ) # Back to the component (eg: 'ICE') directory

        os.chdir( ".." ) # Back to the TestData directory


####################################################################################

if __name__ == "__main__":
    generateGS()
    exit( 0 )
