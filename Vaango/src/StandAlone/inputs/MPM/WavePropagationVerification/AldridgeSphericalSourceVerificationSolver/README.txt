File "README.txt". 

This file contains:
1. Instructions for algorithm execution.
2. Syntax and format of input and output.

Information about background and theory resides in a different file
named "INFO.txt".


-----------
Quick start:
-----------
To compile and execute the sample program, go to a UNIX command tool, enter
directory "AldridgeSphericalSource", and type the command "make".  That's all 
there is to it!  Provided that your input files (*.in) are the same as the
example inputs (*.in-example), the code should compile and run with the
results being the same as explod_for.out-example. See below for more details.


-------------------
System Requirements
-------------------

The algorithm is encoded as a simple Fortran 77 program, named EXPLOD_FOR, and 
is suitable for execution in a Unix/Linux workstation or laptop computational 
environment.  Very minimal input is required:  
 
      a) A text file containing a few algorithm control parameters and a 
         description of the data acquisition geometry (i.e., positions of the 
         source sphere and the receivers).
      b) A text file containing the source activation waveform (or what we 
         call the "source wavelet").
 
Algorithm output consists of an ascii text file containing the calculated 
traces (i.e., time series of elastic/acoustic particle displacement, velocity, 
acceleration, or pressure) at the various receiver locations.  Since the 
algorithm uses closed-form mathematical formulae to calculate responses, 
execution is virtually instantaneous.  Finally, computational memory demand 
is quite small:  there is a single 1D array of real*4 words (currently 
dimensioned at 1,000,001 elements).
 

--------- 
Key files
--------- 
 
1) File "explod_for.f":  Fortran 77 source code.
   The algorithm is entirely self-contained within this source code file.
   All called subroutines follow the main program code.  The prolog comments
   at the start of this file describe input parameters, output trace format, 
   memory allocation issues, compilation command, etc.
 
2) makefile: Associates Fortran source code input and output file numbers with 
   Unix file pathnames, compiles the executable (if necessary), and runs the
   code. If desired, edit the makefile to reflect changes in your input/output
   file names.
 
3) File "explod_for.in":  File containing example set of input parameters.
   Each numerical value of a parameter is followed by a variable name (used 
   within the source code).   There is also a section of commentary describing 
   the input parameters, allowed values, units of measure, etc.

4) File "heaviside.in": Example source waveform consisting of a Heaviside unit
   step function.  See "input syntax" below for details.


-------------------------------------------------
makefile (Execution instructions and other tasks)
-------------------------------------------------
The makefile is an alternative to a script. A makefile is better than a script
because it tracks dependencies. The makefile will ensure, for example, that
attempting to run the executable will be preceded by actually creating the
executable if it does not exist or if the sourcecode has been modified more
recently than the last creation of the executable.

Syntax from the command tool is  
        make <something>

For this makefile, <something> stands for any of the following:
	run
	explod_for.exe
	wipe
	release
Typing make with no argument is equivalent to "make run", which will run
the code.  A small amount of diagnostic info (mainly input parameters) is
echoed to standard output and will appear on the monitor screen.

According to the makefile, running the code (i.e., "make run") requires the 
following files to already exist: 
	explod_for.exe 
	explod_for.in 
	heaviside.in
The last two of these are input files, which means that the user must write
them.  However, the first one is the executable file, and the makefile
contains instructions for making the executable if it does not already 
exist or if its sourcecode file is "younger" than the executable (i.e., if you
modify the sourcecode, then the makefile automatically knows that a new
executable needs to be built).  Because makefiles take care of this dependency 
tracking, you will almost never need to explicitly create a new executable
(i.e., you will never need to type type "make explod_for.exe"). The makefile
knows to do this for you whenever sourcecode is modified.

Typing "make wipe" will delete files that are generated from a run, thereby
keeping your working directory neat and tidy.

Typing "make release" will tar together all files that someone else would
need to run this code themselves (using the makefile to create a release
is a safe approach because then you won't accidentally forget some key
file).

When you type "make" in a command tool, the "makefile" will check if any
changes have been made to the source code (explod_for.f); if so, the makefile 
will re-compile the executable explod_for.exe. For legacy reasons, the 
"makefile" creates symbolic links (shortcuts) from the key input/output files 
to files with much more cryptic names, fort.# (used in Fortran READ or WRITE 
statements).  Finally, the makefile executes the program.

Files input to EXPLOD_FOR:
   1) fort.1: explod_for.in:  algorithm input parameters (you edit to suit).
   2) fort.2: heaviside.out:  sample source waveform (in this case, a
                              Heaviside unit step function, but you write
                              your own).

Files output from EXPLOD_FOR:
   1) fort.11: explod_for.out: traces (time-series) containing
                               calculated responses (i.e., particle
                               displacement, velocity, acceleration, etc.
                               at receivers).


 
---------------------
Syntax of input files
---------------------
Input syntax will be described for the two sample input files that came with 
this release. If you want these input files (or output files) to have a
a different name, change their names in the makefile.

explod_for.in: For the sample input file that came with this release, the 
               pressure is recorded at 4 receivers arrayed along a line at 
               distances of 5, 10, 15, and 20 m from the center of a source 
               sphere of radius 0.3079 m (~1 foot).  The source is "radial 
               traction" with magnitude 10 million Pascals.  Output trace 
               simulation duration is 200 milliseconds and time sample interval 
               is 5 microseconds.  

heaviside.in:  This sample file specifies a heaviside waveform (i.e., an
               instantaneous jump from 0 to 1).  In general, a waveform
               file contains two columns of numbers, written in simple ascii 
               text format:

                  sample time                 sample amplitude

               For this particular source wavelet, which is a Heaviside 
               function, the sample amplitudes are all equal to 1.0.  There 
               are 40,001 lines in the file.  A common problem occurs when the 
               file containing the source waveform has fewer lines (or samples)
               than the requested  number that is determined indirectly from 
               parameters dt and tlen in the input run parameter file 
               "explod_for.in". Search the sourcecode for "nsamps" to see the 
               formulas used to determine the required length of the input 
               waveform.  If in doubt, the required length will be echoed to 
               the screen using the phrase "Number of trace samples requested"


---------------------
Format of output files
---------------------
Traces (time series) for the receivers, written in the same order as the 
order on input, consist of three columns of numbers:
 
  sample time        spatial coordinate        sample amplitude
 
Any "ringyness" at the onset of the waveforms is a Gibbs phenomenon effect, 
associated with the frequency domain calculations.  
 
 
Individual traces are separated by a ">" character.  The spatial coordinate is 
up to the user; it could be the x, y, or z coordinate of the receiver, or the
receiver number, or no coordinate at all (in which case there are only 2 
columns of numbers).  In fact, you should be able to enter the Fortran source 
code, locate the lines that write the output (toward the end of the main
program block), and tailor the output format to your specific desires.  
Similarly, the point in the code where the source wavelet is read can be
located, and the code modified to suit your needs.
 


RECALL: More detailed information about the background of this project and the
underlying physics or problem definition may be found in file "INFO.txt".
 

 
David F. Aldridge 
Geophysics Department 
Sandia National Laboratories 
Albuquerque, New Mexico, USA, 87185-0750 
telephone:   505-284-2823 
facsimile:   505-844-7354 
electronic mail:  dfaldri@sandia.gov 
