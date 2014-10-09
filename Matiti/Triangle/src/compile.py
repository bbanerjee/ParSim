#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#!/usr/bin/python

import os,sys

OFLAGS='OFLAGS=-g -Wall'
CC=''

if os.system('which icpc >/dev/null 2>&1') == 0:
    CC = 'icpc\n'
    OFLAGS = '\nOFLAGS=-O\n'

if os.system('which g++ >/dev/null 2>&1') == 0:
    CC = 'g++\n'
    OFLAGS = '\nOFLAGS=-O3\n'

print ' Compiler Selected      ... ',CC,
CC = 'CC='+CC+'\n'

TFLAGS = 'TFLAGS=-DREDUCED -DANSI_DECLARATORS -DTRILIBRARY -DCDT_ONLY -DLINUX' 

if sys.platform == 'cygwin':
    TFLAGS += ' -DCYGWIN'

    

L = open("./src/makefile.input").readlines()
L.insert(0, '\n#--------------------------\n')
L.insert(0, TFLAGS)
L.insert(0, OFLAGS)
L.insert(0, CC)
L.insert(0, '\n#--------------------------\n')
f = open("./src/makefile.tmp",'w')
f.writelines(L)
f.close()

print ' Compiling Triangle++   ...',
os.system("make --quiet -f ./src/makefile.tmp; rm ./src/makefile.tmp")
print 'Done.'

