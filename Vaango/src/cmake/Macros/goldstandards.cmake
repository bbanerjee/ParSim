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

message("SRC dir = ${Vaango_SOURCE_DIR}")
message("BIN dir = ${Vaango_BINARY_DIR}")
find_package(PythonInterp REQUIRED)
file(MAKE_DIRECTORY ${Vaango_BINARY_DIR}/TestData)
set(components "${Vaango_SOURCE_DIR}/R_Tester/helpers/selectComponents.sh ${Vaango_SOURCE_DIR}/R_Tester")
set(IS_DEBUG "no")
set(SCI_MALLOC_ON "yes")
set(MAKE_PARALLELISM 4)
set(python_command "PYTHONPATH=${Vaango_SOURCE_DIR}/R_Tester/toplevel;${Vaango_SOURCE_DIR}/R_Tester python ${Vaango_SOURCE_DIR}/R_Tester/toplevel/generateGoldStandards.py -b ${Vaango_BINARY_DIR} -s ${Vaango_SOURCE_DIR} -t MPM -v")
message("CMD = ${python_command}")
execute_process(COMMAND ${python_command} RESULT_VARIABLE CMD_RESULT)
#if(CMD_RESULT)
#  message(FATAL_ERROR "Error running gold standards")
#endif()

