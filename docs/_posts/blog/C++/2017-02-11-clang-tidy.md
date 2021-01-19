---
layout: posts
title:  "Auto-modernizing C++ code"
subheadline: "Biswajit Banerjee"
description: "Using clang-tidy with cmake"
date:  2017-02-11 09:30:00
categories:
    - C++
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---
The Clang static analyzer tools come with a handy interface called [clang-tidy](http://releases.llvm.org/3.8.0/tools/clang/tools/extra/docs/clang-tidy/index.html).  

I've been try to set up this tool in my `cmake` toolchain with varying levels of success.  In 
particular I'd like to automate the conversion of old code into a version that uses C++11 and
some C++14 constructs.

---------

#### Attempt 1 ####

---------

In my first attempt I created a `cmake` directory at the top level of my source tree and
added a `clang-tidy.cmake` file that contained the following:

~~~bash

  file(GLOB_RECURSE ALL_MY_SOURCE_FILES *.cpp *.hpp *.cc *.h)
  add_custom_target(
        clang-tidy
        COMMAND /usr/bin/clang-tidy
        ${ALL_MY_SOURCE_FILES}
        -config=''
        -export-fixes='clang-tidy-fixes.dat'
        --
        -std=c++14
        -I${MY_SOURCE_INCLUDE_DIR}
        -I${MY_MPI_INCLUDE_DIR1}
        -I${MY_MPI_INCLUDE_DIR2}
        -I${MY_MPI_INCLUDE_DIR3}
        -I${MY_BOOST_INCLUDE_DIR}
        -I${MY_VTK_INCLUDE_DIR}
  )

~~~

At the bottom of my main `CMakeList.txt` file I added 

~~~bash

  include(cmake/clang-tidy.cmake)

~~~

When I ran `cmake` and then `make clang-tidy`, there were several complaints about
header files not being found, but the suggested fixes were written to the
file `clang-tidy-fixes.dat`.  That file contained numerous suggested fixes in the format

~~~bash

    - FilePath:        /path/to/file/parser.h
      Offset:          9558
      Length:          52
      ReplacementText: '(const auto & token : tokens)'

~~~

Going through that file and then hand-coding the changes
looked daunting, so I looked at automating the code change process.

---------

#### Attempt 2 ####

---------

To make `clang-tidy` apply the fixes automatically, I replaced

~~~ bash

        -export-fixes='clang-tidy-fixes.dat'

~~~

in `clang-tidy.cmake` to 

~~~ bash

        -fix-errors

~~~

When I ran `make clang-tidy` after that change, fixes were automatically applied,
but also applied multiple times.  For example, I got

~~~ c

   std::move(std::move(std::move(std::move(a))));

~~~

instead of

~~~ c

   std::move(a);

~~~

---------

#### Attempt 3 ####

---------
Clearly, attempt 2 wasn't solving the problem, and I went back to the `clang-tidy`
manual and found that I could actually export the compile commands from by build using

~~~ bash

  cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ../src

~~~

That command created a `compile_commands.json` file which contained all the 
path and header information that a build could need.

I then located the python script for running clang-tidy from

~~~ bash

  /usr/lib/llvm-3.8/share/clang/run-clang-tidy.py

~~~

and, from my build directory containing `compile_commands.json`, I ran

~~~ bash

  run-clang-tidy.py ../src -checks=modernize* -fix

~~~


No more complaints about missing header files, and all fixes appeared to have been
applied only once.  It's much easier to use this process than to try to create a
custom target for cmake. 



<a href="https://twitter.com/share" class="twitter-share-button" data-via="parresianz">Tweet</a>
<script>!function(d,s,id){var js,fjs=d.getElementsByTagName(s)[0],p=/^http:/.test(d.location)?'http':'https';if(!d.getElementById(id)){js=d.createElement(s);js.id=id;js.src=p+'://platform.twitter.com/widgets.js';fjs.parentNode.insertBefore(js,fjs);}}(docsument, 'script', 'twitter-wjs');</script>

<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>
