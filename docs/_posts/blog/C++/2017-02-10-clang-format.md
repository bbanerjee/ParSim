---
title:  "Formatting C++ code"
subheadline: "Biswajit Banerjee"
description: "Using clang-format and integration with vim"
date:  2017-02-10 09:30:00
categories:
    - C++
excerpt_separator: <!--more-->
toc: true
toc_label: "Contents"
toc_sticky: true
---

Every once in a while I need to recall how `clang-format` is set up and run on my Ubuntu 16.04 
machine.  This post is meant to act as a reminder.

<!--more-->
If I want to change all the files in my repository to a particular format (I prefer the
Mozilla style), I use:

~~~ bash

  find . -name "*.cc" -exec clang-format -i -style Mozilla {} \;
  find . -name "*.h" -exec clang-format -i -style Mozilla {} \;

~~~

For integration into vim, I create a `.vimrc` file and add the following into it

~~~ vim

  map <C-I> :py3f /path/to/clang-format.py<cr>
  imap <C-I> <c-o> :py3f /path/to/clang-format.py<cr>
  map <C-a> <esc>ggVG<cr>

~~~

The first two commands map `ctrl-i` to clang-format while the third command just
maps `ctrl-a` to "select all the code".

The `clang-format.py` file needs to be a Python3 compatible version and can be downloaded
from [LLVM](https://github.com/llvm-mirror/clang/blob/master/tools/clang-format/clang-format.py).

To make sure that the right format is used by vim, you have to create a 
`.clang-format` file containing the style options you will use.  In my case I use

~~~ bash

  clang-format -style=Mozilla -dump-config > .clang-format

~~~

You can also watch a talk on various clang tools below.

<iframe width='970' height='546' src='//www.youtube.com/embed/cX_GhJ6BuWI' frameborder='0' allowfullscreen></iframe>

