---
title: Checking the build
layout: single
author_profile: true
sidebar:
  nav: "parsim_docs"
---
---

To check the build, you should run a few example problems.

* Create a directory called `runs` under `Vaango`

{% highlight sh %}
     mkdir runs
{% endhighlight %}

* Create links to the optimized and debug executables:

{% highlight sh %}
     cd runs
     ln -s ../dbg/StandAlone/vaango vaango_dbg
     ln -s ../opt/StandAlone/vaango vaango_dbg
{% endhighlight %}

* Create links to the `inputs` directory:

{% highlight sh %}
     ln -s ../src/StandAlone/inputs/MPM MPM_inputs
{% endhighlight %}

* Do a serial test run

{% highlight sh %}
     ./vaango_opt MPM_inputs/const_test_hypo.ups
{% endhighlight %}

* Do a parallel test run

{% highlight sh %}
    mpirun  -np 1 ./vaango_opt MPM_inputs/const_test_hypo.ups
{% endhighlight %}

