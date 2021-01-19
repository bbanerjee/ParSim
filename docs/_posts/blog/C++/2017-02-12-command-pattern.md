---
layout: posts
title:  "Command pattern for regression testing"
subheadline: "Biswajit Banerjee"
description: "Using the command pattern in C++"
date:  2017-02-12 09:30:00
categories:
    - C++
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---
A few days ago I had to refactor a 20,000 line class that was being used as
a regression tester for [discrete element](https://en.wikipedia.org/wiki/Discrete_element_method)
and [peridynamics](https://en.wikipedia.org/wiki/Peridynamics) simulations.  After some thought, I 
realized that the easiest way to achieve what I wanted was by using the [Command 
design pattern](https://en.wikipedia.org/wiki/Command_pattern) (minus the undo option).
My approach was to keep the implementation simple enough that a typical student of
mechanics would be able to follow the details, i.e., no [template metaprogramming](https://en.wikipedia.org/wiki/Template_metaprogramming) a la [Loki](http://loki-lib.sourceforge.net/index.php?n=Main.HomePage).

#### The original code ####
The regression tester in the original code looked something like this:

{% highlight cpp %}
    int main() {
      RegressionTester tester;
      switch(testType) {
      case TEST1:
        tester.doTest1();
        break;
      case TEST2:
        tester.doTest2();
        break;
      ...
      }
      return 0;
    }
{% endhighlight %}

The `RegressionTester` class was huge and contained the code for a large number of tests,
each of which was quite involved.

I wanted to separate out each test into its own self-contained class.

#### The Command pattern approach ####
There are several variations on the command pattern, but the basic idea is to 
use object polymorphism to provide a clean interface for function calls.

##### Step 1: Create a Command interface #####
I created a directory called `TestSimulations` and in that directory created the Command interface
file `Command.h`.  

{% highlight cpp %}
#include <RegressionTester.h>
class Command
{
public:
  virtual ~Command(){};
  virtual void execute(RegressionTester* tester) = 0;
};
{% endhighlight %}

##### Step 2: Create a Command handler #####
The command handler class returns a pointer to the chosen command that is polymorphic
and will return the pointer to the correct test. The header first

{% highlight cpp %}
#include <memory>
class Command;
using CommandP = std::unique_ptr<Command>;

class CommandHandler
{
public:
  CommandP handleCommand(std::string testType);
};
{% endhighlight %}

and then the implementation (`getEnum` translates the string to an enum)

{% highlight cpp %}
#include <TestSimulations/Command.h>
#include <TestSimulations/CommandHandler.h>
#include <TestSimulations/Test1.h>
#include <TestSimulations/Test2.h>
CommandP
CommandHandler::handleCommand(std::string testType)
{
  switch (getEnum(testType)) {
  case TEST1: 
    return std::make_unique<Test1>();
    break;
  case TEST2: 
    return std::make_unique<Test2>();
    break;
  ...
  }
  return nullptr;
}
{% endhighlight %}

##### Step 3: Create the test classes #####
Next we add the actual test classes.  The headers have the form

{% highlight cpp %}
#include <TestSimulations/Command.h>
class Test1 : public Command
{
public:
  virtual void execute(RegressionTester* tester);
};
{% endhighlight %}

while the code contains the detailed algorithm for each test:

{% highlight cpp %}
#include <TestSimulations/Test1.h>
void
Test1::execute(RegressionTester* tester)
{
  // Complex test algorithm
}
{% endhighlight %}

##### Step 4: Call the regression tester #####
This is the entry point.  Even though we have used the `RegressionTester`
class, we haven't explained what it does.  You can use it for generic 
algorithms that some or all tests use, or for anything else that each
test may need.
{% highlight cpp %}
#include <RegressionTester.h>
#include <TestSimulations/Command.h>
#include <TestSimulations/CommandHandler.h>
int main() {
  std::vector<std::string> testTypes = getTestTypesFromInput();
  RegressionTester tester;
  CommandHandler handler;
  // Loop through the tests 
  for (auto testType : testTypes) {
    // The command hadler returns the right type of test object
    CommandP command = handler.handleCommand(testType);
    // Run the test
    command->execute(&tester);
  }
}
{% endhighlight %}

##### Step 5: Run you regression tests #####
Make sure you save a set of "gold-standard" outputs to compare with the
results from your regression tester.  Typically a Python script or a shell
script is used to do the runs, comparisons, and output to a webdocs.

#### Conclusion ####
So, at the expense of an increase in the number of classes and the need for
some look-ups of the virtual table, we have a much cleaner implementation
of the tests and we can add more tests quite easily.


<a href="https://twitter.com/share" class="twitter-share-button" data-via="parresianz">Tweet</a>
<script>!function(d,s,id){var js,fjs=d.getElementsByTagName(s)[0],p=/^http:/.test(d.location)?'http':'https';if(!d.getElementById(id)){js=d.createElement(s);js.id=id;js.src=p+'://platform.twitter.com/widgets.js';fjs.parentNode.insertBefore(js,fjs);}}(docsument, 'script', 'twitter-wjs');</script>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

<!-- <script src="https://cdn.rawgit.com/google/code-prettify/master/loader/run_prettify.js?lang=cpp&amp;skin=sunburst"></script> -->
