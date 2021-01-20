---
title: Downloading ParSim
layout: single
author_profile: true
permalink: "docs/getting-started/"
sidebar:
  nav: "parsim_docs"
---

* Contents
{:toc}

### User version

* Create a clone of the ParSim repository

{% highlight sh %}
   git clone https://github.com/bbanerjee/ParSim
{% endhighlight %}

* Check status

{% highlight sh %}
  git status
  git diff origin/master
{% endhighlight %}

* Just in case the master has changed, update your copy

{% highlight sh %}
  git pull
{% endhighlight %}

### Developer version (for Ubuntu)

If you are a developer, we will suggest the following longer version.

* Create a `GitHub` account

{% highlight sh %}
   https://github.com/signup/free
{% endhighlight %}

* Install `git` on your machine

{% highlight sh %}
   sudo apt-get install git-core
{% endhighlight %}

* Setup git configuration

{% highlight sh %}
   git config --global user.email "your_email@your_address.com"
   git config --global user.name "your_github_username"
{% endhighlight %}

* Email your `github` username to Biswajit @ gmail.com

{% highlight sh %}
   b.banerjee.nz
{% endhighlight %}

* Create a clone of the ParSim repository

{% highlight sh %}
   git clone https://github.com/bbanerjee/ParSim
{% endhighlight %}

* Check branches/master
{% highlight sh %}
  git branch -a
{% endhighlight %}

* Configure `ssh` for secure access (if needed)
{% highlight sh %}
  cd ~/.ssh
  mkdir backup
  cp * backup/
  rm id_rsa*
  ssh-keygen -t rsa -C "your_email@your_address.com"
  sudo apt-get install xclip
  cd ~/ParSim/
  xclip -sel clip < ~/.ssh/id_rsa.pub
  ssh -T git@github.com
  ssh-add
  ssh -T git@github.com
{% endhighlight %}

* For password-less access with `ssh`, use
{% highlight sh %}
  git remote set-url origin git@github.com:bbanerjee/ParSim.git
{% endhighlight %}

* To change the default message editor, do
{% highlight sh %}
  git config --global core.editor "vim"
{% endhighlight %}

#### Getting the latest version, adding files etc.
* Check status

{% highlight sh %}
  git status
  git diff origin/master
{% endhighlight %}

* Update local repository

{% highlight sh %}
  git pull
{% endhighlight %}

* Add file to your local repository

{% highlight sh %}
  git add README
  git commit -m 'Added README file for git'
{% endhighlight %}

* Update main repository with your changes

{% highlight sh %}
  git push origin master
{% endhighlight %}



