# Dependencies

Building NEMO requires following minimum dependencies:

* Unix tools, such as: make, csh, bash, git
* C/C++/Fortran compiler
* pgplot (optional, but highly recommended)

# git

There are two ways you could have obtained the NEMO source code:

1. You cloned from the *upstream* : https://github.com/teuben/nemo . This is fine
if you just want to compile and run, but not ideal if you make modifications and
want to shared them back to the *upstream* via a *pull request* (PR). You
can however "repair" your local repo, discussed below, and still submit a PR.

2. You (forked)[https://guides.github.com/activities/forking/] NEMO
from the *upstream*, and cloned it locally from the repo in your own
github account. This is the ideal method, but you will still need to
set the *upstream* manually if you used github.com. See also the
**gh** command below for an even better way.

3. Sadly on github.com you will also find a **zip** copy of the repo
that does actually work fine, except it's a frozen snapshot and cannot
be efficiently used to collaborate. However, if you cannot install git,
this is probably the only way to bootstrap yourself. For example
https://github.com/teuben/nemo/archive/refs/heads/master.zip, which
will create a directory *nemo-master*. Other branches are available
through similar zip file.

Familiarize yourself with the concept of a pull request on github. There
are some links at the bottom of this document.


## gh:   github CLI

You can safely skip this section if you prefer to work via github.com, though the **gh** command
described here is by far the fastest and easiest way to work with the github ecosystem.

A relatively new addition to github is called "github CLI", which is implemented via the command
**gh** through which many github.com actions can now be run as terminal commands.
Besides [installing](https://cli.github.com/manual/installation)
it once, you also need to authenticate once via your github account:

      gh auth login

after which you could create your own fork, clone locally, and set the *upstream* all
in just one step:

      gh repo fork https://github.com/teuben/nemo

If all is well, the following commands should show the correct *origin* and *upstream*:

      git config --local remote.upstream.url
          >>> git@github.com:teuben/nemo.git
      git config --local remote.origin.url
          >>> git@github.com:YOURNAME/nemo.git

None of these should be blank!

## 2. Cloning your personal fork

If you have cloned your fork

      git clone https://github.com/YOURNAME/nemo

you only need to set the upstream:

      git remote add    upstream https://github.com/teuben/nemo

and you are ready for creating a PR.

## 1. Cloning the official upstream

If you happened to have cloned the official *upstream*

      git clone https://github.com/teuben/nemo

then things are a bit more complicated, because you should have cloned your
fork. However, the following commands will fix this
(assuming you also went to github.com and made that fork):

      git remote rename origin upstream
      git remote add    origin https://github.com/YOURNAME/nemo

Again, the **gh** command now gives a single line shortcut to all this:

      gh repo fork https://github.com/teuben/nemo

## Keeping your clone in sync

You should regularly make sure your local master branch
is in sync with the upstream master branch. This allows you
to work in local branches, and be up to date by branching off the
tip of this upstream master branch.

      git checkout master
      git fetch upstream
      git merge upstream/master
      git status
      git push

## Working in a branch

Assuming your own master is in sync with the upstream master,
here is a typical example, using a branchname **b1**

      git branch b1
      git checkout b1
      
          ## edit/add your files; test them; commit them, e.g.
	  git commit -m "my comment"  existing_file
	  git add new_file
	  git commit -m "my command"  new_file
          
      git push -u origin b1

Now you can issue a pull request on this branch **b1**.  There is a way to do this
via the **gh pr** command sequence. More about that in a future revision of this
document. 

You can even delete a branch, once it has been accepted as a pull request and merged
back in the upstream, it is really not needed anymore:

      git checkout master
      git branch -D b1
      git push origin --delete b1

## Memorable git options

1.  Show all files modified in a branch AAA 

      git diff master...AAA --name-status --diff-filter=M

2.  When was a branch created

      git show --summary `git merge-base AAA master`
      gitk --all --select-commit=`git merge-base foo master`

3. To see which files belonged in a commit, find the sha (e.g. via "git log" or "git log file" if 
   you know the file that belonged to it), then

     git diff-tree --no-commit-id --name-only -r SHA
 

# Tests

From the top level directory in NEMO there are a few basics regression tests and benchmarks:

      make check
      make bench
      make bench3

The **check** target depends on the many **Testfile** files sprinkled throughout NEMO.  Although
you can now find a few **Benchfile** files as well, they have not been put under a top level
target, so the current bench is hardcoded via a script **src/scripts/nemo.bench**. A future
revision of this document may go in a little more details.

# Debugging

Most applications are built using *make* variables defined in the $NEMOLIB/makedefs file.
Although not recommended, 
some hacking is allowed by directly editing this file, with the caveat you are then
bypassing the **configure** step of the install, which also does modify some other files.
Notably the **falcON** package has some extra dependancies that start with configure.

NEMO uses a set of hierarchical Makefiles, and many things (see google) have been written why they
are bad. But we don't recurse!

## Debuggging NEMO applications

Consult the local Makefile of the application. You should be able to override the compiler (CC=, CXX=, FC=)
and related options. If the application was integrated into NEMO, there should be a line

      include $(NEMOLIB)/makedefs

in that Makefile, and check the Makefile which are used.


An example:

      cd $NEMO/usr/trenti/CGS
      make clean install FC=flang FFLAGS=-O3
      make bench2

would create an alternate (and as it happens twice as fast as with gfortran) of the CGS N-body integrator.


# References

Some references on git workflows:

* https://docs.github.com/en/github/getting-started-with-github/fork-a-repo
* http://docs.astropy.org/en/stable/development/workflow/development_workflow.html
* https://www.atlassian.com/git/tutorials/comparing-workflows
* https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow
* http://physics.mnstate.edu/craig/git-novice-pyastro/
* https://www.sitepoint.com/quick-tip-sync-your-fork-with
