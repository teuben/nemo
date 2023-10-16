# Dependencies

Building NEMO requires following minimum dependencies:

* Unix tools, such as: make, csh, bash, git
* C/C++/Fortran compiler
* pgplot (optional, but highly recommended)

# git

There are two ways you could have obtained the NEMO source code:

1. You cloned from the *upstream* : https://github.com/teuben/nemo . This is fine
if you just want to compile and run, but not ideal if you make modifications and
want to share them back to the *upstream* via a *pull request* (PR). You
can however "repair" your local repo, discussed below, and still submit a PR.

2. You [forked](https://guides.github.com/activities/forking/) NEMO
from the *upstream*, and cloned it locally from the repo in your own
github account. This is the ideal method, but you will still need to
set the *upstream* manually if you used github.com. See also the
**gh** command below for an even easier way.

3. Sadly on github.com you will also find a **zip** copy of the repo
that does actually work fine, except it's a frozen snapshot and cannot
be efficiently used to collaborate. However, if you cannot install git,
this is probably the only way to bootstrap yourself. For example
https://github.com/teuben/nemo/archive/refs/heads/master.zip, which
will create a directory *nemo-master*. Other branches are available
through similar zip file construct.

Familiarize yourself with the concept of a pull request on github. There
are some links at the bottom of this document.


## gh:   github CLI

You can safely skip this section if you prefer to work via github.com, though the **gh** command
described here is by far the fastest and easiest way to work with the github ecosystem. You just
have to intall yet another tool for this.

If you can use conda, installation can be done as follows:

      conda install gh --channel conda-forge  

but see also [manual installing instructions](https://cli.github.com/manual/installation),
after this you need to authenticate once via your github account:

      gh auth login

after which you could create your own fork, clone locally, and set the *upstream* all
in just one step:

      gh repo fork https://github.com/teuben/nemo

If all is well, the following commands should show the correct *origin* and *upstream*:

      git remote -v

      origin git@github.com:YOURNAME/nemo.git
      upstream git@github.com:teuben/nemo.git

for both the (fetch) and (push). None of these should be blank! You are now ready for working
in your own branches and issue a pull request (PR). 

## 2. Cloning your personal fork

This section can be skipped if you use the **gh repo fork** method described before this.

If you have cloned your own fork on github.com, you should now clone it

      git clone https://github.com/YOURNAME/nemo

you only need to set the upstream:

      git remote add upstream https://github.com/teuben/nemo

and you are ready for creating a PR (from a branch of course).

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
tip of this upstream master branch. Here's a recipe for that
in your local clone:

      git checkout master
      git fetch upstream
      git merge upstream/master
      git status
      git push

## Working in a branch

Assuming your own master is in sync with the upstream master,
here is a typical example to work in a branch, using a branchname **b1**

      git branch b1
      git checkout b1
      
      ## edit/add your files; test them; commit them, e.g.
      git commit -m "my comment"  existing_file
      git add new_file
      git commit -m "my command"  new_file
          
      git push -u origin b1

Now you can issue a pull request on this branch **b1**.

If you have become a fan on the **gh pr** method, here's the recipe for this:

      git checkout -b b1
        <<edit, test, commit>>
      gh pr create 

once the branch has been merged, you don't need it locally anymore, so delete it

      git branch -d b1
      git push original --delete b1


## Common work in a branch

Lets say there is a new development idea, lets call it "table2", and others will
share the development in this "table2" branch, but until the development is done, this is not
merged to the master yet. 

1. NEMO (*upstream*) itself will start making this new development branch:

      git checkout -b table2
      ...
      git push --set-upstream origin table2
	  
2. others will do then branch off this new branch (user *astroumd* is used as example here):

      gh repo fork https://github.com/teuben/nemo nemo
      cd nemo
      git remote add upstream  https://github.com/teuben/nemo
      git fetch upstream
      git checkout -b myTable3
      git merge upstream/table2
      ...
      git push --set-upstream origin myTable3
	  
   on github.com/astroumd/nemo you can then do a pull request from astroumd::myTable3 to teuben:table2. An 
   alternative (next item) is that the upstream person would pull in myTable3 and tests locally

3. NEMO, optionally, creates an alias site "astroumd" and merge from their myTable3 branch to check out the code

      git checkout table2
      git remote add  astroumd  https://github.com/astroumd/nemo
      git pull astroumd
      git merge [--no-ff] astroumd/myTable3
      .... (resolve conflicts)
      git push

4. after this, all collaborators will need to merge these back:

      git checkout myTable3
      git fetch upstream
      git merge upstream/table2

## Memorable git options

1.  Show all files modified in a branch AAA 

      git diff master...AAA --name-status --diff-filter=M

2.  When was a branch created

      git show --summary `git merge-base AAA master`
      gitk --all --select-commit=`git merge-base foo master`

3. To see which files belonged in a commit, find the sha (e.g. via "git log" or "git log file" if 
   you know the file that belonged to it), then

      git diff-tree --no-commit-id --name-only -r SHA
	 
4. Difference between two SHA's

      git diff <commit-id> <commit-id>

## PR with github CLI

Recapping working with "github CLI" (gh)

Step 1 from the submitter of the PR:

      gh repo fork https://github.com/teuben/nemo
      cd nemo
      git checkout -b teuben1
      $EDITOR index.md
      gh pr create 

Once the branch has been merged by the upstream, there is no need to keep it.
It can be removed as follows:

      git branch -d teuben1
      git push original --delete teuben1

Step 2 by the receiver of the PR:



	

# Tests

From the top level directory in NEMO there are a few basics regression tests and benchmarks:

      make check
      make bench5
      make bench

The **check** target depends on the many **Testfile** files sprinkled throughout NEMO.

The **bench5** target currently runs 6 NEMO programs all designed to run 5 seconds on
some particular processor that your NEMO compilation can now be compared to. The default
*speed* for that processor is defined to be 1000, anything larger is a faster processor.
Currently in 2023 we're approaching a single core score of 2000!


The **bench** target is still under development, and has a more
diverse set of programs, all controlled by sprinkled **Benchfile**
files. See the script **src/scripts/nemo.bench**.

# Debugging

Most applications are built using *make* variables defined in the **$NEMOLIB/makedefs** file.
Although not recommended, 
hacking is always possible by directly editing this file, with the caveat you are then
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
* https://how-to.dev/how-git-stores-data
