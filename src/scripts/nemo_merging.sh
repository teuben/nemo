#! /bin/bash
#
#      create two NEMO trees, nemo1 and nemo2, and prepare some merging scenarios
#      preferably do this in a clean directory
#      set the mode= to control what scenario fits you best.
#
#      Usage:   the url1/2=, b1/2= and mode= can all be changed via cmdline
#

#    Set the two repo's
url1=https://github.com/teuben/nemo
url2=https://github.com/psanth/nemo

#    Set the branche names in these repos to compare/merge
b1=c99
b2=c++_nemo

#    Select a mode how to prepare the merging
#
#       0. do nothing
#       1. merge in your given branch
#       2. make a local branch of the remote with the same name
#
mode=0


# ============================================================================================

#    Poor man's command line parser
for arg in $*; do\
  export $arg
done


#   ---                  grab both repos
git clone $url1 nemo1
git clone $url2 nemo2

#   ---                  add the other as a remote, and switch to their respective branches
cd nemo1
git remote add nemo2 $url2
git fetch nemo2
if test $mode = 1; then
    git checkout $b1
    git merge nemo2/$b2
elif test $mode = 2; then
    git checkout -b $b2    
    git merge nemo2/$b2
fi
cd ..

cd nemo2
git remote add nemo1 $url1
git fetch nemo1
if test $mode = 1; then
    git checkout $b2
    git merge nemo1/$b1
elif test $mode = 2; then
    git checkout -b $b1
    git merge nemo1/$b1
fi
    
cd ..
