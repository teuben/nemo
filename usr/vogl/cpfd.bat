echo off
rem
rem  Copys vogl directories onto floppy disk
rem
rem  Usage: cd to vogl directory, then cpfd A:
rem
echo MAKING DIRECTORIES
md %1\vogl
md %1\vogl\src
md %1\vogl\src\msfort
md %1\vogl\src\sunfort
md %1\vogl\drivers
md %1\vogl\drivers\ibmpc
md %1\vogl\examples
md %1\vogl\hershey
md %1\vogl\hershey\fonts
md %1\vogl\hershey\docs
md %1\vogl\hershey\src
md %1\vogl\hershey\data
md %1\vogl\docs
echo "DONE MAKING DIRECTORIES"
rem
rem  Do the copying
rem
echo BEGIN COPYING
copy *.* %1\vogl
copy src\*.* %1\vogl\src
copy src\msfort\*.* %1\vogl\src\msfort
copy src\sunfort\*.* %1\vogl\src\sunfort
rem
copy drivers\*.* %1\vogl\drivers
copy drivers\ibmpc\*.* %1\vogl\drivers\ibmpc
rem
copy examples\*.* %1\vogl\examples
rem
copy hershey\*.* %1\vogl\hershey
copy hershey\src\*.* %1\vogl\hershey\src
copy hershey\docs\*.* %1\vogl\hershey\docs
copy hershey\data\*.* %1\vogl\hershey\data
copy hershey\fonts\*.* %1\vogl\hershey\fonts
copy docs\*.* %1\vogl\docs
echo FINISHED COPYING
