cd src
make -fmakefile.tc > c:src.err
if errorlevel > 0 goto end
cd ..\hershey\src
make -fmakefile.tc > c:her.err
if errorlevel > 0 goto end
md c:\lib\hershey
command /c mkfnts c:\lib\hershey
cd ..\..\examples
make -fmakefile.tc > c:exa.err
:end
