cd src
make -fmakefile.dj %1
if errorlevel > 0 goto end
cd ..\hershey\src
make -fmakefile.dj %1
if errorlevel > 0 goto end
if "%1" == "" goto ex
:ex
cd ..\..\examples
make -fmakefile.dj %1
:end
