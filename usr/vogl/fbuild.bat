cd src\msfort
make F=%1 makefile.msc
if errorlevel > 0 goto end
cd ..
make F=%1 makefile.msc
if errorlevel > 0 goto end
make F=%1 fmakefil.msc
if errorlevel > 0 goto end
if "%2" == "nofonts" goto nof
cd ..\hershey\src
make makefile.msc
if errorlevel > 0 goto end
make fmakefil.msc
if errorlevel > 0 goto end
md c:\lib\hershey
command /c mkfnts c:\lib\hershey
cd ..\..\examples
make F=%1 fmakefil.msc
goto end
:nof
cd ..\examples
make F=%1 fmakefil.msc
if errorlevel > 0 goto end
:end
