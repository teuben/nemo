cd src
make -fmakefile.bcc > c:\bernie\vogl\src.err
if errorlevel > 0 goto end
cd ..\hershey\src
make -fmakefile.bcc > c:\bernie\vogl\her.err
if errorlevel > 0 goto end
cd ..\..\examples\mswin
make -fmakefile.bcc > c:\bernie\vogl\exa.err
cd ..\hershey\src
md c:\lib\hershey
mkfnts c:\lib\hershey
:end
