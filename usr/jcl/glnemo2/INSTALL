
                     G l n e m o 2      I n s t a l l a t i o n
                     = = = = = = = = = = = = = = = = = = = = = = 
                     see http://projets.lam.fr/projects/glnemo2/wiki/Wiki#Installation
  
  
0) Requirements:
------------------  
 Glnemo2 compiles and runs fine on Linux, Windows and MacOSX platform.
  In order to compile it, you need QT development  library (qt > 4.6), qt5 is preferred.
Go to http://qt-project.org/downloads to download QT, click on the link for
your appropriate platform.
 You also need a decent video card with a fast opengl driver. Glnemo2
has been successfully tested on Nvidia and  ATI GPU cards. If the
installed OpenGL driver is not fully compatible with the requirement, it
can also run in software mode, but very slowly.


1) Compilation :
------------------
  Check that you have 'qmake' utility in your path ('qmake' is normally included 
with qt package, but on for example fedora9 it is called qmake-qt4 since both qt3 and
qt4 are still installed there). Change to 'glnemo2' directory and enter the command :

    a) Linux and Windows only
       - - - - - - - - - - - -
       qmake -recursive
          then
       make

       
    b) Mac OS X
       - - - - -
       qmake -spec macx-g++ -recursive
          then
       make

    c) cmake (prefered method for Linux and Mac)
       - - -
    You must have cmake utility installed (> 2.8)
    cd build
    cmake ..
    make
    make install

 if the compilation fails, please send me a full report by e-mail : 
  jean-charles.lambert@lam.fr 
   
it should take a while to compile, at the end you should have a 'glnemo2' binary located in
"PATH/glnemo2/bin/ARCH/debug" directory. PATH is where is located glnemo2 directory. ARCH depends on
you architecture, and could be "unix", "win32" or "macx".

2) Installation :
-------------------
   just enter the command :

   make install

3) GyrfalcON runtime manipulator
--------------------------------

see $NEMO/usr/jcl/glnemo2/gyrfalcon/README

4) Binaries
-----------

You can download and install directly glnemo2 binary for different linux disrtributions, MacOSX and windows 32 et 64 bits platform.

See http://projets.lam.fr/projects/glnemo2/wiki/Download

