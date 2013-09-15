/*
 *
c     This program calculates the rotation curve for an oblate spheroidal 
c     bulge. It assumes a projected surface density that follows a Sersic 
c     profile. The central surface magnitude M_0, characteristic radius R_0 and
c     Sersic index n_sers must be given. Also, the inclination and intrinsic 
c     axis ratio of the bulge must be given. Finally, the distance and the 
c     Galactic foreground extinction of the galaxy is needed.

c     Copyright (C) 2007, Edo Noordermeer
c     E-mail: edo.noordermeer@gmail.com
c
c     If you use this program for a scientific publication, I would appreciate 
c     a reference to the paper 'The rotation curves of flattened Sersic bulges
c     (Noordermeer 2008)'.
 */


#include <nemo.h>

string defv[] = {
  "out=bulgerot.dat\n         output log file",
  "galaxy=thisgalaxy\n        Identification",
  "mu0=15\n                   apparent central R-band surface magnitude M_0",
  "r0=1\n                     characteristic radius R_0 in arcseconds",
  "n=1\n                      Sersic index",
  "inc=0\n                    Inclination",
  "q=1\n                      Axis ratio b/a ",
  "dist=10\n                  Distance in Mpc",
  "ar=0\n                     galactic foreground extinction",
  "radii=0:2.5:0.01\n         Radii to calculate for (kpc)",
  "VERSION=1.1\n              12-sep-2013 PJT",
  NULL,
};


string usage="rotation curve of oblate spheroidal bulge";

string cvsid="$Id:";
 

nemo_main()
{
  warning("work in progress");
}
