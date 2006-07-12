// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// Rotcurve.cc                                                                 |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2004                                          |
//           Paul McMillan, 2005-2006                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
//           paul.mcmillan@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
//     find the rotation cuve of a disk in a given potential                   |
//                                                                             |
//                                                                             |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v1.0    16/01/2006   PJM wrote this                                         |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "2.3"
#define falcON_VERSION_D "24-jun-2005 Paul McMillan                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile Rotcurve
#endif
#include <public/io.h>                             // WD I/O utilities          
#include <main.h>                                  // main & NEMO stuff        
#include <externacc.h>                             // external fields
#include <utils/tupel.h>                           // tupels
//-----------------------------------------------------------------------------+
string defv[] = {
  "out=???\n          output file                                        ",
  "rmax=6\n           maximum r in table                                 ",
  "nvalue=100\n       number of values in the table                      ",
  "potname=\n         name of external acceleration field (required)     ",
  "potpars=\n         parameters of external acceleration field          ",
  "potfile=\n         file required by external acceleration field       ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string usage =
  "Rotcurve:   gives the disk rotation curve given the total potential field\n";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  typedef tupel<3,double> vectd;
  const nemo_acc*aex = hasvalue("potname")?           // IF(potname given) THEN 
    new nemo_acc(getparam  ("potname"),               //   initialize external  
		 getparam_z("potpars"),               //   accelerations        
		 getparam_z("potfile")) : 0;          // ELSE: no potential     
  int n=getiparam("nvalue");
  double rm=getdparam("rmax"),
    dr=rm/double(n);
  double *pot   = new double[n];
  vectd *pos_e  = new vectd[n];
  vectd *acc_e  = new vectd[n];

  for(int i=0;i!=n;++i) {
   pos_e[i]    = zero;
   pos_e[i][0] = dr*double(i+1);
 }
  aex->set(0.,n,0,pos_e,0,0,pot,acc_e,0);
  double vc[n];
  for (int i=0;i!=n;++i) vc[i]=sqrt(-pos_e[i][0] * acc_e[i][0]);

  delete[] pot;
  delete[] acc_e;
  
  output TAB;
  TAB.open(getparam("out"),1);
  if(!TAB) falcON_THROW("cannot open tabfile\n");
  TAB << "#\n"
      << "# \""<< (*(ask_history())) <<"\"\n"
      << "#\n#   R      v_c\n";
  for(int i=0; i!=n; ++i)
    TAB << std::setw(6) << pos_e[i][0] << ' ' << vc[i] << "\n";
  TAB << std::endl;

  delete[] pos_e;
}
