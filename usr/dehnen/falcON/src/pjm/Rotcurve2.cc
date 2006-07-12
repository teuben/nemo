// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// Rotcurve2.cc                                                                |
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
//     find the rotation curve of a disk in a given potential                  |
//                                                                             |
//                                                                             |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v1.0    08/03/2006   PJM wrote this                                         |
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
#include <utils/tupel.h>                          // tupels
//-----------------------------------------------------------------------------+
string defv[] = {
  "out=???\n          output file                                        ",
  "rmax=6\n           maximum r in table                                 ",
  "nvalue=100\n       number of values in the table                      ",
  "dpotname=\n        name of disk external acceleration field (required)",
  "dpotpars=\n        parameters of disk external acceleration field     ",
  "dpotfile=\n        file required by disk external acceleration field  ",
  "bpotname=\n        name of bulge acceleration field (required)        ",
  "bpotpars=\n        parameters of bulge acceleration field             ",
  "bpotfile=\n        file required by bulge acceleration field          ",
  "hpotname=\n        name of halo acceleration field (required)         ",
  "hpotpars=\n        parameters of halo acceleration field              ",
  "hpotfile=\n        file required by halo acceleration field           ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string usage =
  "Rotcurve2:   gives the disk rotation curve given the total potential field\n";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  typedef tupel<3,double> vectd;
  const nemo_acc*aex1 = hasvalue("dpotname")?           // IF(potname given) THEN 
    new nemo_acc(getparam  ("dpotname"),               //   initialize external  
		 getparam_z("dpotpars"),               //   accelerations        
		 getparam_z("dpotfile")) : 0;          // ELSE: no potential     

  const nemo_acc*aex2 = hasvalue("bpotname")?           // IF(potname given) THEN 
    new nemo_acc(getparam  ("bpotname"),               //   initialize external  
		 getparam_z("bpotpars"),               //   accelerations        
		 getparam_z("bpotfile")) : 0;          // ELSE: no potential     

  const nemo_acc*aex3 = hasvalue("hpotname")?           // IF(potname given) THEN 
    new nemo_acc(getparam  ("hpotname"),               //   initialize external  
		 getparam_z("hpotpars"),               //   accelerations        
		 getparam_z("hpotfile")) : 0;          // ELSE: no potential     

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
  aex1->set(0.,n,0,pos_e,0,0,pot,acc_e,0);
  double vcd[n];
  double vcb[n];
  double vch[n];
  double vc[n];

  for (int i=0;i!=n;++i) vcd[i]=(acc_e[i][0] < 0.)? sqrt(-pos_e[i][0] * acc_e[i][0]): 0.;

  aex2->set(0.,n,0,pos_e,0,0,pot,acc_e,0);
  for (int i=0;i!=n;++i) vcb[i]=(acc_e[i][0] < 0.)? sqrt(-pos_e[i][0] * acc_e[i][0]): 0.;

  aex3->set(0.,n,0,pos_e,0,0,pot,acc_e,0);
  for (int i=0;i!=n;++i) vch[i]=(acc_e[i][0] < 0.)? sqrt(-pos_e[i][0] * acc_e[i][0]): 0.;

  for (int i=0;i!=n;++i) 
    vc[i] =sqrt(vcd[i]*vcd[i] + vcb[i]*vcb[i] + vch[i]*vch[i]);

  delete[] pot;
  delete[] acc_e;
  
  output TAB;
  TAB.open(getparam("out"),1);
  if(!TAB) falcON_THROW("cannot open tabfile\n");
  TAB << "#\n"
      << "# \""<< (*(ask_history())) <<"\"\n"
      << "#\n# R     v_c    v_c,d    v_c,b    v_c,h\n";
  for(int i=0; i!=n; ++i)
    TAB << pos_e[i][0] << ' ' << vc[i] << ' ' 
	<< vcd[i] << ' ' << vcb[i] << ' ' << vch[i] << "\n";
  TAB << std::endl;

  delete[] pos_e;
}
