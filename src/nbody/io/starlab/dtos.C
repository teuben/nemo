/*
 *  DTOS: convert dyn to NEMO snapshot 
 *		22-mar-94
 *		12-apr-97	SINGLEPREC
 *		18-jul-01	read more times?
 *              14-oct-03       ieck, keeping up with many starlab changes...
 *              29-dec-03       also added phi into the output stream, key as well
 *               2-jan-04       attempt to get some UBV info for nvodemo2004
 *               5-jan-04  1.6b fixed taking log's the right way !!!
 *               8-jan-04  1.7  options to select what to get and where it goes
 *              30-jul-05  2.1  for new starlab dir.str.
 *
 * TODO:  somehow the dynamic_cast() in ubvri.C is not recognized in this mode,
 *        but a simple HelloWorld program has.
 */

// #include <iostream>
#include <nemo_stdinc.h>                 /* NEMO */
#include <getparam.h>
//#include <history.h>
//#include <extstring.h>

#include "stod_subs.h"	             /* NEMO to STARLAB interface */
#include <snapshot/snapshot.h>       /* merely for setting the attribute bits */

#include "pdyn.h"	             /* STARLAB */

#if 1
/* my own test version, that doesn't need -lsstar ; which gave linking problems */
#include "ubvri.C"
#endif

nemo_string defv[] = {
    "out=???\n          Output snapshot file (input dyn from stdin)",
    "headline=\n        Random verbiage for user",
    "options=\n         Optional starlab things to write out ",
    "VERSION=2.1\n      30-jul-05 PJT",
    NULL,
};

nemo_string usage = "convert STARLAB dyn to NEMO snapshot";

void nemo_main(void)
{
  pdyn *proot;
  int k, i=0, n=0, nbody=0;
  real *mass, *pos, *vel, *aux, *phi, *mptr, *pptr, *vptr, *aptr, *hptr, tsnap, tsnap0;
  real *acc, *pacc, u_out, b_out, v_out;
  double umag,bmag,vmag,rmag,imag;
  double m,logl,logt;
  int *key, *kptr;
  vec temp;
  int bits, first = 1;
  char *fname = getparam("out");
  char *hline = getparam("headline");
  char *options = getparam("options");
    
  check_real(sizeof(real));   /* make sure real==double */
#if 1
  bits = TimeBit | MassBit | PhaseSpaceBit;
#endif
  if (scanopt(options,"aux")) bits |= AuxBit;
  if (scanopt(options,"phi")) bits |= PotentialBit;
  if (scanopt(options,"acc")) bits |= AccelerationBit;
  if (scanopt(options,"key")) bits |= KeyBit;

  while ((proot = getpdyn(cin)) != NULL) {
    
    if (nbody == 0) {
      for_all_leaves(pdyn, proot, bi) {
	nbody++;
      }
      dprintf(0,"Found %d bodies\n",nbody);
      
      mptr = mass = (real *) allocate(nbody*sizeof(real));
      pptr = pos  = (real *) allocate(3*nbody*sizeof(real));
      vptr = vel  = (real *) allocate(3*nbody*sizeof(real));
      aptr = aux  = (real *) allocate(nbody*sizeof(real));
      hptr = phi  = (real *) allocate(nbody*sizeof(real));
      pacc = acc  = (real *) allocate(3*nbody*sizeof(real));
      kptr = key  = (int *)  allocate(nbody*sizeof(int));
    } else {
      mptr = mass;
      pptr = pos; 
      vptr = vel;
      pacc = acc;
      aptr = aux;
      hptr = phi;
      kptr = key;
    }

    for_all_leaves(pdyn, proot, bi) {
      
      m = *mptr++ = bi->get_mass();                  // mass
	
      temp = bi->get_pos();                          // position
      for (k = 0; k < 3; k++) *pptr++ = temp[k];
      
      temp = bi->get_vel();                          // velocity
      for (k = 0; k < 3; k++) *vptr++ = temp[k];
      logl = *aptr++ = bi->get_luminosity();         // luminosity
      logt = *hptr++ = bi->get_temperature();        // temperature
      *kptr++ = bi->get_stellar_type();              // stellar type
      logl = log10(logl);
      logt = log10(logt);
      nemo_ltm_to_ubvri(logl,logt,m, umag,bmag,vmag,rmag,imag);
      dprintf(1,"logT,logR,UBVRI: %g %g %g %g %g %g %g\n",
	      logt,logl,
	      umag,bmag,vmag,rmag,imag);
      u_out = umag;
      b_out = bmag;
      v_out = vmag;
      *pacc++ = u_out;            /*   U, V, B-V */
      *pacc++ = v_out;
      *pacc++ = b_out-v_out;
    }
    tsnap = proot->get_system_time();
    cerr << "output now "<< tsnap << endl;
    put_snap_c(fname,hline,nbody,tsnap,mass,pos,vel,acc,aux,phi,key);
    rmtree(proot);
    if (first) {
      first = 0;
    } else {
      if (tsnap == tsnap0) 
	error("Snapshots with same time not allowed: t=%g; probably bad input type from kira",tsnap);
      tsnap0 = tsnap;
    }
  }
}


