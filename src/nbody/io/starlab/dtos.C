/*
 *  DTOS: convert dyn to NEMO snapshot 
 *		22-mar-94
 *		12-apr-97	SINGLEPREC
 *		18-jul-01	read more times?
 *              14-oct-03       ieck, keeping up with many starlab changes...
 */

#include <stdinc.h>                 /* NEMO */
#include <getparam.h>
#include <history.h>
#include <extstring.h>

#include "stod_subs.h"	            /* NEMO to STARLAB interface */

#include "pdyn.h"	            /* STARLAB */

typedef char *nemo_string;

nemo_string defv[] = {
    "out=???\n          Output snapshot file (input dyn from stdin)",
    "VERSION=1.3\n      14-oct-03 PJT",
    NULL,
};

nemo_string usage = "convert STARLAB dyn to NEMO snapshot";

void nemo_main(void)
{
  pdyn *proot;
  int k, i=0, n=0, nbody=0;
  real *mass, *pos, *vel, *aux, *mptr, *pptr, *vptr, *aptr, tsnap;
  vec temp;
    
  check_real(sizeof(real));   /* make sure real==double */

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
    } else {
      mptr = mass;
      pptr = pos; 
      vptr = vel; 
      aptr = aux;
    }

    for_all_leaves(pdyn, proot, bi) {
      
      *mptr++ = bi->get_mass();                      // mass
	
      temp = bi->get_pos();                          // position
      for (k = 0; k < 3; k++) *pptr++ = temp[k];
      
      temp = bi->get_vel();                          // velocity
      for (k = 0; k < 3; k++) *vptr++ = temp[k];
#if 0
      *aptr++ = bi->get_luminosity();                // luminosity
#else
      *aptr++ = bi->get_temperature();
#endif
    }
    tsnap = proot->get_system_time();
    cerr << "output now "<< tsnap << endl;
    put_snap_c(getparam("out"),nbody,tsnap,mass,pos,vel,aux);
    rmtree(proot);
  }
}
