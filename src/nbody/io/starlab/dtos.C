/*
 *  DTOS: convert dyn to NEMO snapshot 
 *		22-mar-94
 *		12-apr-97	SINGLEPREC
 */

#include <stdinc.h>                 /* NEMO */
#include <getparam.h>
#include <history.h>
#include <extstring.h>

#include "stod_subs.h"	            /* NEMO to STARLAB interface */

#include "dyn.h"	            /* STARLAB */

string defv[] = {
    "out=???\n          Output snapshot file (input dyn from stdin)",
    "VERSION=1.1\n      31-jul-96 PJT",
    NULL,
};

string usage = "convert STARLAB dyn to NEMO snapshot";

void nemo_main(void)
{
    dyn *root, *ni;
    int k, i=0, n=0, nbody=0;
    real *mass, *pos, *vel, *mptr, *pptr, *vptr;
    vector temp;

    check_real(sizeof(real));

    root = get_dyn(cin);

    for (ni = root->get_oldest_daughter();
         ni != NULL;
         ni = ni->get_younger_sister())       nbody++;      // count total N

    mptr = mass = (real *) allocate(nbody*sizeof(real));
    pptr = pos  = (real *) allocate(3*nbody*sizeof(real));
    vptr = vel  = (real *) allocate(3*nbody*sizeof(real));
    
    for (ni = root->get_oldest_daughter();
         ni != NULL;
         ni = ni->get_younger_sister()) {

	 cerr << "Working on " << ++i << endl;

         *mptr++ = ni->get_mass();                      // mass
         
         temp = ni->get_pos();                          // position
         for (k = 0; k < 3; k++) *pptr++ = temp[k];
         
         temp = ni->get_vel();                          // velocity
         for (k = 0; k < 3; k++) *vptr++ = temp[k];
     }

     cerr << "output now"<< endl;

     put_snap_c(getparam("out"),nbody,mass,pos,vel);

}
