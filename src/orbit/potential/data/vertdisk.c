/*******************************************************************************
*                                                                              *
* vertdisk.c                                                                   *
*                                                                              *
* copyright by Walter Dehnen 2001                                              *
*                                                                              *
* Vertical potential for a plane-parallel stellar disk.                        *
* Useful for simulations of disk-shocking of, say, globular clusters           *
*                                                                              *
* The following three density models are available                             *
*                                                                              *
* 1. thin disk:                                                                *
*
* rho(z) = Sig * delta(z)                                                      *
*                                                                              *
* 2. exponential disk:                                                         *
*                                                                              *
*           Sig        -|z|                                                    *
* rho(z) = ----- * exp ----                                                    *
*           2*h          h                                                     *
*                                                                              *
* 3. sech^2 disk:                                                              *
*                                                                              *
*           Sig        2  z                                                    *
* rho(z) = ----- * sech  ---                                                   *
*           4*h          2*h                                                   *
*                                                                              *
* Parameters (to be given by potpars=...) are:                                 *
*                                                                              *
* par[0] = not used (reserved for pattern speed in NEMO)                       *
* par[1] = h scale-height  par[1] = 0 -> thin disk                             *
*                          par[1] > 0 -> vertically exponential disk           *
*                          par[1] < 0 -> sech^2 disk with h=|par[1]|           *
* par[2] = Sig: disk surface density                                           *
*                                                                              *
* We always assume G=1.                                                        *
*                                                                              *
********************************************************************************
*                                                                              *
* Note that this potential is plane-parallel. If you are only interested in    *
* tidal field exerted from a disk on, say, a globular cluster, you may better  *
* use tidaldisk.c, which does not accelerate the globular cluster as a whole,  *
* but assumes a constant vertical velocity.                                    *
*                                                                              *
*******************************************************************************/
#include "vertdisk.h"

void inipotential(int *npar, double *par, string name)
{
  int n = *npar;
  if(n>0) omega = par[0];
  if(n>2) warning("vertdisk: npar=%d only 3 parameters accepted",n);
  ini_vertdisk(n-1,par+1);
  dprintf(1," inipotential: vertdisk\n");
  dprintf(1," parameters: h, Sig = %f %f\n",h,Sig);
  par[0] = omega;
}
/*----------------------------------------------------------------------------*/
void potential_double(int   *ndim,            	    /* I: # dims, ignored     */
		      double*pos,             	    /* I: (x,y,z)             */
		      double*acc,            	    /* O: (ax,ay,az)          */
		      double*pot,            	    /* O: potential           */
		      double*time)           	    /* I: time                */
{
  DISKPOTENTIAL_VARS                                /* vertdisk.h: aux vars   */
  acc[0] = 0.0;                                     /* no force in x and      */
  acc[1] = 0.0;                                     /* no force in y          */
  DISKPOTENTIAL_FUNC(pos[2],acc[2],*pot)
}
/*----------------------------------------------------------------------------*/
void potential_float(int  *ndim,            	    /* I: # dims, ignored     */
		     float*pos,             	    /* I: (x,y,z)             */
		     float*acc,            	    /* O: (ax,ay,az)          */
		     float*pot,            	    /* O: potential           */
		     float*time)           	    /* I: time                */
{
  DISKPOTENTIAL_VARS                                /* vertdisk.h: aux vars   */
  acc[0] = 0.0;                                     /* no force in x and      */
  acc[1] = 0.0;                                     /* no force in y          */
  DISKPOTENTIAL_FUNC(pos[2],acc[2],*pot)
}
/******************************************************************************/
