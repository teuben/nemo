/*******************************************************************************
*                                                                              *
* tidaldisk.c                                                                  *
*                                                                              *
* copyright by Walter Dehnen 2001                                              *
*                                                                              *
*******************************************************************************/
/*CTEX
*                                                                              
* Tidal field exerted by a (plane-parallel) stellar disk on a cluster passing  
* through with constant vertical velocity.                                     
* Useful for simulations of disk-shocking of, say, globular clusters           
*                                                                              
* The following three density models are available                             
*                             
* 1. thin disk:               
*                             
* $$                                                                             
* \rho(z) = \Sigma * \delta(z)
* $$                                                                             
*                                                                              
* 2. exponential disk:                                                         
*
* $$                                                                             
* \rho(z) = {\Sigma \over {2h}} \exp{ { -|z|} \over h}
* $$
*   
* 3. sech$^2$ disk:                                                       
*
* $$                                                                    
*  \rho(z) = {\Sigma \over {4h}} sech^2{ { z \over {2h}}}
* $$                                                                    
*                                                                       
* Parameters (to be given by potpars=...) are:                         
* \begin{verbatim}
* par[0] = not used (reserved for pattern speed in NEMO)               
* par[1] = h    scale-height  
* .                            par[1] = 0 -> thin disk                  
* .                            par[1] > 0 -> vertically exponential disk
* .                            par[1] < 0 -> sech$^2$ disk with h=|par[1]|
* par[2] = Sig  disk surface density                        
* par[3] = Vz   constant vertical velocity of cluster center
* par[4] = Z0   cluster center z-position at t=0     
* par[5] = add  boolean: add tidal potential or not? 
* \end{verbatim}
*
* We always assume G=1.
*
* If you want to include the acceleration of the disk on the cluster as a      
* whole, rather than assume a constant velocity, use {\bf vertdisk}
*
* Some words on the mechanics                         
*
* Assume that the plane-parallel disk potential and force are given by         
* $$                                                                             
*     \Phi(Z)  
* $$                                                                             
* and
* $$                                                                             
*     F(Z) = -\Phi'(Z).                                         
* $$                                                                             
* Then, the tidal force exerted on a star at position z w.r.t. to cluster      
* center, which in turn is at absolute height Zc = Z0 + t Vz, is simply        
* $$                                                                             
*     F_t(z) = F(Zc+z) - F(Zc).                                             
* $$                                                                             
* Integrating this from z=0 to z gives the associated tidal potential as       
* $$                                                                             
*     \Phi_t(z) = \Phi(Zc+z) - \Phi(Zc) + z * F(Zc).                           
* $$                                                                             
* Whenever the tidal force \& potential are desired at a new time t, we         
* pre-compute $Zc$ and the plane-parallel potential and force at $Z=Zc$.
* Note that when both $Zc$ and $Zc+z$ are outside of the mass of the disk (and 
* $Z=0$
* is not between them), both tidal force and potential vanish identically.     
* The standard treatment of tidal forces corresponds to approximating (2) by   
* $F(Zc) + z * F'(Zc)$. This method, however, breaks down for disks that are     
* thin compared to the cluster, while our method is always valid, even for a
* razor thin disk.
*/

#include "vertdisk.h"

local double Vz  = 1.0;                     /* vertical v of cluster center   */
local double Z0  = 0.0;                     /* height of cluster center at t=0*/
local bool   add = 1;                       /* add tidal potential            */
local double tnow,Znow,Anow,Pnow;           /* auxiliary variables            */

/*----------------------------------------------------------------------------*/
void ini_potential(int *npar, double *par, string name)
{
  DISKPOTENTIAL_VARS                                /* vertdisk.h: aux vars   */
  int n = *npar;

  if(n>0) omega = par[0];
  ini_vertdisk(n-1,par+1);
  if(n<2) warning("tidaldisk: h defaulting to %f",h);
  if(n<3) warning("tidaldisk: Sig defaulting to %f",Sig);
  if(n>3) Vz = par[3];
  else    warning("tidaldisk: Vz defaulting to %f",Vz);
  if(n>4) Z0 = par[4];
  else    warning("tidaldisk: Z0 defaulting to %f",Z0);
  if(n>5) add= par[5] > 0.;
  if(n>6) warning("tidaldisk: npar=%d only 5 parameters accepted",n);
  dprintf(1," inipotential: tidaldisk\n");
  dprintf(1," parameters: h, Sig, Vz, Z0, add= %f %f %f %f %d\n",h,Sig,Vz,Z0,add);
  par[0] = omega;
  tnow = 0.0;                                       /* current time           */
  Znow = Z0+Vz*tnow;                                /* current z_c            */
  DISKPOTENTIAL_FUNC(Znow,Anow,Pnow);               /* current acc_c & pot_c  */
}
/*----------------------------------------------------------------------------*/
void potential_double(int   *ndim,            	    /* I: # dims, ignored     */
		      double*pos,             	    /* I: (x,y,z)             */
		      double*acc,            	    /* O: (ax,ay,az)          */
		      double*pot,            	    /* O: potential           */
		      double*time)           	    /* I: time                */
{
  DISKPOTENTIAL_VARS                                /* vertdisk.h: aux vars   */
  register double Z;                                /* position w.r.t. disk   */
  acc[0] = 0.0;                                     /* no force in x and      */
  acc[1] = 0.0;                                     /* no force in y          */
  if(*time != tnow) {                               /* IF(new time)         > */
    tnow  =*time;                                   /*   current time         */
    Znow  = Z0+Vz*tnow;                             /*   current z_c          */
    DISKPOTENTIAL_FUNC(Znow,Anow,Pnow)              /*   current acc_c & pot_c*/
  }                                                 /* <                      */
  Z=Znow+pos[2];                                    /* position    w.r.t. disk*/
  DISKPOTENTIAL_FUNC(Z,acc[2],*pot);                /* force & pot w.r.t. disk*/
  if(add) *pot+= pos[2]*Anow - Pnow;                /* pot   w.r.t. center    */
  else    *pot = 0.0;                               /* pot   not given        */
  acc[2] -= Anow;                                   /* force w.r.t. center    */
}
/*----------------------------------------------------------------------------*/
void potential_float(int  *ndim,            	    /* I: # dims, ignored     */
		     float*pos,             	    /* I: (x,y,z)             */
		     float*acc,            	    /* O: (ax,ay,az)          */
		     float*pot,            	    /* O: potential           */
		     float*time)           	    /* I: time                */
{
  DISKPOTENTIAL_VARS                                /* vertdisk.h: aux vars   */
  register double Z;                                /* position w.r.t. disk   */
  acc[0] = 0.0;                                     /* no force in x and      */
  acc[1] = 0.0;                                     /* no force in y          */
  if(*time != tnow) {                               /* IF(new time)         > */
    tnow  =*time;                                   /*   current time         */
    Znow  = Z0+Vz*tnow;                             /*   current z_c          */
    DISKPOTENTIAL_FUNC(Znow,Anow,Pnow)              /*   current acc_c & pot_c*/
  }                                                 /* <                      */
  Z=Znow+pos[2];                                    /* position    w.r.t. disk*/
  DISKPOTENTIAL_FUNC(Z,acc[2],*pot);                /* force & pot w.r.t. disk*/
  if(add) *pot+= pos[2]*Anow - Pnow;                /* pot   w.r.t. center    */
  else    *pot = 0.0;                               /* pot   not given        */
  acc[2] -= Anow;                                   /* force w.r.t. center    */
}
/*----------------------------------------------------------------------------*/
  
  
