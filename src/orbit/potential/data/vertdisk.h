/*******************************************************************************
*                                                                              *
* vertdisk.h                                                                   *
*                                                                              *
* copyright by Walter Dehnen 2001                                              *
* file to be included by vertdisk.c and tidaldisk.c                            *
*                                                                              *
********************************************************************************
*                                                                              *
* The following three density models are available                             *
*                                                                              *
* 1. thin disk:                                                                *
*                                                                              *
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
* par[2] = Sig: surface density                                                *
*                                                                              *
* We always assume G=1.                                                        *
*                                                                              *
*******************************************************************************/
#ifndef  included_vertdisk_h
#define  included_vertdisk_h

#include <stdinc.h>
#include <vectmath.h>                   /* define DIMensionality              */

local double omega = 0.0;               /* pattern speed, not used            */
local double h     = 1.0;               /* scale height                       */
local double Sig   = 1.0;               /* surface density                    */
local double ih,fpot,facc;              /* auxiliary variables                */
local enum { thin=0, expo=1, sech=2}    /* enum for disk type                 */
             type;                      /* type of disk                       */

#if defined(TWODIM)
#error"vertical_disk potential needs three dimensions"
#endif

local void ini_vertdisk(int n, double *par)
{
  if(n>0) {
    if(par[0] == 0.0) {
      type = thin;
      h    = 0.0;
    } else if(par[0] > 0.0) {
      type = expo;
      h    = par[0];
    } else {
      type = sech;
      h    =-par[0];
    }
  }
  if(n>1) Sig = par[1];
  switch(type) {
  case thin:
    ih   = 1.0;
    fpot = 6.2831853071795864770 * Sig;
    facc = fpot;
    break;
  case expo:
    ih   = 1.0/h;
    facc = 6.2831853071795864770 * Sig;
    fpot = facc * h;
    break;
  case sech:
    ih   = 1.0/h;
    facc = 6.2831853071795864770 * Sig;
    fpot = 2 * h * facc;
    break;
  default:
    warning("unknown disk type in vertical_disk potential");
  }
}
/*
 * We now define the macro DISKPOTENTIAL_FUNC that computes the vertical      *
 * force and potential given the vertical height z.                           *
 * We also need an extra macro DISKPOTENTIAL_VARS defining the auxiliary      *
 * variables, since in C and in contrast to C++ we cannot define variables    *
 * where needed, but must define them at the beginning of any routine.        *
 * So DISKPOTENTIAL_VARS must be included in the header of any routine        *
 * containing DISKPOTENTIAL_FUNC, even if the latter is not used (e.g. if     *
 * in some if clause).                                                        *
 *                                                                            *
 * In C++ we would have defined an inline function instead of macros.         *
 */
#define DISKPOTENTIAL_VARS							\
  register double tmp, emp;                        /* auxiliary variables    */
#define DISKPOTENTIAL_FUNC(Z,AZ,PZ)						\
  if(Z == 0.0) {                                   /* IF(in middle of disk)  */	\
    PZ = 0.0;                                      /*   potential and force  */	\
    AZ = 0.0;                                      /*   do vanish by symmetry*/	\
  } else {                                         /* ELSE(not in disk plane)*/	\
    switch(type) {                                 /*   switch type of disk  */	\
    case thin:                                     /*   thin disk:           */	\
      if(Z > 0.0) {                                /*     above plane:       */	\
        PZ = fpot * Z;                             /*       potential and    */	\
        AZ =-facc;                                 /*       force            */	\
      } else {                                     /*     below plane:       */	\
        PZ =-fpot * Z;                             /*       potential and    */	\
        AZ = facc;                                 /*       force            */	\
      } break;                                     /*     done               */	\
    case expo:                                     /*   exponential disk:    */	\
      tmp = Z * ih;                                /*     pre-compute: z/h   */	\
      if(tmp > 0.0) {                              /*     above plane:       */	\
        emp = exp(-tmp);                           /*       exp(-z/h)        */	\
        PZ  = fpot * ( emp + tmp - 1 );            /*       potential        */	\
        AZ  =-facc * ( 1 - emp );                  /*       and force        */	\
      } else {                                     /*     below plane:       */	\
        emp = exp(tmp);                            /*       exp(-|z|/h)      */	\
        PZ  = fpot * ( emp - tmp - 1 );            /*       potential        */	\
        AZ  = facc * ( 1 - emp );                  /*       and force        */	\
      } break;                                     /*     done               */	\
    case sech:                                     /*   sech^2 disk:         */	\
      tmp = Z > 0.0 ? Z*ih : -Z*ih;                /*     precompute |z|/h   */	\
      emp = exp(-tmp);                             /*     exp(-|z|/h)        */	\
      PZ  = fpot * (0.5*tmp+log(1+emp));           /*     potential and      */	\
      if(Z > 0.0) AZ = facc*(emp-1)/(1+emp);       /*     force, depending   */	\
      else        AZ = facc*(1-emp)/(1+emp);       /*     on sign of z       */	\
      break;                                       /*   done                 */	\
    default:									\
      warning("unknown disk type in vertical disk potential");			\
    }										\
  }

#endif /* included_vertdisk_h */
