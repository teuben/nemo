/*
 * plummerv.c:  (spherical) plummer potential with variable mass/radius drive via file
 *
 *	18-jun-2003    toy model for Walter and Lia               PJT
 *	
 */

/*CTEX
 *	{\bf potname=plummerv
 *       potpars={\it $\Omega,M,R,dtcheck$}}
 *
 *  Plummer potential (BT, pp.42, eq. 2.47, see also MNRAS 71, 460 (1911))
 *
 * $$
 *    \Phi(r) = -  {  M(t)  \over
 *                    {   {(R(t)_c^2 + r^2)}^{1/2} }  }
 * $$
 */                     

  
#include <stdinc.h>
#include <filestruct.h>
#include <vectmath.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
 
local double omega = 0.0;
local double plummer_mass = 1.0;
local double plummer_radius = 1.0;
local double plummer_dt = 0.1;

local string tmr_file = NULL;

local double r2;
local double last_t = 0.0;

static void re_check(void);

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) plummer_mass = par[1];
    if (n>2) plummer_radius = par[2];
    if (n>3) plummer_dt = par[3];
    if (n>4) warning("plummer: npar=%d only 4 parameters accepted",n);

    dprintf (1,"INIPOTENTIAL Plummer: [3d version]\n");
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  mass, radius = %f %f %f \n",plummer_mass,plummer_radius,plummer_dt);
    if (name == NULL) 
      warning("plummerv: No potfile used");
    else {
      tmr_file = name;
      dprintf (1,"  potfile: %s\n",tmr_file);
    }
	
    r2 = sqr(plummer_radius);
    par[0] = omega;
}


void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    register double tmp;

    if (*time > last_t) re_check();

    tmp = 1.0/(r2 + pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    *pot = -sqrt(tmp);
    tmp *= (*pot) * plummer_mass;
    *pot *= plummer_mass;
    acc[0] = tmp*pos[0];
    acc[1] = tmp*pos[1];
    acc[2] = tmp*pos[2];
}

void potential_float (int *ndim,float *pos,float *acc,float *pot,float *time)
{
    register float tmp;

    if (*time > last_t) re_check();


    tmp = 1.0/(r2 + pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    *pot = -sqrt(tmp);
    tmp *= (*pot) * plummer_mass;
    *pot *= plummer_mass;
    acc[0] = tmp*pos[0];
    acc[1] = tmp*pos[1];
    acc[2] = tmp*pos[2];
}

static int fexist(string name) 
{
#if 0
  /* somehow, this stat64 cannot be retrieved by loadobj() */
  /* thus always make a dummy plummerv tmr file */

  struct stat buf;
  if (stat(name,&buf)==0) return 1;
  return 0;
#else
  return 1;
#endif
}

/*
 * advance time by another "dt", and check the "tmr" file
 * for new values of the mass and radius.   
   See snapcore.c for a matching writer for this reader!
 *
 */

static void re_check(void)
{
  double t1, m1, r1;
  stream f;

  last_t +=  plummer_dt;

  if (tmr_file) {
    if (fexist(tmr_file)) {
      f = stropen(tmr_file,"r");
      get_set(f,"plummerv");
      get_data_coerced(f,"Time",DoubleType,&t1,0);
      get_data_coerced(f,"Mass",DoubleType,&m1,0);
      get_data_coerced(f,"Radius",DoubleType,&r1,0);
      get_tes(f,"plummerv");
      strclose(f);
      dprintf(0,"PLUMMERV: %g %g %g\n",t1,m1,r1);

      plummer_mass = m1;
      plummer_radius = r1;
      r2 = sqr(plummer_radius);
    } else
      warning("File %s does not exist yet",tmr_file);
  }
}
