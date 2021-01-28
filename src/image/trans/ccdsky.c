/* 
 * CCDSKY: lazy sky scaling calculator to help converting theory data
 *         into degrees and m/s for FITS
 *         Can optionally convert an image. 
 *
 *      17-aug-2012     created        Peter Teuben
 *      22-aug-2012     added sdv=     PJT
 *      28-aug-2012     implemented iscale=
 *      29-jan-2013     allow distance to be in cosmological 'z'
 *      16-mar-2013     3.0 added Wright's CosmoCalculator math, verbose
 *                      (cosmocalc in ASCL?)
 *
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <extstring.h>
#include <image.h>

#include <mks.h>

string defv[] = {
  "in=\n            Optional input image file",
  "out=\n           Output image file",
  "d=1,pc\n         Distance to object, optional unit [pc]",
  "r=1,AU\n         Length scale of object, optional unit [AU]",
  "v=1,km/s\n       Velocity scale of object, optional unit [km/s]",
  "sdv=1\n          Integrated Flux (must be in Jy.km/s)",
  "scale=1\n        Scale image values",
  "H=71,0.27,0.73\n Hubble Constant, in case [d] is 'z', with optional WM and WV",
  "nsteps=1000\n    Integrations steps in the cosmo code (accuracy)",
  "VERSION=3.0a\n   16-mar-2013 PJT",
  NULL,
};

string usage = "lazy sky scaling (cosmology) calculator";

string cvsid = "$Id$";


#define CVI(x,y,z)  CubeValue(iptr,x,y,z)


#define MAXU 64


/*   
 * See also:  http://www.cv.nrao.edu/course/astr534/HILine.html
 */

static real HI_factor = 2.3595e5;

/* 
 * Xco=2e20 cm-2/K km/s, and alpha_co=4.3 Msun/K km/s 
 * also included 1.36  factor due to Helium contribution to the mass. 
 * Use Xco=14e20 for 13CO, about 7 times that of 12CO
 *
 * Older material is also: http://ned.ipac.caltech.edu/level5/March09/Solomon/Solomon2.html
 */

static real CO_factor = 1.05e4;



/* 
 * CC (cosmology) parameters
 */

static real H0;
static real WM;
static real WV;
static real DA_Mpc;

void setCC(int nc, real *cpar);
void CC(real z, int nsteps, int verbose);



/*
 *  convert a dimensionless 'z' to a distance in Mpc
 *  given a Hubble constant in km/s/Mpc
 */

real z_to_d(real z, real H)
{
  real z1 = (z+1)*(z+1);
  real d = (z1-1)/(z1+1) * c_MKS / H / 1000.0;    /*   d is now Mpc  */
  return d;
}


/*
 *     get_nu :   parse a   NUMBER,UNIT    string
 *     input:   kv         e.g.   "1,pc"
 *              defunit    e.g.   "pc"
 *     output:  value      
 *              unit
 *
 */

void get_nu(string kv, real *value, string unit, string defunit)
{
  string *sp = burststring(kv,",");
  int nsp = xstrlen(sp,sizeof(string)) - 1;
  int nv;

  *value = 1.0;
  strcpy(unit,defunit);
  if (nsp == 0) return;

  nv = nemoinpr(sp[0],value,1);
  if (nv != 1) error("error parsing %s",sp[0]);
  if (nsp == 1) return;

  strcpy(unit,sp[1]);
}

/*
 * efactor:  scaling factor between two common units
 *           u1:    first unit
 *           u2:    second unit
 *           returned value:   how many u1's in a u2  (i.e.   u2/u1)      
 *
 *           @todo need a more general unit + prefix conversion, this
 *                 was just a quick hack
 */

real efactor(string u1, string u2)
{
  real s = -1.0;
  if (streq(u1,"AU")) {
    if (streq(u2,"AU")) 
      s = 1.0;
    else if (streq(u2,"pc")) 
      s = PC/AU;
    else if (streq(u2,"km")) 
      s = 1000.0/AU;
    else if (streq(u2,"Kpc") || streq(u2,"kpc"))
      s = 1e3*PC/AU;
    else if (streq(u2,"Mpc") || streq(u2,"mpc"))
      s = 1e6*PC/AU;
    else
      error("Comparison_1 unit %s not understood for %s",u2,u1);
  } else if (streq(u1,"pc")) {
    if (streq(u2,"pc")) 
      s = 1.0;
    else if (streq(u2,"Kpc") || streq(u2,"kpc"))
      s = 1e3;
    else if (streq(u2,"Mpc") || streq(u2,"mpc"))
      s = 1e6;
    else if (streq(u2,"Gpc") || streq(u2,"Gpc"))
      s = 1e9;
    else if (streq(u2,"AU"))
      s = AU/PC;
    else if (streq(u2,"z")) {
      s = -1.0;  /* to be ammended later */
    } else
      error("Comparison_2 unit %s not understood for %s",u2,u1);
  } else if (streq(u1,"Mpc")) {
    if (streq(u2,"pc")) 
      s = 1e-6;
    else if (streq(u2,"Kpc") || streq(u2,"kpc"))
      s = 1e-3;
    else if (streq(u2,"Mpc") || streq(u2,"mpc"))
      s = 1;
    else if (streq(u2,"Gpc") || streq(u2,"gpc"))
      s = 1e3;
    else if (streq(u2,"z")) {
      s = -1.0;  /* to be ammended later */
    } else
      error("Comparison_3 unit %s not understood for %s",u2,u1);
  } else if (streq(u1,"m/s")) {
    if (streq(u2,"m/s")) 
      s = 1;
    else if  (streq(u2,"km/s")) 
      s = 1e3;
    else
      error("Comparison_4 unit %s not understood for %s",u2,u1);      
  } else
    error("Comparison_5 unit %s not understood",u1);
  dprintf(1,"u1=%s u2=%s  s=%g\n",u1,u2,s);
  return s;
}



void nemo_main(void)
{
    stream  instr, outstr;
    int     ix, iy, iz, nx, ny, nz, nc;
    int     nsteps = getiparam("nsteps");
    imageptr iptr=NULL;
    real d,r,v,  rscale, vscale, iscale = getrparam("scale");
    real mass, sdv_scale,  sdv = getrparam("sdv");
    real cpar[3];
    char ds[MAXU], rs[MAXU], vs[MAXU];

    nc = nemoinpr(getparam("H"), cpar, 3);
    if (nc < 0) error("Parsing cosmo H=%s",getparam("H"));
    setCC(nc, cpar);   /* set H0, WM, WV */

    /*  get the distance and scales      */

    get_nu(getparam("d"),&d,ds,"pc");
    get_nu(getparam("r"),&r,rs,"AU");
    get_nu(getparam("v"),&v,vs,"km/s");

    if (streq(ds,"z")) {
      CC(d,nsteps,0);
      CC(d,nsteps,1);
      printf("d=%g %s  [%g Mpc]\n",d,ds,z_to_d(d,H0));
      // d = z_to_d(d,H0);
      d = DA_Mpc;
      strcpy(ds,"Mpc");
    } else
      printf("d=%g %s\n",d,ds);
    printf("r=%g %s\n",r,rs);
    printf("v=%g %s\n",v,vs);
    printf("SdV=%g Jy.km/s\n",sdv);

    /* angles: convert to deg for FITS */
    rscale = (r * AU) / ( d * PC) * 180.0/PI;   
    rscale *= efactor("AU",rs);
    rscale /= efactor("pc",ds);
    if (rscale < 0) error("z-machine d");
    if (rscale > 0.01666) 
      printf("rscale=%g  [ %g arcsec   %g arcmin]\n",rscale,rscale*3600.0,rscale*60);
    else if (rscale < 2.777e-7)
      printf("rscale=%g  [ %g arcsec   %g mas]\n",rscale,rscale*3600.0,rscale*3600000.0);
    else
      printf("rscale=%g  (%g arcsec)\n",rscale,rscale*3600.0);
    /* velocities:  convert to m/s for FITS */
    vscale = v;
    vscale *= efactor("m/s",vs);
    printf("vscale=%g\n",vscale);

    /* mapvalues */
    printf("iscale=%g\n",iscale);

    /* sdv calculation */
    sdv_scale = efactor("Mpc",ds);
    if (sdv_scale < 0) error("z-machine sdv");
    mass = HI_factor * sqr( d * sdv_scale ) * sdv;
    printf("Mass(HI) = %g  \n",mass);
    mass = CO_factor * sqr( d * sdv_scale ) * sdv;
    printf("Mass(H2) = %g  (alpha=4.3; includes 1.36 He contribution)\n",mass);

    if (hasvalue("in") && hasvalue("out")) {      /* patch image if needed */
      instr = stropen(getparam("in"), "r");
      read_image( instr, &iptr);
      Xmin(iptr) *= rscale;   // @todo   review these when Xref != 0
      Ymin(iptr) *= rscale;
      Zmin(iptr) *= vscale;
      Dx(iptr) *= rscale;
      Dy(iptr) *= rscale;
      Dz(iptr) *= vscale;
      Beamx(iptr) *= rscale;
      Beamy(iptr) *= rscale;
      Beamz(iptr) *= vscale;
      // Xref,Yref,Zref no need to change
      if (iscale != 1.0) {
	nx = Nx(iptr);
	ny = Ny(iptr);
	nz = Nz(iptr);
	for (iz=0; iz<nz; iz++)
	  for (iy=0; iy<ny; iy++)
	    for (ix=0; ix<nx; ix++)
	      CubeValue(iptr,ix,iy,iz) *= iscale;
      }

      outstr = stropen(getparam("out"), "w");
      write_image(outstr, iptr);
    }
}



/*
 * taken from Wright's Cosmo Calculator (Schombert's python module)
 * See http://www.astro.ucla.edu/~wright

example: 

CC.py 0.1 71 0.27 0.73
  T      D kpc/arcsec  m-M
  12.38 413.49 1.82 38.29


 */

void setCC(int nc, real *cpar)
{
  if (nc==0) {              // no values, assume Benchmark Model
    H0 = 75.0;
    WM = 0.3;
    WV = 1.0 - WM - 0.4165/(H0*H0);
  } else if (nc == 1) {     // one value, assume Benchmark Model with given Ho
    H0 = cpar[0];
    WM = 0.3;
    WV = 1.0 - WM - 0.4165/(H0*H0);
  } else if (nc == 2) {     // Universe is Open, use Ho, Wm and set Wv to 0.
    H0 = cpar[0];
    WM = cpar[1];
    WV = 0.0;
  } else if (nc == 3) {     // Universe is General, use Ho, Wm and given Wv
    H0 = cpar[0];
    WM = cpar[1];
    WV = cpar[2];
  }

}

void CC(real z, int n, int verbose)
{
  int i;
  real WR, WK, Tyr, DTT, DTT_Gyr, age, age_Gyr, zage, zage_Gyr,
    DCMR, DCMR_Mpc, DCMR_Gyr, DA, DA_Gyr, kpc_DA, DL, DL_Mpc, DL_Gyr,
    V_Gpc, a, az, h, adot, ratio, x, y, DCMT, VCM, pi, c;
  
  
  /* H0, WM, WV have been set in setCC() */

  /* initialize constants */
  
  WR = 0.0;         // Omega(radiation)
  WK = 0.0;         // Omega curvaturve = 1-Omega(total)
  c = c_MKS/1000.0; // velocity of light in km/sec
  Tyr = 977.8;      // coefficent for converting 1/H into Gyr
  DTT = 0.5;        // time from z to now in units of 1/H0
  DTT_Gyr = 0.0;    // value of DTT in Gyr
  age = 0.5;        // age of Universe in units of 1/H0
  age_Gyr = 0.0;    // value of age in Gyr
  zage = 0.1;       // age of Universe at redshift z in units of 1/H0
  zage_Gyr = 0.0;   // value of zage in Gyr
  DCMR = 0.0;       // comoving radial distance in units of c/H0
  DCMR_Mpc = 0.0;  
  DCMR_Gyr = 0.0; 
  DA = 0.0;         // angular size distance
  DA_Mpc = 0.0;
  DA_Gyr = 0.0;
  kpc_DA = 0.0;
  DL = 0.0;         // luminosity distance
  DL_Mpc = 0.0;
  DL_Gyr = 0.0;     // DL in units of billions of light years
  V_Gpc = 0.0;
  a = 1.0;          // 1/(1+z), the scale factor of the Universe
  az = 0.5;         //  1/(1+z(object))

  h = H0/100;       // 
  WR = 4.165E-5/(h*h);   //includes 3 massless neutrino species, T0 = 2.72528
  WK = 1-WM-WR-WV;
  az = 1.0/(1+1.0*z);
  age = 0.;
  for (i=0; i<n; i++) {
    a = az*(i+0.5)/n;
    adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a));
    age = age + 1./adot;
  }

  zage = az*age/n;
  zage_Gyr = (Tyr/H0)*zage;
  DTT = 0.0;
  DCMR = 0.0;

  /* do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule */

  for (i=0; i<n; i++) {
    a = az+(1-az)*(i+0.5)/n;
    adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a));
    DTT = DTT + 1./adot;
    DCMR = DCMR + 1./(a*adot);
  }

  DTT = (1.-az)*DTT/n;
  DCMR = (1.-az)*DCMR/n;
  age = DTT+zage;
  age_Gyr = age*(Tyr/H0);
  DTT_Gyr = (Tyr/H0)*DTT;
  DCMR_Gyr = (Tyr/H0)*DCMR;
  DCMR_Mpc = (c/H0)*DCMR;

  /*  tangential comoving distance */

  ratio = 1.00;
  x = sqrt(abs(WK))*DCMR;
  if (x > 0.1) {
    if (WK > 0) 
      ratio =  0.5*(exp(x)-exp(-x))/x ;
    else
      ratio = sin(x)/x;
  } else {
    y = x*x;
    if (WK < 0) y = -y;
    ratio = 1. + y/6. + y*y/120.;
  }
  DCMT = ratio*DCMR;
  DA = az*DCMT;
  DA_Mpc = (c/H0)*DA;
  kpc_DA = DA_Mpc/206.264806;
  DA_Gyr = (Tyr/H0)*DA;
  DL = DA/(az*az);
  DL_Mpc = (c/H0)*DL;
  DL_Gyr = (Tyr/H0)*DL;

  /* comoving volume computation */

  ratio = 1.00;
  x = sqrt(abs(WK))*DCMR;
  if (x > 0.1) {
    if (WK > 0) 
      ratio = (0.125*(exp(2.*x)-exp(-2.*x))-x/2.)/(x*x*x/3.);
    else
      ratio = (x/2. - sin(2.*x)/4.)/(x*x*x/3.);
  } else {
    y = x*x;
    if (WK < 0) y = -y;
    ratio = 1. + y/5. + (2./105.)*y*y;
  }
  VCM = ratio*DCMR*DCMR*DCMR/3.;
  V_Gpc = FOUR_PI*qbe(0.001*c/H0)*VCM;

  if (verbose) {
    printf("-------------------------------------------------------------\n");
    printf("For H_o = %1.1f  Omega_M = %1.2f Omega_vac = %1.2f z = %1.3f\n",H0,WM,WV,z);
    printf("It is now %1.3f Gyr since the Big Bang.\n", age_Gyr);
    printf("The age at redshift z was %1.3f Gyr.\n",zage_Gyr);
    printf("The light travel time was %1.3f Gyr.\n",DTT_Gyr);
    printf("The comoving radial distance, which goes into Hubbles law, is ");
    printf("%1.1f Mpc or %1.3f Gly\n",DCMR_Mpc ,DCMR_Gyr);
    printf("The comoving volume within redshift z is %1.3f Gpc^3.\n" ,V_Gpc );
    printf("The angular size distance D_A is %1.3f Mpc or %1.3f Gly.\n",DA_Mpc, DA_Gyr);
    printf("This gives a scale of %.3f  kpc/arcsec\n" ,kpc_DA );
    printf("The luminosity distance D_L is %1.1f Mpc or %1.3f Gly.\n",DL_Mpc, DL_Gyr);
    printf("The distance modulus, m-M, is %1.2f\n" , 5*log10(DL_Mpc*1e6)-5 );
    printf("-------------------------------------------------------------\n");
  } else {
    printf("%1.2f " , zage_Gyr);
    printf("%1.2f " , DCMR_Mpc);
    printf("%1.2f " , kpc_DA);
    printf("%1.2f " , (5*log10(DL_Mpc*1e6)-5));
    printf("\n");
  }

}
