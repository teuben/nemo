/*
 * MKEXPHOT:    set up an exponential disk in a halo
 *		
 *		This is the C version of a fortran program written by
 *		Stefano Casertano, back in 198xxx
 *
 *	xx-sep-90	created from SC's fortran source
 *	16-oct-90	use default NEMO/dat
 *	22-oct-90       experiment with faking different halo fraction
 *	22-jul-91	few extra NEMO V2.x and fixed the IMSL interface
 *	23-mar-97  1.2a fix protos and SINGLEPREC  fix
 *       9-sep-01       b    gsl/xrandom
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filefn.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/barebody.h>
#include <snapshot/put_snap.c>
#include <snapshot/get_snap.c>

#include <math.h>

extern char *getmem(int);     /* memory allocator with space checker */

extern double xrandom(double,double), grandom(double,double);

/* 
        This program constructs a model galaxy with an exponential
     density distribution.  The gravitational field is read from
     a file, which has been produced by the big program.  This
     way the distribution function is closer to an equilibrium
     solution.  The gravitational field will be that due to the disk
     (which is fed in a noise-reduced form to program for the
     purpose of field computation) plus a number of non-disk
     components, which might include a flattened halo (Richstone-like)
     and other things.
        We will compute, for a number of selected radii, the radial
     force, the rotation velocity (cold & hot), the epicyclic frequency,
     and the mass-averaged vertical frequency.  The radial velocity
     dispersion will be assigned as a constant * the vertical, which
     in turn is determined by the constant-thickness condition.

        The system of units used here has  G=1,  {l} = 1 kpc, {t} = 1 Myr.
     The physical quantities output by MASSMOD have to be translated in
     this system.  The unit of masses is {m} = 2.22e11 Msolar, and for
     velocities it is {v} = 977.806 km/s.

        We need some extra arrays to store the disk model parameters
     and the radius of all particles

 */

string defv[] = {	/* DEFAULT INPUT PARAMETERS */
    "in=???\n             Input table (w/ radius,density,forces,sigmas)",
    "out=???\n		  Output file name (snapshot)",
    "nbody=1024\n	  Number of particles for disk",
    "test=\n              Input file with positions to find new velocities",
    "rmin=0.0\n	          Inner cutoff radius",
    "rmax=25.0\n          Outer cutoff radius",
    "ratio=2.0\n	  Z-R Velocity dispersion ratio (at large radii)", 
    "rcsig=2.0\n	  Z-R Velocity dispersion conversion scalelength", 
    "seed=0\n	          Random number seed",
    "headline=\n	  Extra text headline for output file",
    "fhalo=\n             Different halo/disk mass fraction from default",
    "VERSION=1.2c\n	  12-apr-04 PJT",
    NULL,
};

string usage = "set up an exponential disk in a halo";

	/* 'free' parameters, input to program by user */
local int nbody=0;
local real rmin, rmax, ratio, rcsig;
local int seed;

	/* 'fixed' parameters, read from table */
local int nrad=0;
local real rout,  h,  z0,  dens0, vhmax,  rhcore,  q2i,  dkmass,  hamass;


	/* arrays needed for the particles in the disk */
local Body *btab=NULL;
local double *rpart, *sigma1, *sigma2, *sigma3, *vrsq0, *zmax, *vrsq2, *temp;

	/* some other stuff */
local double u1[1], s1[1], u2[1], s2[1];      /* used in IMSL emulators */
local double tzero = 0.0;                     /* time of snapshot */



local double *radius , *densit;               /* used in table */
local double *fraddk , *fradha , *frad;       /* used in table */
local double *frdis0 , *frhal0 , *frtot0;     /* used in table */
local double *frdis2 , *frhal2 , *frtot2;     /* used in table */
local double *vsq , *vsq2;                    /* used in table */
local double *vrotsq , *vdisk , *vhalo;       /* used in table */
local double *vrot;                           /* used in table */
local double *sigmad , *sigmah , *sigmaz;     /* used in table */
local double *dervel , *sigmar , *sigthe;     /* used in table */
local double *sigr2 , *dersig , *qtoom;       /* used in table */
local double *velhot , *omerat , *epicyc;     /* used in table */
local double *work;                           /* used in table */


nemo_main()
{
    nbody = getiparam("nbody");
    rmin = getdparam("rmin");
    rmax = getdparam("rmax");
    ratio = getdparam("ratio");
    rcsig = getdparam("rcsig");

    read_table(getparam("in"),getparam("test"));
    set_model(getparam("test"));
    write_snap(getparam("out"), getparam("headline"));
}


read_table(tabfile,testfile)
string tabfile,testfile;
{
    stream tabstr;
    char line[256];
    int  i;
    real smolen, fhalo, fhalosqrt;
    string fn;

    sprintf(line,".:%s/data:%s/dat",getenv("NEMO"),getenv("NEMO"));
    fn = pathfind(line,tabfile);
    if (fn==NULL) error("Table %s not found",tabfile);
    tabstr = stropen(fn,"r");

/* ASCII table version */

    if (fgets(line,256,tabstr) == NULL) error("Error reading %s\n",tabfile);
    sscanf(line,"%d %lf %lf %lf %lf", &nrad, &rout, &h, &z0, &dens0);
    if (fgets(line,256,tabstr) == NULL) error("Error reading %s\n",tabfile);
    sscanf(line,"%lf %lf %lf %lf %lf", &vhmax, &rhcore, &q2i, &dkmass, &hamass);

    dprintf(0,"Model parameters read from table %s:\n",tabfile);
    dprintf(0," disk: h=%5.2f, z0=%5.2f, central density=%g, total mass=%10.7f\n",
               h, z0, dens0, dkmass); 
    dprintf(0," halo: vmax=%f core radius=%5.2f flattening=%5.2f total mass=%f\n",
               vhmax, rhcore, q2i, hamass);

    fn = getparam("fhalo");
    if (*fn) {                             /* fake a different halo/disk mass */
        fhalo = getdparam("fhalo");
        dprintf(0,"=> halo/disk mass ratio changed from %f to %f\n",
		hamass/dkmass,fhalo);
        if (fhalo <= 0.0) error("Invalid fhalo=%s",fn);
        fhalo = fhalo / (hamass/dkmass);        /* now: correction factor */
        fhalosqrt = sqrt(fhalo);
        hamass *= fhalo;
        vhmax *= fhalosqrt;
	dprintf(0,"=> halo: vmax=%f mass=%f\n",vhmax,hamass);
    } else {
        fhalo = 1.0;
        fhalosqrt = 1.0;
    }
    dprintf(0,"  (halo+disk)/disk = %f\n",(dkmass+hamass)/dkmass);

    alloc_tables(testfile);             /* allocate all the tables etc. */

    for (i=0; i<nrad; i++) {
        if (fgets(line,256,tabstr) == NULL) {
            dprintf(0,"Error reading from table %s\n",tabfile);
            dprintf(0,"Should have read %d, could only read %d entries\n",
                nrad, i-1);
            nrad = i-1;
            break;
        }
        sscanf(line," %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &radius[i], &densit[i], &frdis0[i], &frhal0[i], &frtot0[i],
            &frdis2[i], &frhal2[i], &frtot2[i], &vdisk[i], &vhalo[i], 
            &vrot[i], &sigmad[i], &sigmah[i], &sigmaz[i]);
        dprintf(2,"reading(%d): radius=%f\n",i,radius[i]);
        if (fhalo != 1.0) {
            frhal0[i] *= fhalo;
            frhal2[i] *= fhalo;
            vhalo[i] *= fhalosqrt;
            frtot0[i] = frdis0[i] + frhal0[i];
            frtot2[i] = frdis2[i] + frhal2[i];
            vrot[i] = sqrt(sqr(vdisk[i])+sqr(vhalo[i]));
        }
    }
    strclose(tabstr);
    if (nrad < 3) error("Table file %s too small\n",tabfile);
    dprintf(0,"Table has %d radii between radius = %f and %f\n",
             nrad, radius[0], radius[nrad-1]);

/*
 *       Before doing anything else with these, we will smooth
 *    the data (essential, as later we shall differentiate them!)
 */
      smolen = 0.5;
      smooth (radius, sigmaz, nrad, smolen, work);
      /* smooth (radius, frtot0, nrad, smolen, work); */
      /* smooth (radius, frtot2, nrad, smolen, work); */
}

alloc_tables(testfile)
string testfile;
{
    int   bits, nbpr = sizeof(real);
    stream teststr;

    radius = (double *) getmem(nbpr * nrad);
    densit = (double *) getmem(nbpr * nrad);
    fraddk = (double *) getmem(nbpr * nrad);
    fradha = (double *) getmem(nbpr * nrad);
    frad   = (double *) getmem(nbpr * nrad);
    frdis0 = (double *) getmem(nbpr * nrad);
    frhal0 = (double *) getmem(nbpr * nrad);
    frtot0 = (double *) getmem(nbpr * nrad);
    frdis2 = (double *) getmem(nbpr * nrad);
    frhal2 = (double *) getmem(nbpr * nrad);
    frtot2 = (double *) getmem(nbpr * nrad);
    vsq    = (double *) getmem(nbpr * nrad);
    vsq2   = (double *) getmem(nbpr * nrad);
    vrotsq = (double *) getmem(nbpr * nrad);
    vdisk  = (double *) getmem(nbpr * nrad);
    vhalo  = (double *) getmem(nbpr * nrad);
    vrot   = (double *) getmem(nbpr * nrad);
    sigmad = (double *) getmem(nbpr * nrad);
    sigmah = (double *) getmem(nbpr * nrad);
    sigmaz = (double *) getmem(nbpr * nrad);
    dervel = (double *) getmem(nbpr * nrad);
    sigmar = (double *) getmem(nbpr * nrad);
    sigthe = (double *) getmem(nbpr * nrad);
    sigr2  = (double *) getmem(nbpr * nrad);
    dersig = (double *) getmem(nbpr * nrad);
    qtoom  = (double *) getmem(nbpr * nrad);
    velhot = (double *) getmem(nbpr * nrad);
    omerat = (double *) getmem(nbpr * nrad);
    epicyc = (double *) getmem(nbpr * nrad);
    work   = (double *) getmem(nbpr * nrad * 3);

    if (*testfile) {
        teststr = stropen(testfile,"r");
	get_history(teststr);
        get_snap(teststr,&btab,&nbody,&tzero,&bits);
        strclose(teststr);
        dprintf(0,"[Read %d bodies from testfile %s]\n",nbody,testfile);
    } else
        btab = (Body *) getmem(sizeof(Body) * nbody);


    rpart  = (double *) getmem(nbpr * nbody);
    sigma1 = (double *) getmem(nbpr * nbody);
    sigma2 = (double *) getmem(nbpr * nbody);
    sigma3 = (double *) getmem(nbpr * nbody);
    vrsq0  = (double *) getmem(nbpr * nbody);
    zmax   = (double *) getmem(nbpr * nbody);
    vrsq2  = (double *) getmem(nbpr * nbody);
    temp   = (double *) getmem(nbpr * nbody);

    
}

write_snap(name, headline)
string name;
string headline;
{
    stream outstr;
    int bits = MassBit | PhaseSpaceBit;

    if (! streq(headline, ""))
	set_headline(headline);
    outstr = stropen(name, "w");
    put_history(outstr);
    put_snap(outstr, &btab, &nbody, &tzero, &bits);
    strclose(outstr);
}

set_model(testfile)
string testfile;
{
    int    i, m1, m2, ier, n, ie, nbad1, nbad2;
    double a1, a2, e1, e2, dlsigz, r1sig,r,fff, vvcorr, othcor;
    double dum, con, yy, ey, yy1, del, ang, cang, sang, yran;		/* dum not used */
    double vr, vt, vz, zmaxsq, squvel, vcirc;
    Body   *bp;

/*
        The file called holds all the basic information we
     need on the mass model, namely the z-averaged radial and vertical
     force.  It also has already the values of the disk mass, etcetera.
     The disk mass will be corrected for the finite extent of the disk.
*/
      seed = init_xrandom(getparam("seed"));

      a1 = rmin / h;
      a2 = rmax / h;
      e1 = exp( -a1) * (1.0 + a1);
      e2 = exp( -a2) * (1.0 + a2);
      dkmass = dkmass * (e1 - e2);
/*
        Now set up for the interpolation.  We want to calculate
     the circular velocity (for the hot disk) and the velocity
     dispersion at any radius, so that we can store the result
     in the arrays  VRSQ,  VRSQ2,  SIGMA.  This will be done by interpolating
     in the relevant arrays.

        The circular velocity for the hot disk is given by the formula

     velhot^2 = vtotal^2 + vvcorr + (sigmar^2 - sigthe^2)

        with

     vvcorr = sigmar^2 * d ln(mu*sigmar^2) / d ln(r)

        We will obtain  SIGMAR  and its derivative by interpolation, whereas
     we will just use  -r/h  as the logarithmic derivative of the density.
        Therefore:

     vvcorr = r * [ d (sigmar^2) /dr  -  sigmar^2 / h ]

     SIGMAR  (square stored in SIGR2) is defined as a constant *
     the vertical velocity dispersion
     VZDISP,  which in turn is fixed from the thickness and surface density.
     However, the constant is lowered in the innermost parts, to avoid
     the rotation velocity to be imaginary.  The actual form is

          sigmar = sigmaz * ratio * (r + r1)/(r + rc)

     Here,  RC  is a free input parameter, and  R1  is determined so
     that the  SIGMAR  has zero derivative in the center.  This gives
     the condition:
           r1 = rc / (1 - rc * d(ln(sigmaz))/dr)

     which would translate into   r1 = 2 h rc / (2 h + rc)  for the
     `theoretical' scaling   SIGMAZ \propto exp(-r/2h)
     
     Therefore we first calculate the logarithmic derivative of  SIGMAZ
     in the center:
 */  

      icsccu (radius, sigmaz, nrad, work, nrad, ier);   /* was 'ier' set ??? */
      m1 = 1;
      u1[0] = 0.0;
      m2 = 0;
      dcsevu (radius, sigmaz, nrad, work, nrad, u1, s1, m1, s2, m2, ier);
      dlsigz = s1[0] / sigmaz[0];
/*
 *     This assumes that the first radius is 0!
 */
      r1sig = rcsig / (1. - rcsig * dlsigz);
      dprintf(0,"At r=0 ratio of radial to vertical vel.disp. = %f\n",
                ratio*r1sig/rcsig);
/*
 *        Calculate VROTSQ, the rotation velocity at  z=z0
 */
      for (i=0; i<nrad; i++)
        vrotsq[i] = -radius[i] * (frtot0[i] + 0.5 * sqr(z0) * frtot2[i]);
/*
 *       Calculate DERVEL  (actually 2*dervel; the division by two
 *    is done later
 */
      icsccu (radius, vrotsq, nrad, work, nrad, ier);
      m2 = 0;
      dcsevu (radius, vrotsq, nrad, work, nrad, radius,
                                      dervel, nrad, s2, m2, ier);

      nbad1 = 0;
      for (i=0; i<nrad; i++) {
            r = radius[i];
/* 
 *        First correct DERVEL (was too big by a factor two), and
 *                    calculate OMERAT (2 Omega/kappa) and EPICYC
 */
            if (vrotsq[i] > 0.0) {
                  vrot[i] = sqrt (vrotsq[i]);
                  dervel[i] = 0.5 * radius[i] * dervel[i] / vrotsq[i];
            } else {
                  vrot[i] = 0.0;
                  dervel[i] = 1.0;
                  nbad1++;
		  dprintf(1,"i=%d: vrotsq = %f at radius %f\n",
				i,vrotsq[i],radius[i]);
            }

            omerat[i] = sqrt( 2.0 / (1.0 + dervel[i]) );
            epicyc[i] = 0.0;
            if (r > 0.0) epicyc[i] = 2.0 * vrot[i] / r / omerat[i];
/*
 *       Then calculate the velocity dispersions
 */
            fff = (r + r1sig) / (r + rcsig);
            sigmar[i] = fff * ratio * sigmaz[i];
            sigthe[i] = sigmar[i] / omerat[i];
            sigr2[i] = sqr(sigmar[i]);
/*
 *       And now Toomre's  Q
 */
            qtoom[i] = 0.0;
            densit[i] = dens0 * exp(-r/h);
            if (densit[i] > dens0/300.0)
                  qtoom[i] = sigmar[i] * epicyc[i] /
                                  (3.360 * densit[i]);   /* G==1 */
      }
      if (nbad1) dprintf(0,"There were %d radii with vrotsq <= 0\n",nbad1);
/*
 *          Correct the `cold' rotation curve
 */
      m2 = 0;

      icsccu (radius, sigr2, nrad, work, nrad, ier);
      dcsevu (radius, sigr2, nrad, work, nrad, radius,
          dersig, nrad, s2, m2, ier);
/*
 *
 */
      nbad2=0;
      for (i=0; i<nrad; i++) {
          r = radius[i];
          vvcorr = r * (dersig[i] - sigr2[i] / h);
          othcor = sqr(sigmar[i]) - sqr(sigthe[i]);

          vsq[i] = -r * frtot0[i] + vvcorr + othcor;
          vsq2[i] = -r * frtot2[i];

          if (vsq[i] >= 0.0)
              velhot[i] = sqrt( vsq[i] );
          else {
              velhot[i] = 0.0;
              nbad2++;
              dprintf(1,"i=%d: vsq = %f at radius %f\n",
                        i, vsq[i], radius[i]);
	  }

          dprintf(3,"AAA %f %f %f %f %f %f %f %f\n", r, r*frtot0[i],
               sigr2[i], dersig[i], vvcorr, othcor, vsq[i], velhot[i]);
      }
      if (nbad2) dprintf(0,"There were %d radii with vsq <= 0\n",nbad2);

      for (i=0; i<nrad; i++)
          dprintf(2,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                radius[i], densit[i], epicyc[i], vdisk[i],
                vhalo[i], vrot[i], sigmad[i], sigmah[i], sigmaz[i],
                sigmar[i], sigthe[i], velhot[i], omerat[i],
                dervel[i], qtoom[i]);

/*
 *        Now the vectors with the model properties have been filled.
 */


      if (*testfile) {		/* get radii from old file */
          for (i=0, bp=btab; i<nbody; i++, bp++)
              rpart[i] = sqrt(sqr(Pos(bp)[0]) + sqr(Pos(bp)[1]));
      } else {			/* generate radii between rmin and rmax */
          a1 = rmin / h;
          a2 = rmax / h;
          e1 = exp( -a1) * (1.0 + a1);
          e2 = exp( -a2) * (1.0 + a2);

          for (i=0; i<nbody; i++) {
              con = xrandom(e1, e2);	/* uniform in mass-space */
              yy = 1.0;			/* initial guess r for M(r)=con */
              do {			/* loop until r found */
                  ey = exp( -yy);
                  yy1 = yy + (ey * (1.0 + yy) - con) / (yy * ey);
                  del = ABS(yy - yy1);
                  yy = yy1;
              } while (del > 1.0e-5);
              if(n == 2) yy = yy + 0.02;   /* was n set ??? */
              if(n == 3) yy = yy - 0.01;
              rpart[i] = h * yy;
          }
/*
 *        set up positions
 */
          for (i=0, bp=btab; i<nbody; i++, bp++) {
              r = rpart[i];
              ang = TWO_PI * xrandom(0.0,1.0);
              cang = cos(ang);
              sang = sin(ang);
              Pos(bp)[0] = r * cang;
              Pos(bp)[1] = r * sang;

              if (z0 > 0.0) {
                  yran = xrandom(-1.0, 1.0);
		  Pos(bp)[2] = 0.5*z0 * log10((1.0 + yran) / (1.0 - yran));
              } else
                  Pos(bp)[2] = 0.0;
          }
      } /* testfile */

/*
 *
 *      Now we know the radii for all particles.  Let us fill in
 *   the vectors with velocity dispersions and circular velocity
 *    at each radius
 */

/*
 *   1) square of circular velocity, orders 0 and 2  (VRSQ0, VRSQ2)
 */
      icsccu(radius, vsq, nrad, work, nrad, ie);  /* was 'ie' set ?? */
      icsevu(radius, vsq, nrad, work, nrad,
	                 rpart, vrsq0, nbody, ie);


      icsccu(radius, vsq2, nrad, work, nrad, ie);
      icsevu(radius, vsq2, nrad, work, nrad,
        rpart, vrsq2, nbody, ie);
/*
 *     2) radial velocity dispersion (SIGMA1)
 */
      icsccu (radius, sigmar, nrad, work, nrad, ie);
      icsevu (radius, sigmar, nrad, work, nrad,
          rpart, sigma1, nbody, ie);

/*
 *     3) tangential velocity dispersion (SIGMA2)
 */
      icsccu (radius, sigthe, nrad, work, nrad, ie);
      icsevu (radius, sigthe, nrad, work, nrad,
        rpart, sigma2, nbody, ie);
/*
 *    4) vertical velocity dispersion (SIGMA3)
 */
      icsccu (radius, sigmaz, nrad, work, nrad, ie);
      icsevu (radius, sigmaz, nrad, work, nrad,
         rpart, sigma3, nbody, ie);

/*
 *     This should do it!
 *          now, velocities
 */
       for (i=0, bp=btab; i<nbody; i++, bp++) {
            vr = grandom(0.0,sigma1[i]);
            vt = grandom(0.0,sigma2[i]);
            if (z0 > 0.0) 
                  vz = grandom(0.0,sigma3[i]);
            else
                  vz = 0.0;
            if (rpart[i] > 0.0) {
                  cang = Pos(bp)[0] / rpart[i];
                  sang = Pos(bp)[1] / rpart[i];
            } else {
                  cang = 0.0;
                  sang = 0.0;
            }
/*
 *      correction for vertical energy in the rotation velocity
 */
            zmaxsq = sqr(Pos(bp)[2]) + 0.5 * sqr(vz * z0 / sigma3[i]);;
            squvel = vrsq0[i] + 0.25 * vrsq2[i] * zmaxsq;
            vcirc = sqrt (MAX (0.0, squvel));
            Vel(bp)[0] = -(vcirc + vt) * sang + vr * cang;
            Vel(bp)[1] = (vcirc + vt) * cang + vr * sang;
            Vel(bp)[2] = vz;
            Mass(bp) = dkmass / nbody;	/* all bodies equal mass */
      }
}

smooth (x, y, n, smolen, work)
int n;
real *x, *y, *work, smolen;
{
    int i, j;
    real weight, r1, r2, xxx, f;	/* r2 not used */

    for (i=0; i<n; i++) {
        r1 = x[i];
        work[i] = 0.0;
        weight = 0.0;

        for (j=0; j<n; j++) {
            xxx = 0.50 * sqr((x[j]-r1)/smolen);
            if (xxx < 10) {
                f = exp (-xxx);
                weight  += f;
                work[i] += f * y[j];
            }
        }
        if (weight > 0.0) work[i] /= weight;
    }
    for (i=0; i<n; i++)
        y[i] = work[i];
}
