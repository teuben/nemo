/* 
 * CCDVEL: create a velocity field: either grid or retrace
 *
 *  todo
 *      - general epicyclic orbits (re-introduce 'q' and 'alpha'
 *      - check 1-based and 0-based versions of rotcur/tables
 *        for the center position
 *      - allow vz motions
 *
 *  One could either invert the sky and find out in which ring
 *  that pixel lies, and compute the radial velocity, which means
 *  for warps that given a certain width one could get holes in
 *  the coverage, (GIPSY mode)
 *      or,
 *  one can fake a more or less uniform sampling of the disk plane,
 *  and project it. This means *two* images are allocated, one for
 *  velocities, the other to count the number of occurences per cell.
 *  (NEMO mode)
 *
 *  Standard bench on sparc1:
 *      time ccdvel junk1 0:50 0:100:10 inc=45 pa=30::10,30:60
 *
 * Timing: 23.2" original GIPSY program (C version)
 *         16.6" after radius() is called less, break loop when ring found
 *         10.5" after storing ring sinp, cosp, cosi in lookup table
 *          5.6" after deleteing sqrt() in radius; now callee's need repair
 *               since their interpolation is not good anymore
 *          6.3" correct interpolation using 'e/(r+sqrt(r*r-e))' instead of 'e'
 *
 *  
 *
 *  18-may-91   V1.0  Created                                           PJT 
 *  12-jun-92   V1.2  rotcurfit= to be output table from rotcur         PJT
 *  31-jul-92   V1.3  rings have now become pure radii, like rotcurfit= PJT
 *                    removed q= and alpha= as keywords for now
 *   1-aug-92   V1.3a optimized GIPSY mode (_1)
 *   3-aug-92       b fixed size bug, override center= allowed          PJT
 *		    c another shortcut to gain speed
 *  12-aug-92       d added headline=                           PJT
 *                  e write Dx/Dy into image (was never done)   PJT
 *  22-nov-94       f fixed -DSINGLEPREC declarations    	PJT
 *  11-jan-95       g ooops, bug in cosi[0] as cosi[n]		pjt
 *  29-aug-95   V1.4  added toy spiral perturbations            pjt
 *  22-feb-97	    a made sure it worked with SINGLEPREC	pjt
 *		      and kept variables local
 *  12-apr-01   V1.5  make some crude (no neighbor voxel pollution) 
 *                    cubes optionally
 *  20-apr-01       a cube creation not worked too well yet, but fixed
 *                    some backwards compatibility problems
 *   1-may-01   V1.6  Added amp= keyword for fudzing density maps    PJT/LGM
 *   3-may-01   V1.7  Added vexp= keyword                            PJT
 *   8-sep-01      a  init_xrandom
 *  29-jun-02   V1.8  add blc= for "BLC"                            PJT
 *  20-jul-02      a  fix error message if too many array elements   pjt
 *   6-feb-03   V2.0  add phases and fourier components to vel's for Rahul    pjt
 *  11-sep-03   V2.1  add intpol=    PJT
 *  30-nov-20   V2.2  fixes when in= template given; cleanup for ansi     PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <image.h>
#include <table.h>

string defv[] = {
        "out=???\n      Output file name (an image)",
        "radii=\n       Radii in disk",
        "vrot=\n        Rotation velocities (or intensity if amp=1)",
	"vexp=\n        Expansion velocities",
        "inc=\n         Inclinations",
        "pa=\n          Position angle: C.C.W. from +Y-axis",
#if 0
        /* to be implemented later e.g. Bosma's thesis 1978 method */
        "q=1\n          Axis ratio of orbit (**not implemented**)",
        "alpha=1\n      Axis ratio of epicycle (1<alpha<2)",
#endif

        "aspiral=0\n    Amplitude of spiral perurbation",
        "pspiral=45\n   phase of vr/vt perturbations",
        "kspiral=1\n    spiral wavenumber, inverse pitch angle",
        "nspiral=2\n    measure of inverse width of spiral",
	"tspiral=0\n	phase of spiral w.r.t. line of nodes", 

	"phrot=0\n      phases of vrot fourier component (per ring)",
	"phexp=0\n      phases of vexp fourier component (per ring)",
	"mrot=0\n       mode of fourier component for vrot (one number)",
	"mexp=0\n       mode of fourier component for vexp (one number)",

        "size=128\n     Map/Cube size",
        "cell=1\n       Cell size",
        "center=\n      Rotation center of disk [default: map center]",
	"blc=0\n        Set bottom-left corner (normally 0,0)",
        "vsys=\n        Systemic velocity",
        "rotcurfit=\n   initial conditions from rotcur output table",
        "fixring=1\n    ring to use for center, and vsys",
        "noise=0\n      Gaussian noise added to velocities",
	"amp=f\n        Create density maps instead of velocity?",
        "seed=0\n       Initial random seed",
	"in=\n          Template 2D image for cube generation",
	"intpol=linear\n Interpolation: linear or constant",
	"headline=\n	Optional random verbiage",
        "VERSION=2.2\n  30-nov-2020 PJT",
        NULL,
};

string usage="create a velocity field from a rotation curve";

#define RPD  PI/180
#ifndef HUGE
#define HUGE 1.0e20
#endif

local stream  outstr;                    /* output file */

local int  size[3];                /* 2D/3D size of map/cube   */
local real cell[3];                /* cell sizes */
local real center[3];              /* rot center */
local real blc[3];                 /* bottom left corner */
local real rmin, rmax;             /* spatial extent of disk  */
local real vsys = 0.0;             /* systemic velocity */
local real noise;                  /* some noise added to the velocities */
local int  out_mode;               /* output mode       */
local imageptr iptr;               /* template image */

                             /* 1 = sky -> gal */
                             /* 2 = gal -> sky */
local real aspiral, pspiral, kspiral, nspiral;    /* spiral parameters */
local real vrotfac, vexpfac, theta0;

local real Qamp;   /* use amplitudes instead of velocities */

local real undef = 0.0;            /* could also use IEEE NaN ??? set_fblank ? */
local bool Qcube = FALSE;
local bool Qlinear = TRUE;         /* really: intpol=linear or const */

#define MAXCOL    7     /* columns to read from rotcur table file */
#if !defined(MAXRAD)
#define MAXRAD 1024     /* radii to store info of                 */
#endif

/* arrays which store info from table/key: */
local real rad_i[MAXRAD],                                      /* radii */
           vrot_i[MAXRAD], vexp_i[MAXRAD],                /* velocities */
           inc_i[MAXRAD], theta_i[MAXRAD],                   /* viewing */
           xpos_i[MAXRAD], ypos_i[MAXRAD], vsys_i[MAXRAD],    /* center */
           phrot_i[MAXRAD], phexp_i[MAXRAD];                  /* phases */
local int  mrot_i, mexp_i;                            /* fourier moment */
local int  nrad;

local int gridx(real), gridy(real), gridz(real);
local void vel_create_1(stream), vel_create_2(stream);
local real radius(real, real , real , real , real, real);
local real linfit(real *, real *, real, int, int);
local real my_image(real, real);

extern double grandom(double,double);

void setparams(void);


void nemo_main(void)
{
    setparams();
    outstr = stropen(getparam("out"),"w");
    switch (out_mode) {
       case 1:  vel_create_1(outstr);  break;
       case 2:  vel_create_2(outstr);  break;
       default: error("No valid mode");
    }
    strclose(outstr);
}

void setparams(void)
{
    int i,n, ifix, colnr[MAXCOL];
    real *coldat[MAXCOL];
    stream tstr;
    string intpol;

    Qcube = hasvalue("in");
 
    out_mode = 1;       /* fix to 'sky -> gal plane' mode */

    switch ( nemoinpi(getparam("size"),size,3) ) {  /* MAPSIZE */
       case 1:                  /*  nx[,ny] */
                size[1] = size[0];
		size[2] = 1;
                break;
       case 2:                  /*  nx,ny */
		size[2] = 1;
                break;
       case 3:                  /*  nx,ny */
                break;
       default:                 /* --- some error --- */
                error("Need size= keyword");
    }

    switch ( nemoinpr(getparam("cell"),cell,3) ) {   /* CELL SIZE */
       case 1:                  /*  dx[,dy] */
                cell[2] = cell[1] = cell[0];
                break;
       case 2:                  /*  dx,dy */
                break;
       case 3:                  /*  dx,dy */
                break;
       default:                 /* --- some error --- */
                error("Need cell= keyword");
    }

    switch ( nemoinpr(getparam("center"),center,2) ) {   /* ROT CENTER */
       case 1:                  /*  x0[,y0] */
                center[1] = center[0];
                break;
       case 2:                  /*  x0,y0 */
                break;
       case 0:                  /*  [0,0,0] */
                center[0] = 0.5*(size[0]-1.0);
                center[1] = 0.5*(size[1]-1.0);
                break;
       default:                 /* --- some error --- */
                error("Syntax error in offset= keyword");
    }

    switch ( nemoinpr(getparam("blc"),blc,2) ) {         /* BLC edge */
       case 1:                  /*  xmin[,ymin] */
                blc[1] = blc[0];
                break;
       case 2:                  /*  xmin,ymin */
                break;
       case 0:                  /*  [0,0,0] */
	        blc[0] = blc[0] = 0.0;
                break;
       default:                 /* --- some error --- */
                error("Syntax error in blc= keyword");
    }


    if (hasvalue("rotcurfit")) {  /* if table given: override all other */

        /* set column numbers and pointers to data to be received */
        colnr[0] = 1;     coldat[0] = rad_i;            /* radius */
        colnr[1] = 2;     coldat[1] = vsys_i;           /* systemic vel */
        colnr[2] = 4;     coldat[2] = vrot_i;           /* rotation vel */
        colnr[3] = 6;     coldat[3] = theta_i;          /* position angle */
        colnr[4] = 8;     coldat[4] = inc_i;            /* inclination */
        colnr[5] = 10;    coldat[5] = xpos_i;           /* x center */
        colnr[6] = 12;    coldat[6] = ypos_i;           /* y center */
        /* get data from table */
        nrad = get_atable(stropen(getparam("rotcurfit"),"r"),
                          7, colnr, coldat, MAXRAD);

        if (hasvalue("radii")) warning("Cannot override radii");
        if (hasvalue("vrot")) {
            n = nemoinpr(getparam("vrot"),vrot_i,MAXRAD);
            if (n<1) error("vrot=: Need at least one (%d) ",n);
            dprintf(0,"Overriding %d velocitie(s); setting remaining %d to %g\n",
                    n,nrad-n,vrot_i[n-1]);
            for (i=n; i<nrad; i++) vrot_i[i] = vrot_i[n-1];
        }
        if (hasvalue("inc")) {
           n = nemoinpr(getparam("inc"),inc_i,MAXRAD);
           if (n<1) error("inc=: Need at least one (%d)",n);
           dprintf(0,"Overriding %d inclinations\n",n);
           for (i=n; i<nrad; i++) inc_i[i] = inc_i[n-1];
        }
        if (hasvalue("pa")) {
           n = nemoinpr(getparam("pa"),theta_i,nrad);
           if (n<1) error("pa=: Need at least one (%d)",n);
           dprintf(0,"Overriding %d position angles\n",n);
           for (i=n; i<nrad; i++) theta_i[i] = theta_i[n-1];
        }

        for (i=0; i<nrad; i++) {        /* store angles in radians */
            theta_i[i] *= RPD;
            inc_i[i]   *= RPD;
        }

        ifix = getiparam("fixring") - 1;
        if (ifix < 0 || ifix >= nrad) {
    	    warning("ReSetting fixring (%d) to the first ring",ifix);
            ifix = 0;
        }
        if (hasvalue("vsys"))           /* get vsys from table, or key */
            vsys = getdparam("vsys");
        else
            vsys = vsys_i[ifix];
        
	if (hasvalue("center")) {    
            if (nemoinpr(getparam("center"),center,2) != 2)
                error("center= needs two values if you want to override");
        } else {
            center[0] = xpos_i[ifix];       /* always get center from table ? */
            center[1] = ypos_i[ifix];
	}
        
    } else {

        nrad = nemoinpr(getparam("radii"),rad_i,MAXRAD);
        if (nrad<2) error("radii=: Need at least two radii (%d)",nrad);
        dprintf(0,"Found %d radii\n",nrad);

        n = nemoinpr(getparam("vrot"),vrot_i,nrad);
        if (n<1) error("vrot=: Need at least one (%d) or at most %d",n,nrad);
        dprintf(0,"Found %d rot velocities\n",n);
        for (i=n; i<nrad; i++) vrot_i[i] = vrot_i[n-1];

        n = nemoinpr(getparam("vexp"),vexp_i,nrad);
        if (n<0) error("vexp=: Need at least one (%d) or at most %d",n,nrad);
	if (n==0) {
	  vexp_i[0] = 0.0;
	  n = 1;
	} else
	  dprintf(0,"Found %d exp velocities\n",n);
        for (i=n; i<nrad; i++) vexp_i[i] = vexp_i[n-1];

        n = nemoinpr(getparam("inc"),inc_i,nrad);
        if (n<1) error("inc=: Need at least one (%d) or at most %d",n,nrad);
        dprintf(0,"Found %d inclinations\n",n);
        for (i=n; i<nrad; i++) inc_i[i] = inc_i[n-1];

        n = nemoinpr(getparam("pa"),theta_i,nrad);
        if (n<1) error("pa=: Need at least one (%d) or at most %d",n,nrad);
        dprintf(0,"Found %d position angles\n",n);
        for (i=n; i<nrad; i++) theta_i[i] = theta_i[n-1];



        n = nemoinpr(getparam("phrot"),phrot_i,nrad);
        if (n<1) error("phrot=: Need at least one (%d) or at most %d",n,nrad);
        dprintf(0,"Found %d phase angles for rotation\n",n);
        for (i=n; i<nrad; i++) phrot_i[i] = phrot_i[n-1];

        n = nemoinpr(getparam("phexp"),phexp_i,nrad);
        if (n<1) error("phexp=: Need at least one (%d) or at most %d",n,nrad);
        dprintf(0,"Found %d phase angles for expansion\n",n);
        for (i=n; i<nrad; i++) phexp_i[i] = phexp_i[n-1];

        mrot_i = getiparam("mrot");
	mexp_i = getiparam("mexp");

	if (mrot_i || mexp_i)
	  warning("You are using a new option mrot= and/or mexp=; be aware");

        vsys = getdparam("vsys");
   
        for (i=0; i<nrad; i++) {     /* store angles in radians */
            theta_i[i] *= RPD;
            inc_i[i]   *= RPD;
	    phrot_i[i] *= RPD;
	    phexp_i[i] *= RPD;
        }
    }

    if (Qcube) {
      warning("New experimental CUBE mode");
      tstr = stropen(getparam("in"),"r");
      read_image(tstr,&iptr);
      strclose(tstr);

      size[0] = Nx(iptr);
      size[1] = Ny(iptr);
      size[2] = 1;

      cell[0] = -Dx(iptr);
      cell[1] = Dy(iptr);
      cell[2] = Dz(iptr);
      Qcube = FALSE;
    }
    


    noise = getdparam("noise");
    init_xrandom(getparam("seed"));
    aspiral = getdparam("aspiral");
    pspiral = getdparam("pspiral");
    kspiral = -getdparam("kspiral");	/* positive = trailing arms */
    nspiral = getdparam("nspiral");
    vrotfac = aspiral * cos(pspiral*RPD);
    vexpfac = aspiral * sin(pspiral*RPD);
    theta0 = getdparam("tspiral") * RPD;

    /* set for convenience */

    rmin = rad_i[0];
    rmax = rad_i[nrad-1];

    /* output to user */
    printf("Mapsize: %d * %d * %d\n",size[0],size[1],size[2]);
    printf("Cell size: %g * %g * %g\n",cell[0],cell[1],cell[2]);
    printf("Center: %g * %g pixels\n",center[0],center[1]);
    printf("BLC: %g %g\n",blc[0],blc[1]);
    printf("Systemic velocity: %g\n",vsys);
    
    for (n=1; n<nrad; n++)
      if (rad_i[n] < rad_i[n-1]) error("Radii not sorted (@%d: %g < %g)",
				       n,rad_i[n] < rad_i[n-1]);

    if (hasvalue("headline")) set_headline(getparam("headline"));

    Qamp = getbparam("amp");
    
    intpol = getparam("intpol");
    if (*intpol == 'l')
      Qlinear = TRUE;
    else if (*intpol == 'c')
      Qlinear = FALSE;
    else
      error("intpol=%s should be l(inear) or c(onstant)",intpol);
} /* setparams */


/*
 * create a velocity field (GIPSY method)
 *              0..nx-1 and 0..ny-1 
 *      start from pixel, work back to gal plane and interpolate
 *      (a.k.a. retracing method)
 *
 */
local void vel_create_1(stream outstr)
{
    int  i, j, k, n, nx, ny, nz, ir, nr, ip, np;
    real x, y, t, vrot, vexp, ival, fexp, frot, phexp, phrot;
    real cost, sint, cosk, sink, rad, drad, phi, dphi, sinth, costh;
    real omega, kappa, bigx, bigy, p, r, m_min, m_max, delta1, delta2;
    real x0, y0, z0, dx, dy, dz, eps, theta, inc, e1, e2, vspi, vobs;
    imageptr vptr;
    real e[MAXRAD], sinp[MAXRAD], cosp[MAXRAD], cosi[MAXRAD];

    dprintf(0,"gridding mode: sky->gal plane\n");

    m_min = HUGE; m_max = -HUGE;
    nx = size[0];
    ny = size[1];
    nz = size[2];
    dx = cell[0];   /* ??? perhaps -cell[0] = astro. convention */
    dy = cell[1];
    dz = cell[2];
    x0 = center[0];
    y0 = center[1];
    z0 = center[2];
    eps = 0.1 * sqrt(dx*dx+dy*dy);       /* some small fraction of a pixel */

    if (Qcube) {
      if (!create_cube(&vptr, nx, ny, nz))   /* image cube */
        error("Could not create cube from scratch");
    } else {
      if (!create_image(&vptr, nx, ny))   /* velocity image */
        error("Could not create image from scratch");
    }

    for (n=0; n<nrad; n++) {
        sinp[n] = sin(theta_i[n]);
        cosp[n] = cos(theta_i[n]);
        cosi[n] = cos(inc_i[n]);
    }

    for (j=0; j<ny; j++) {          /* Loop over all pixels */
        y = dy*(j-y0) + blc[1];
        for (i=0; i<nx; i++) {      
	    x = dx*(i-x0) + blc[0];        /* (x,y) are w.r.t. center */
	    for (k=0; k<nz; k++)
	      CubeValue(vptr,i,j,k) = undef;       /* set to 'undefined' */

            r = sqrt(sqr(x)+sqr(y));        /* get projected radius on sky */
	    // dprintf(0,"PJT %d %d %g %g  %g %g\n",i,j,x,y,r,rmax);
            if (r > rmax) continue;          /*  certainly outside disk */

            delta1 = e[0] = radius(rad_i[0], sinp[0], cosp[0], cosi[0], x, y);
            n = nrad-1;
            delta2 = e[n] = radius(rad_i[n], sinp[n], cosp[n], cosi[n], x, y);
            if (delta1*delta2 > 0) continue;        /* not in disk at all ? */

	    dprintf(1,"PJT %d %d %g %g  %g %g   %g %g\n",i,j,x,y,r,rmax, delta1, delta2);

            for (n=1, nr=-1; n<nrad; n++) {    /* find in which 'ring' (x,y) is */
                e[n] = radius(rad_i[n], sinp[n], cosp[n], cosi[n], x, y);
                if (e[n-1]*e[n] <= 0.0) {   /* found it, interpolate to rad */
                    e1 = e[n-1]/(rad_i[n-1]+sqrt(sqr(rad_i[n-1])-e[n-1]));
                    e2 = e[n]  /(rad_i[n]  +sqrt(sqr(rad_i[n])  -e[n]));
                    rad = (rad_i[n]*e1 - rad_i[n-1]*e2)/(e1-e2);
                    nr = n;
                    break;    /* 'got it: rad' is the innermost galactic radius */
                }
            }   
            if (nr<0) {     /* should not happen? perhaps with strange warps? */
                warning("? Warp programming mistake at cell (%d,%d) ?",i,j);
                continue;
            }
            if (rad > eps) {
                theta = linfit(theta_i, rad_i, rad, nr-1, nr);
                inc   = linfit(inc_i,   rad_i, rad, nr-1, nr);
                vrot  = linfit(vrot_i,  rad_i, rad, nr-1, nr);
                vexp  = linfit(vexp_i,  rad_i, rad, nr-1, nr);
		sinth = sin(theta);
                costh = cos(theta);
                cost = (-x*sinth + y*costh) / rad;
                sint = (-x*costh - y*sinth) / (rad * cos(inc));
		if (mrot_i || mexp_i) {         /* modulated orbit */
		  phrot = mrot_i ? linfit(phrot_i, rad_i, rad, nr-1, nr) : 0.0;
		  phexp = mexp_i ? linfit(phexp_i, rad_i, rad, nr-1, nr) : 0.0;
		  t = atan2(sint,cost);
		  if (mrot_i) vrot *= cos(mrot_i*(t-phrot));
		  if (mexp_i) vexp *= cos(mexp_i*(t-phexp));
		}
                if (aspiral>0) {                    /* spiral addition ? */
                    theta = atan2(sint,cost) - kspiral * rad * TWO_PI - theta0;
                    vspi = pow(sin(theta),nspiral) * cos(theta);
                    vrot += vrotfac * vspi;
                    vexp += vexpfac * vspi;
                }
		if (Qamp)
		  vobs = vsys + vrot;
		else
		  vobs = vsys + (vrot*cost+vexp*sint)*sin(inc);
            } else 
	        vobs = vsys;   /* near the center math goes wrong, set to vsys */

	    if (Qcube) {
	      k = gridz(vobs);
	      if (k >= 0 && k < nz) {
		ival = my_image(x,y);
		CubeValue(vptr,i,j,k) = ival;
		dprintf(1,"%d %d %d %g\n",i,j,k,vobs);
	      }
	    } else {
	      CubeValue(vptr,i,j,0) = vobs;
	      if (noise>0.0)
                CubeValue(vptr,i,j,0) += grandom(0.0,(double)noise);
	    }
        }
    }

    n=0;
    for (k=0; k<nz; k++)        /* get min and max in map */
      for (j=0; j<ny; j++)
        for (i=0; i<nx; i++) 
            if (CubeValue(vptr,i,j,k) != undef) {
                m_min = MIN(CubeValue(vptr,i,j,k),m_min);
                m_max = MAX(CubeValue(vptr,i,j,k),m_max);
            } else
                n++;
    MapMin(vptr) = m_min;
    MapMax(vptr) = m_max;
    Dx(vptr) = dx;
    Dy(vptr) = dy;
    Dz(vptr) = dz;
    Xmin(vptr) = blc[0];
    Ymin(vptr) = blc[1];
    Zmin(vptr) = blc[2];

    if (n>0) warning("%d/%d cells with no signal",n,nx*ny*nz);
    printf("Min and max in map: %g %g\n",m_min,m_max);
    write_image (outstr,vptr);      /* write out the velocity field/ data cube */
}
/*
 * RADIUS:  find out how much a point is inside or outside a
 *          projected ellipse in the deprojected plane
 *
 *  notes for improvement:  return r*r-xd*xd-yd*yd, since we're
 *                only interested in the sign
 *      make temporary arrays sin(pa), cos(pa), cos(inc)
 *      which should speed up the code quite a bit
 *  Later, since 'r' is known, the official 'd' return value can be
 *  computed
 */
 
local real radius(real r, real sinp, real cosp, real cosi, real x, real y)
{
    real xd = (-x*sinp + y*cosp);
    real yd = (-x*cosp - y*sinp)/cosi; 
#if 0
    return r-sqrt(xd*xd+yd*yd);
#else
    return r*r - xd*xd - yd*yd; /* danger: parent code needs mods */
#endif
}

/*
 * LINFIT: linear interpolation into an array with known offsets (m,n)
 *          
 */
 
local real linfit(real *f, real *x, real x0, int m, int n)
{
  if (Qlinear)
    return  ((x0-x[n])*f[m]-(x0-x[m])*f[n])/(x[m]-x[n]);
  else
    return  f[m];
}

/*
 * create a velocity field (NEMO method)
 *              0..nx-1 and 0..ny-1 
 *      start from gal plane, and project; needs to do a particle simulation
 *   
 * NOTE: we don't use this method: see snapgrid for something closer
 */

local void  vel_create_2(stream outstr)
{
    int    i, j, n, iobs, jobs, nx, ny, ir, nr, ip, np;
    real x, y, vx, vy, vrad, vtan, vobs;
    real cost, sint, cosp, sinp, cosi, sini, cosk, sink, rad, drad, phi, dphi;
    real omega, kappa, bigx, bigy, p, r, m_min, m_max, alpha = 1, q=1;
    real theta, inc;
    imageptr optr, vptr;

    warning("vel_create_2: This routine has not been upgraded");

    dprintf(0,"gridding mode: gal->sky plane\n");

    m_min = HUGE; m_max = -HUGE;
    nx = size[0];
    ny = size[1];
    cost = cos(theta); sint = sin(theta);   /* fixed for whole disk */
    cosi = cos(inc);   sini = sin(inc);

    if (!create_image(&optr, nx, ny))   /* density image */
        error("Could not create image from scratch");
    if (!create_image(&vptr, nx, ny))   /* velocity image */
        error("Could not create image from scratch");

    for (j=0; j<ny; j++)                /* init maps to zero */
        for (i=0; i<ny; i++)
            MapValue(optr,i,j) = MapValue(vptr,i,j) = undef;

    drad = MIN(cell[0],cell[1]);              /* guess for stepsize in radius */
    nr = (int) (2*(rmax-rmin)/drad);              /* use 2 times oversampling */

    for (ir=0, rad=rmin; ir<nr; ir++, rad+=drad) {     /* loop over all radii */
        if (rad==0.0) continue;                      /*  but avoid the origin */
        dphi = drad/rad;                    /* set step in angle for this rad */
        np = (int) (2*TWO_PI/dphi);               /* use 2 times oversampling */
        bigx = rad*(1-q)/(1+q);
        bigy = alpha*bigx;        /* fake for now */
        for (ip=0, phi=0; ip<np; ip++, phi+=dphi) {   /* loop over all angles */
#if 0
            cosp = cos(phi);   
            sinp = sin(phi);       /* phase along circular orbit */
            x = rad*cosp;   /* circular orbit */
            y = rad*sinp;
            vtan = 100.0;   /* circular orbit, flat rotcur */
            vrad = 0.0;
#else
            sink = sin(2*phi);
            cosk = cos(2*phi);
            p = phi-bigy*sink/rad;
            r = rad + bigx*cosk;
            cosp = cos(p);   
            sinp = sin(p);

            omega = 100/rad;        /* fake */
            kappa = 2*omega;

            vtan = omega*rad - bigy*kappa*cosk;
            vrad = -kappa*bigx*sink;

            x = r*cosp;
            y = r*sinp;
            
#endif 

            vx = vrad*cosp - vtan*sinp;             /* cartesian velocities */
            vy = vrad*sinp + vtan*cosp;

            iobs = gridx(x*cost-y*sint);            /* grid them */
            jobs = gridy(cosi*(x*sint+y*cost));
            if (iobs>=0 && jobs>=0) {               /* add if inside grid */
                vobs = -sini*(vx*sint+vy*cost);
                MapValue(optr,iobs,jobs) += 1.0;
                MapValue(vptr,iobs,jobs) += vobs;
                /* if more neigbhoring cells to be polluted, 
                 * lot of extra coding will go here 
                 */
            }
        } /* angles */
    } /* radii */

    n=0;
    for (j=0; j<ny; j++)        /* get mean velocity for cells with info */
        for (i=0; i<ny; i++) 
            if (MapValue(optr,i,j) > 0) {
                MapValue(vptr,i,j) /= MapValue(optr,i,j);
                m_min = MIN(MapValue(vptr,i,j),m_min);
                m_max = MAX(MapValue(vptr,i,j),m_max);
            } else
                n++;

    if (n>0) warning("%d cells with no signal",n);
    printf("Min and max in map: %g %g\n",m_min,m_max);
    MapMin(vptr) = m_min;
    MapMax(vptr) = m_max;
    Dx(vptr) = cell[0];
    Dy(vptr) = cell[1];
    Xmin(vptr) = blc[0];
    Ymin(vptr) = blc[1];
    write_image (outstr,vptr);      /* write out the velocity field */
}

local int gridx(real x)
{
     int i;

     i = (int) floor(0.5+center[0]+x/cell[0]) ;
     if (i<0 || i>=size[0]) 
        return -1;
     else
        return i;
}

local int gridy(real y)
{
     int i;

     i = (int) floor(0.5+center[1]+y/cell[1]) ;
     if (i<0 || i>=size[1]) 
        return -1;
     else
        return i;
}

local int gridz(real z)
{
     int i;

     i = (int) floor(0.5+center[2]+z/cell[2]) ;
#if 1
     i += size[2]/2;
#endif
     if (i<0 || i>=size[2]) 
        return -1;
     else
        return i;
}

/* ... and now for something different ....  */

real  xim = 2517;
real  yim = 2573;
real dxim = 0.025;
real dyim = 0.025;

real my_image(real x, real y)
{
  int i, j;

  i = (int) floor(0.5 + xim/2 + x/dxim);
  j = (int) floor(0.5 + yim/2 + y/dyim);

  if (i<0 || i>=Nx(iptr) || j<0 || j>=Ny(iptr)) return 0.0;

  return 255-MapValue(iptr,i,j);
}


