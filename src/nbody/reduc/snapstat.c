/*
 * STAT.C: determines various statistics from an N-body system
 *
 *       9-jan-87  P.J.Teuben  created
 *      21-jan-87  V1.1 : added options
 *      30-jul-87  V1.2 : dirty cleanup, better check it - PJT
 *       1-Nov-87  V1.3 : some cosmetic changes PJT
 *                 << database i/o is really becoming obsolete by now>>
 *      10-nov-87  V1.3b : clausius Virial added
 *       7-jun-88  V1.4: finally new filestruct
 *      20-jul-88  V1.41: --small bug in logic--
 *              BUG: Qacc/Qphi not allacted - yet wants to get Qacc - apr 89
 *	19-oct-90  V1.5 Cleaned up various - (re)allocate
 *	16-nov-90   1.5a  (re)allocation bug removed
 *       20-aug-91  correct pi/acc... allocation 
 *      20-jan-94   1.5c sqr() decl for solaris
 *      15-mar-06   1.5f use statics to hide names
 *       1-aug-06   1.5g make it listen to the times= keyword
 *      11-feb-19   1.6  add crossing time estimate
 *       8-apr-19   1.6c   fix times= bug
 *      11-apr-19   1.6d   add virial ration 2T/W
 */

/**************** INCLUDE FILES ********************************/ 

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <snapshot/snapshot.h>  

#ifndef HUGE
# define  HUGE  1e20
#endif

#ifndef TIMEFUZZ
# define TIMEFUZZ        0.0001  /* tolerance in time comparisons */
#endif

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n                   input file name ",
    "times=all\n                range of times to plot ",
    "exact=f\n                  use hackforces or exact forces ",
    "minradfrac=0.1\n           <1, for core radius det. ",
    "all=false\n                want to do all?",
    "pot=false\n                don't  all N*N calcu's ",
    "eps=0.025\n                Softening length, if needed ",
    "r_h=false\n                Want half-mass radius",
    "r_v=false\n                Want Virial radius",
    "r_c=false\n                Want core radius",
    "rms=false\n                Want rms",
    "ecutoff=0.0\n              Cutoff for bound particles",
    "verbose=t\n                verbose mode?",
    "VERSION=1.6d\n             11-apr-2019 PJT",
    NULL
};

string usage="determine various statistics of an N-body system";

string cvsid="$Id$";


/**************** SOME GLOBAL VARIABLES ************************/

local real tsnap;             /* current time of snapshot */
local int    nbody;           /* current number of used bodies */

local int mbody=0;            /* declared space for arrays */
local real *mass=NULL;        /* this will to masses array */
local real *phase=NULL;       /* will point to phases block */
local real *phi=NULL;         /* -- to potential (and eventually energies) */
local real *acc=NULL;         /* point to forces block */
local real *ax=NULL;          /* -- forces in x */
local real *ay=NULL;          /* -- forces in y */
local real *az=NULL;          /* -- forces in z */
/*-------------------------------------------------------*/
/*              Accessor macros for above vectors        */
/*-------------------------------------------------------*/
#define Pos(i)          (phase+NDIM*2*i)
#define Vel(i)          (phase+NDIM*2*i+NDIM)
#define Acc(i)          (acc+NDIM*i)
/*-------------------------------------------------------*/
/*  The following vectors also have allocated lenght 'mbody' */
local real **xp, **yp, **zp;  /* pointers to positions; for sorting */
local real **up, **vp, **wp;  /* pointers to velocities; for sorting */

local string times;                           /* input parameters */
local real minradfrac;
local real eps, sqreps;
local bool Qpot, Qr_v, Qr_c, Qr_h, Qrms, Qexact;
local bool verbose;
local real Ecutoff;
local bool need_phi, need_acc, need_rad, Qacc;

local real x1,testy1,z1,x2,y2,z2;               /* buffers used in stat analysis */
local real u1,v1,w1,u2,v2,w2;
local int  n1;
local real r2min = -1;

local real ucm, vcm, wcm;                           /* center of mass motion */
local real etot;                                    /* total energy */
local real mtot;                                    /* total mass     */
local real r_v;                                     /* Virial radius    */
local real r_c;                                     /* Core radius      */
local real r_h;                                     /* Half-mass radius */
local real r_hpx, r_hpy, r_hpz;                     /* could be done in one r_hp */
local real mass_fraction[]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,-1.0};
local real mass_radius[]  ={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

local real *epot=NULL;                              /* potential energies */
local real *rad=NULL;                               /* radii */
local int  *idr=NULL;                               /* index array for sorting */


/****************************** START OF PROGRAM **********************/

void nemo_main()
{
    stream instr;
    int    gs;
    
    times = getparam("times");
    minradfrac = getdparam("minradfrac");
    eps = getdparam("eps");
    sqreps = eps*eps;
    Qpot = getbparam("pot");
    Qr_v = getbparam("r_v");
    Qr_c = getbparam("r_c");
    Qrms = getbparam("rms");
    Qr_h = getbparam("r_h");
    Ecutoff = getdparam("ecutoff");
    verbose = getbparam("verbose");
    Qexact = getbparam("exact");
    if (getbparam("all"))
        Qpot=Qr_v=Qr_h=Qr_c=Qrms=TRUE;
    need_phi = FALSE | Qexact | Qpot;
    need_acc = FALSE | Qexact | Qpot;
    need_rad = TRUE;
    
    instr = stropen(getparam("in"), "r");

    for(;;) {
      if ( (gs = get_snap(instr)) > 0) {
	if (!streq(times,"all") && !within(tsnap,times,TIMEFUZZ)) {
	  dprintf(1,"Skipping time %g\n", tsnap);
	  continue;
	}
	analysis( nbody );          /* this scans through all particles */
	                            /* in worst case O(N*N) method */
	radii( nbody );             /* various radii of system */
	shape ();                   /* shape analysis */
      } else if (gs < 0)
	break;                      /* not snapshots found */
    } 
}

/*
 *  get_snap:  get next valid snapshot from data stream
 *      returns:   -1:  not a Snapshot, take it as EOF
 *                  0:  no Particles here, but not EOF yet
 *		    1:  yes, Particles here
 */

       /* NOTE: coordinate system is assumed to be cartesian */

get_snap(stream instr) /* returns:   -1: not a snapshot   0: no ParticlesTag */
{
    int i, j;
    real *p;

    dprintf(1,"get_snap\n");
    
    get_history(instr);         /* just to be safe */

    if (!get_tag_ok(instr, SnapShotTag)) {    /* must be a snapshot */
      dprintf(1,"not a SnapShot: \n");
      return -1;
    }

    get_set(instr, SnapShotTag);
      get_set(instr, ParametersTag);
        get_data(instr, NobjTag, IntType, &nbody, 0);
        if (get_tag_ok(instr, TimeTag))
           get_data_coerced(instr, TimeTag, RealType, &tsnap, 0);
        else
           tsnap=0.0;
        dprintf (2,"SnapShot : Time=%f\n",tsnap);
      get_tes(instr, ParametersTag);

      if (!get_tag_ok(instr,ParticlesTag)) {    /* if no particles, return */
         get_tes(instr,SnapShotTag);            /* with OK status anyhow   */
         return(0);                             /* data is probably diags  */
      }
      get_set(instr, ParticlesTag);             /* OK, we have particles   */
         snap_alloc();                         /* see if more space needed */
         /* BUG in obtaining Masses from the snapshot                      */
                        /* this mechanism is not fool proof: if snapshots  */
                        /* have different size, and only the first one has */
                        /* masses, subsequents runs may have faul masses   */
         if (get_tag_ok(instr,MassTag)) {
            dprintf(2,"Reading %d masses\n",nbody);
            get_data_coerced(instr, MassTag, RealType, mass, nbody, 0);
         }

         if (get_tag_ok(instr,PhaseSpaceTag)) {
            dprintf (2,"Reading %d phasespace \n",nbody);
            get_data_coerced(instr,PhaseSpaceTag,RealType,phase,nbody,2,NDIM,0);
         } else {
            warning("No phasespace in snapshot: time=%f",tsnap);            
            get_tes(instr,ParticlesTag);
            get_tes(instr,SnapShotTag);
            mbody = MAX(mbody,nbody);    /* still update the alloc counter */
            return(0);
         }

         for (i=0; i<nbody; i++) {          /* fill in radii */
            rad[i] = 0.0;
            for (j=0; j<NDIM; j++)
                rad[i] += sqr(Pos(i)[j]);
            rad[i] = sqrt(rad[i]);
         }
         if (Qexact)
            exact(); /* timeconsuming - calculate exact forces/potential */
         else {
            if (need_phi) {
                if (get_tag_ok(instr,PotentialTag)) {
                    dprintf (2,"Reading %d potentials\n",nbody);
                    get_data(instr, PotentialTag, RealType, phi, nbody, 0);
                } else
                    error("Need potentials in this snapshot or use exact=t");
            }
            if (need_acc) {
               if (!get_tag_ok(instr,AccelerationTag)) {
                  error("Need forces in this snapshot");
                  dprintf (2,"Reading %d accelarations\n",nbody);
                  get_data(instr, AccelerationTag, RealType, acc, nbody, NDIM, 0);
                  for (i=0, p=acc; i<nbody; i++) {        /* kludge */
                    ax[i] = *p++;
                    ay[i] = *p++;
                    az[i] = *p++;
                  }
	       }
            }
         }
      get_tes(instr, ParticlesTag);
    get_tes(instr,SnapShotTag);

    if (NDIM!=3) error("Program only works with 3D data");

    for (i=0, p=phase; i<nbody; i++) {   /* (re)set pointers to phases */
        xp[i] = p++;    yp[i] = p++;   zp[i] = p++;
        up[i] = p++;    vp[i] = p++;   wp[i] = p++;
    }
    mbody = MAX(mbody,nbody);           /* reset allocate space counter */
    return 1;
}

snap_alloc()
{
    mass = (real *) reallocate(mass, nbody*sizeof(real));
    phase = (real *)reallocate(phase,nbody*2*NDIM*sizeof(real));
    xp = (real **)  reallocate(xp,   nbody*sizeof(real *));
    yp = (real **)  reallocate(yp,   nbody*sizeof(real *));
    zp = (real **)  reallocate(zp,   nbody*sizeof(real *));
    up = (real **)  reallocate(up,   nbody*sizeof(real *));
    vp = (real **)  reallocate(vp,   nbody*sizeof(real *));
    wp = (real **)  reallocate(wp,   nbody*sizeof(real *));
    rad = (real *)  reallocate(rad,  nbody*sizeof(real));
    idr = (int *)   reallocate(idr,  nbody*sizeof(int));
    if (need_phi) {
        phi  = (real *) reallocate(phi, nbody*sizeof(real));
        epot = (real *) reallocate(epot,nbody*sizeof(real));
    }
    if (need_acc) {
        acc = (real *)reallocate(acc,nbody*NDIM*sizeof(real));
        ax = (real *) reallocate(ax, nbody*sizeof(real));
        ay = (real *) reallocate(ay, nbody*sizeof(real));
        az = (real *) reallocate(az, nbody*sizeof(real));
    }
}

exact()                         /* exact potential and forces */
{
    int i,j,k;
    real rij;

    dprintf (2,"Doing an exact potential calculation\n");    
    for (i=0; i<nbody; i++)
        phi[i] = 0;
    for (i=1; i<nbody; i++) {
        for (j=0; j<i; j++) {
            rij = 0.0;
            for (k=0; k<NDIM; k++)
                rij += sqr(Pos(i)[k] - Pos(j)[k]);
            rij = 1.0/sqrt(sqreps + rij);
            phi[i] -= mass[j]*rij;              /* G==1 */
            phi[j] -= mass[i]*rij;              /* G==1 */
        }
    }
    rij=0;
    for (i=0; i<nbody; i++)
        rij +=  phi[i];
    dprintf(2,"Total of body potentials = %f\n",rij);
}



/*
 *  ANALYSIS:
 *    -  mean position & velocities as well as their dispersion
 *    -  total energy = kinetic + potential energy
 */

analysis(int nbody)
{
        int i,j,i1,i2;
        real x,y,z,u,v,w;
        real inv_rad, pmass, tmp, r2;
        real xdir, ydir, zdir, artmp, axtmp, aytmp, aztmp;

        ini_analysis();
        if (nbody>200 && Qpot && verbose)
                dprintf (1,"Be patient...this operation takes a while\n");
                
        for (i=0; i<nbody; i++) {
                x = *xp[i]; y = *yp[i]; z = *zp[i];     /* positions  */
                u = *up[i]; v = *vp[i]; w = *wp[i];     /* velocities */
                pmass = mass[i];                              /* mass */
                add_analysis (i,pmass,x,y,z,u,v,w);   /* add to analysis */

                if (Qpot || Qr_v) 
                   for (j=i+1; j<nbody; j++) {
                        xdir = x - *xp[j];
                        ydir = y - *yp[j];
                        zdir = z - *zp[j];
                        r2 = sqr(xdir) + sqr(ydir) + sqr(zdir);
                        dprintf (2,"%d %d %f\n",i,j,sqrt(r2));
                        if (r2<r2min) {
			  i1=i;
			  i2=j;
			  r2min=r2;
			  dprintf(1,"rmin=%g for (%d,%d) mass (%g,%g)\n",sqrt(r2min),i1,i2,mass[i1],mass[i2]);
			}
                        if (!Qacc) {
                          artmp = mass[j]/((r2+sqreps)*sqrt(r2));  /* tmp force */
                          axtmp = xdir*artmp;
                          aytmp = ydir*artmp;
                          aztmp = zdir*artmp;
                          ax[i] -= axtmp;               ax[j] += axtmp;
                          ay[i] -= aytmp;               ay[j] += aytmp;
                          az[i] -= aztmp;               az[j] += aztmp;
                        }
                        if (Qpot) {
                           tmp = mass[j] * pmass / sqrt(r2 + sqreps);
                           epot[i] -= tmp;
                           epot[j] -= tmp;
                        }
                        if (Qr_v)
                           r_v += 1.0/sqrt(r2);
                   }
        }
        report_analysis();
}

/*
 *  Some routines to calculate mean positions and velocities & 
 *  accompanying utilities
 */
 
ini_analysis()
{                               /*   initialize some common stuff   */
        int  i;
        
        n1=0;                                           /* total part. */
        x1=x2=testy1=y2=z1=z2=u1=u2=v1=v2=w1=w2=0.0;        /* pos & vel */
        if (Qpot)
            for (i=0; i<nbody; i++) {
                epot[i]=0.0;                            /* energy */
                ax[i]=0.0;                        /* forces */
                ay[i]=0.0;                        /* note: 3b1 compiler */
                az[i]=0.0;                        /* error here prevented */
            }
        r_v=0.0;                                        /* virial radius */
        mtot=0.0;                                       /* mass */
        r2min = HUGE;                           /* min interparticle distance */
}

add_analysis(ipart,pmass,x,y,z,u,v,w)
real pmass,x,y,z,u,v,w;
int    ipart;
{
        real r2;

        n1 += 1;                                        /* statistics */
        x1 += x;        testy1 += y;        z1 += z;        /* on whole   */
        x2 += x * x;    y2 += y * y;    z2 += z * z;    /* system     */
        u1 += u;        v1 += v;        w1 += w;
        u2 += u * u;    v2 += v * v;    w2 += w * w;
        
        mtot += pmass;                                  /* tot mass  */
}       

report_analysis()
{
        real r,v;
        real xm,ym,zm,xs,ys,zs;
        real um,vm,wm,us,vs,ws;
        real rmsvel, t_cr;
        real ecc_x, ecc_y, ecc_z;
        real ekintot, epottot, ecomtot, e, emin, emax, cv, ecm, ekin;
        int    i, nplus;
                                        
        if (n1<2) {
                printf ("REPORT_ANALYSIS: n1=%d\n",n1);
                return(0);
        }
        if (verbose) {
            printf ("Nobj= %d\n\nTime of snapshot= %f\n",nbody,tsnap);

           printf ("\n\nTotal mass = %f\n", mtot);
           printf ("Mean position & velocities of %d particles: \n",n1);
        }
        ms (n1,x1,x2,&xm,&xs);
        ms (n1,testy1,y2,&ym,&ys);
        ms (n1,z1,z2,&zm,&zs);

        ms (n1,u1,u2,&um,&us);  ucm=um;
        ms (n1,v1,v2,&vm,&vs);  vcm=vm;
        ms (n1,w1,w2,&wm,&ws);  wcm=wm;
        if (verbose) {  

        
        printf ("pos:  %f +/- %f    %f +/- %f    %f +/- %f\n",
                       xm  ,  xs,   ym,    ys,   zm,    zs);
        printf ("vel:  %f +/- %f    %f +/- %f    %f +/- %f\n\n",
                       um  ,  us,   vm,    vs,   wm,    ws);
        }
        if (verbose) { 
           printf ("Smallest interparticle distance = %f\n",sqrt(r2min));
        } 
        rmsvel = sqrt(us*us+vs*vs+ws*ws);
        if (Qrms) {
           if (verbose)
                printf ("RMS velocity = ");
           printf ("%f\n",rmsvel);
        }

        if (!Qpot)                      /* if no potentials available */
           return(0);                   /* return now */
                                        /* otherwise do more here */
        ekintot = epottot = ecomtot = 0.0;      /* reset totals */
        ecm = 0.5*(sqr(ucm)+sqr(vcm)+sqr(wcm)); /* kin energy of COM */
        nplus=0;                        /* # of stars with positive energy */
        for (i=0; i<nbody; i++) {       /* determine min and max of energies */
                ekin=0.5*mass[i] * 
                   ( sqr(*up[i]-ucm) + sqr(*vp[i]-vcm) + sqr(*wp[i]-wcm) );
                e=ekin + epot[i];
                if (i==0) emin=emax=e;
                emin = MIN(emin,e);
                emax = MAX(emax,e);
                if (e>0.0) nplus++;     /* count 'unbound' particles */
                ekintot += ekin;        /* accumulate total energies */
                epottot += epot[i]; 
                ecomtot += mass[i]*ecm;
        }
        epottot *= 0.5;         /* correct  (ij) with (ji) for double count */
        printf ("E (int energy) = T (kinetic) + U (potential)  E(kin com): \n");
        printf ("%f   =   %f  +  %f  +  %f \n",
                ekintot+epottot, ekintot, epottot,ecomtot);
        printf ("%d stars with positive energy in COM frame; emin=%f emax=%f\n",
                nplus, emin, emax);

        cv = 0.0;               /* Clausius Virial */
        for (i=0;  i<nbody; i++) 
                cv += mass[i] * (*xp[i]*ax[i] + *yp[i]*ay[i] + *zp[i]*az[i]);
        printf ("Clausius energy = %f\n",cv);
        printf ("Virial = %f (Clausius => %f)    2T/W = %f  (Clausius => %f)\n",
                epottot+2*ekintot,cv+2*ekintot,-2*ekintot/epottot, -2*ekintot/cv);

	t_cr = pow(mtot,2.5) / pow(2*ABS(ekintot+epottot),1.5);
	printf("Crossing time = %f\n", t_cr);
	
}


ms (sumn1, sumx1, sumx2, mean, sigma)
int sumn1;                      /*  input:  number of values used to sum.. */
real  sumx1,sumx2;            /*  input:  sum of X and X^2               */
real  *mean,*sigma;           /*  output: mean X and dispersion in X     */
{
        *mean  = sumx1/(real)sumn1;                           
        *sigma = sqrt(sumx2/(real)sumn1 - sqr(*mean));
}


/*
 *  Following routines deal with calculation of typical radii,
 *  half_mass_radii, either intrinsic or projected, virial radius,
 *  etc.
 *
 */
 
radii(int n)
{
/*      real halfmass_radius();                               */
        real radius, dr, drmin, r2, sum, half_sur_den_0;
        int    ntot, i, low, mid, high;
        
        sortptr (rad, idr, n);          /* sort radii in index array */
        if (Qr_h) {
           mass_radii (mass_fraction, mass_radius);
           if (verbose)
                printf ("Time ");
           printf ("%f ",tsnap);
           i=0;
           if (verbose)
                printf ("Radii at massfractions 0.1,0.2,...0.9\n");
           while (mass_fraction[i]>0)
                printf ("%8.5f ",mass_radius[i++]);
           printf ("\n");
        }



/*   
 *    Generalized Surface density test  
 */
 
    if (Qr_c) {
        
        /*  first compute smallest interparticle distance  */
 
        drmin = rad[idr[n-1]];                  /* largest distance */
        for (i=1; i<n; i++) {           /* because arrays were sorted */
                dr = rad[idr[i]] - rad[idr[i-1]];
                if (dr<drmin)
                        drmin=dr;
        }
        drmin *= minradfrac;            /* and take a fraction of that */
	/*      printf ("SD: drmin*minradfrac=%20.10e\n",drmin);                */


        /*  find central surface density */
        sum=0.0;
        for (i=0; i<n; i++)
                sum += mass[idr[i]] / ( sqr(rad[idr[i]]) );
        half_sur_den_0 = 0.5 * sum;
        printf ("SD: central surface brightness = %f\n",sum);   

	if (n>2) {
          /*  then subdivide interval until surface density half the central   */
	  low = 0; high = n-1;
	  while ((high-low)>1) {
                mid = (high+low)/2;
                radius=rad[idr[mid]] + drmin;
                /*              printf ("RC: Radius #%d = %f ",mid,radius);     */
                sum=0.0;
                for (i=mid+1; i<n; i++)
                   sum += mass[idr[i]] / ( rad[idr[i]] * sqrt (
                          (rad[idr[i]] - radius)*(rad[idr[i]] + radius) ) );
		/*              printf (" sur_den = %f\n",sum);                 */
                if (sum>half_sur_den_0)
                        low = mid;
                else
                        high = mid;
	  }       
	  r_c = rad[idr[mid]];                    /* core radius */
	  printf ("\nr_c = %f\n",r_c);
	}

    }  /* end Qr_c  */
                

#if 0
        project_radius (xp, yp, n);
        r_hpz = sqrt( halfmass_radius() );
        printf ("XY-projected half_mass radius = %f\n", r_hpz);

        project_radius (xp, zp, n);
        r_hpy = sqrt( halfmass_radius() );
        printf ("XZ-projected half_mass radius = %f\n", r_hpy);

        project_radius (yp, zp, n);
        r_hpx = sqrt( halfmass_radius() );
        printf ("YZ-projected half_mass radius = %f\n", r_hpx);
#endif

        if (Qr_v) {             
           r_v = 0.5 * n * (n-1) / r_v;
           printf ("Virial radius r_v = %f\n",r_v);
        }
}

/*
 *   assumes that idx[] and rad[] have been properly filles up to 'n'
 *   such that rad[idx[i]] 1=0..n-1 are sorted in ascending order
 *
 *   returns then the halfmass radius (either intrinsic or projected,
 *   depending how rad[] was filled in (see: project_radius)
 */
 
mass_radii (real *mf, real *mr)
{
        real cmass, fmass;
        int i, k;
        
        k=0;                    /* radii */
        cmass=0.0;              /* cumulative */
        fmass = mf[k] * mtot;   /* first cumulative mass */
        i=0;                    /* particle counter */
        do {
              cmass += mass[idr[i++]];          /* add up sorted mass */
              if (cmass>fmass) {
                  mr[k++] = rad[idr[i-1]];
                  fmass = mf[k] * mtot;         /* next cumul. mass */
              }
        }
        while (i<nbody && fmass>0.0);
}

/*
 *    combines an *[x] and *[y] array to form a projected squared 
 *    radius in rad[]
 */
  
project_radius(x,y,n)
real *x[], *y[];
int n;
{
        int i;
        double sqr();
        
        for (i=0; i<n; i++)
                rad[i] =  sqr(*x[i]) + sqr(*y[i]);

        sortptr(rad,idr,n);

}

/*
 *  SHAPE: shape analysis on particles with energy less than a certain 
 *         (upper limit) cutoff energy (getdparam("ecutoff"))
 */

shape()
{
        int i;
        real x,y,z;
        real xm,ym,zm,xs,ys,zs;
        real e, ecm, ekin;
        
        if (!Qpot)                      /* we do need energies here */
           return(0);
          
        n1=0;                                           /* total part. */
        x1=x2=testy1=y2=z1=z2=0.0;                          /* pos & vel */
        ecm = 0.5*(sqr(ucm)+sqr(vcm)+sqr(wcm));
                /* now a shape analysis on particles with e<Ecutoff */

        for (i=0; i<nbody; i++) {
                ekin=0.5*mass[i] * 
                   ( sqr(*up[i]-ucm) + sqr(*vp[i]-vcm) + sqr(*wp[i]-wcm) );
                e = ekin + epot[i];
                if (e<Ecutoff) {
                        x = *xp[i]; y = *yp[i]; z = *zp[i];
                        n1 += 1;
                        x1 += x;        testy1 += y;        z1 += z;
                        x2 += x * x;    y2 += y * y;    z2 += z * z;
                }
        }
        if (n1>0) {     /* form mean and dispersion in all three coordinates */
            ms (n1,x1,x2,&xm,&xs);
            ms (n1,testy1,y2,&ym,&ys);
            ms (n1,z1,z2,&zm,&zs);

            printf ("Shape analysis for E<%f: (%d particles)\n",Ecutoff,n1);
            printf (" %f +/- %f    %f +/- %f    %f +/- %f \n",
                xm,xs, ym,ys, zm,zs);
        } else
            printf ("No particles with energy below cutoff %f\n",Ecutoff);
}
