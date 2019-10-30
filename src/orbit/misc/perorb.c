/*
 *  PERORB:    find periodic orbits from a near guessed starting one
 *
 *      The periodic orbit (PO) must be planar and symmetric to at least 
 *      one principle plane.
 *
 *      It iterates by converging in the Surface Of Section (SOS) plane
 *      X-VX or Y-VY. (VX and VY are in rotating frame of reference)
 *
 *      Technique has been used in various forms by various researchers, 
 *	and is described in its most original form in Henon (1965,
 *	Ann.Astr. xxx, yyy.
 * 	See also:
 *	  - El-Sabaa & Sherief (1990) Ap&SS 167,305.
 *
 *      Cross-correllated SOS's may be implemented later.
 *      no real plotting yet, just a table
 *
 *	   mar-82         original created as POSOS on Cyber 7600 (Univ.Gron)
 *	                  (Sheltran)
 *         nov-82         apparently improved...
 *      22-may-90   V1.0  Kludged for NEMO for Pittsburgh 1990 Workshop    PJT
 *      27-may-90   V1.0a enhanced with SOS plot options                   PJT
 *      24-may-91   V1.1  rotating frame of reference kludge --            PJT
 *      13-jun-91   V1.2  SOS option                                       PJT
 *			  plus lots of debugging + better Leapfrog method  PJT
 *      19-jun-91   V1.3  optional input file initial conditions           PJT
 *	24-may-92	b add_history -> app_history; added <potential.h>  PJT
 *	22-jul-93   V1.4a fix a_potential declaration bug
 *	23-mar-95       b get_pattern
 *      28-mar-95   V1.5  added mode= to select integration                PJT
 *      22-apr-96   V1.5b allowed orbits which retrace themselves          pjt
 *                        e.g. R-z periodic orbits
 *			-- oops, not neeeded, just use period=1	--
 *      10-apr-97       c no more 'this is beta' message                   pjt
 *      10-jan-03       d SGN -> SIGN                                      pjt
 *       1-feb-03       e report energy conservation, use new IOM errors   pjt
 *      12-apr-04       f SIGN -> SNG (since NumRec uses  SIGN)            pjt
 *      29-oct-2019 V1.7  allow phase= to have 2 or 3 values               pjt
 *
 *  TODO:   check why rk2 and rk4 both seem to converge at the same rate
 *          perhaps SOS is linearly interpolated, higher order rk4 defeats this?
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>

#include <filestruct.h>
#include <potential.h>
#include <orbit.h>

string defv[] = {
    "in=\n                     Optional input orbit file for initial P.O.",
    "out=\n                    Optional output orbit file",
    "freqout=1\n               Frequency of output in steps",
    "maxsteps=5000\n           Maximum number of steps for one period",
    "dt=0.01\n                 Timestep",
    "phase=0,1,0,0.6,0,0\n     Initial phase-space coordinates[,energy]",
    "step=0.1,0.01\n           Step in {pos|energy}, perturbing pos_step",
    "dir=y\n                   Direction to step in (x or y)",
    "norbit=1\n                Number of orbits to do",
    "accuracy=0.0001\n         Relative accuracy in p.o.",
    "maxiter=50\n              Maximum iterations allowed to find p.o.",
    "ncross=1\n                Number of SOS crossings per iteration",
    "period=2\n                Number of SOS points per period",
    "potname=\n                Name of external potential(5)",
    "potpars=\n                Parameters of potential(5)",
    "potfile=\n                Extra data-file for potential(5)",
    "tab=\n                    Optional table with x,vy,y,vx,nsteps,T,E",
    "mode=rk4\n                integration method (euler,leapfrog,rk2,rk4)",
    "last=f\n                  Use the last orbit in the in= orbit file?",
    "headline=\n               Random verbiage for output file",
    "VERSION=1.7\n             29-oct-2019 PJT",
    NULL,
};

string usage="search for periodic orbits in a potential (SOS method)";

string  infile, outfile;                /* file names */
stream  instr, outstr;                  /* file streams */
    
orbitptr o1=NULL, o2=NULL;              /* pointer to orbits */

real   dt, dt2;                         /* stepping */
proc   pot;
a_potential p;

double phase[2*NDIM+1];                 /* initial phase coord's + optional energy */
double phasestep, smallstep;
double *perorb=NULL;                    /* point to final answers*/
double *xsos, *ysos, *tsos, *esos;      /* store all SOS coord's here */
int norb, iorb;                         /* keep track which orbit */
int icross, ncross;                     /* SOS crossings */
int maxiter;                            /* max. iterations allowed */
int dirint, maxsteps, period, freqout, sign;
double eps_min, eps_max;
double omega, omegasq, omegato;         /* pattern speed, and related ^2,*2 */
FILE *tab;
bool  Qfix;                             /* flag if to step in Energy instead */
bool  Qlast;                            /* last orbit in in= ?? */
iproc cycle;                            /* integrator */

string integration_modes = "euler leapfrog rk2 rk4";
int   cycle_euler(orbitptr), 
      cycle_leapfrog(orbitptr),
      cycle_rk2(orbitptr), 
      cycle_rk4(orbitptr);

/*----------------------------------------------------------------------------*/

nemo_main()
{
  int i, n, maxout, nold;

    if (NDIM!=3) error("Program can only run with NDIM=3");

    setparams();
    if (*outfile && ncross==1)
        outstr = stropen (outfile,"w");
    else {
        warning("DRY RUN: No output orbit file created");
        outstr = NULL;
    }

    if (outstr) {
	app_history(getparam("headline"));
        put_history(outstr);
    }
    for (iorb=0; iorb<norb; iorb++) {
        if (iorb>0) Qfix = FALSE;       /* only use Energy once */
        prepare();                   /* initialize or get a good estimate */
        n = iterate();      /* iterate and return number of steps in orbit */
        if (n==0) {
            warning("Bad convergence or incomplete orbit %d - skipped",iorb);
            continue;
        }
        maxout = Nsteps(o1);               /* save how much space we have */
        nold   = Nsteps(o1);               /* save how much space we have */
        maxout = MAXsteps(o1);               /* save how much space we have */
        Nsteps(o1) = n;                           /* set these for output */
        if (outstr) my_write_orbit(outstr, o1, maxout, freqout, period);
        Nsteps(o1) = nold;
        if (ncross>1) {
            if (dirint==1)
                dprintf(0,"Time     Ycross    VYcross   Energy\n");
            else if (dirint==0)
                dprintf(0,"Time     Xcross    VXcross   Energy\n");
            else
                dprintf(0,"Tsos,Xsos,Ysos,Esos: bogus SOS table\n");
            for (i=0; i<ncross; i++) {
                printf("%g %g %g %g\n",tsos[i], xsos[i], ysos[i], esos[i]);
            }
        } 
    }
    if (outstr) strclose(outstr);
}

setparams()
{
    char *cp;
    int   n, ndim, s, imode;
    double acc[3], e_pot, time, vel;
        
    maxsteps = getiparam("maxsteps");
    if (allocate_orbit (&o1,NDIM,maxsteps)==0)
        error ("Error in allocating orbit 1");
    if (allocate_orbit (&o2,NDIM,maxsteps)==0)
        error ("Error in allocating orbit 2");

    dt = getdparam("dt");
    dt2 = 0.5*dt;

    maxiter = getiparam("maxiter");
    Qlast = getbparam("last");

    eps_min = getdparam("accuracy");
    eps_max = 1/eps_min;    /* just a guess */

    cp = getparam("dir");
    switch (*cp) {
         case 'x': dirint =0; break;
         case 'y': dirint =1; break;
         default:  error("Invalid stepper direction %s (must be x or y)\n",cp);
    }

    n = nemoinpd(getparam("phase"),phase,2*NDIM+1);
    if (n==2 || n==3) {
      if (n==3) phase[6] = phase[2];
      if (dirint==0) {         /* (x0,v0[,E]) */
	phase[4] = phase[1];
	phase[1] = phase[2] = phase[3] = phase[5] = 0.0;
      } else if (dirint==1) {  /* (y0,u0[,E]) */
	phase[3] = phase[1];
	phase[1] = phase[0];
	phase[0] = phase[2] = phase[4] = phase[5] = 0.0;
      } else error("illegal dirint");
      if (n==2)
	Qfix = FALSE;
      else
	Qfix = TRUE;
    } else if (n!=2*NDIM) {        /* if not just pos&vel given: */
        if (n==2*NDIM+1) {  /* accept the last one as a fixed energy */
            dprintf(0,"Using %g as (fixed) energy for first orbit\n",phase[2*NDIM]);
            Qfix = TRUE;
        } else              /* otherwise some error */
            error("Wrong number (%d) of phase space coordinates",n);
    } else
        Qfix = FALSE;

    if (hasvalue("in")) {	    /* initial conditions from an old orbit */
    	infile = getparam("in");
        instr = stropen(infile,"r");
        p.name = NULL;
        p.pars = NULL;
        p.file = NULL;
        while (read_orbit(instr,&o1)) {
            phase[0] = Xorb(o1,0);
            phase[1] = Yorb(o1,0);
            phase[2] = Zorb(o1,0);
            phase[3] = Uorb(o1,0);
            phase[4] = Vorb(o1,0);
            phase[5] = Worb(o1,0);
            phase[6] = I1(o1);
            p.name = PotName(o1);
            p.pars = PotPars(o1);
            p.file = PotFile(o1);
	    if (!Qlast) break;
	}
	dprintf(0,"[Read initial conditions from %s, maxsteps=%d\n]",infile,MAXsteps(o1));
        strclose(instr);
    }
    if(hasvalue("potname"))p.name = getparam("potname");
    if(hasvalue("potpars"))p.pars = getparam("potpars");
    if(hasvalue("potfile"))p.file = getparam("potfile");
    if(p.name==NULL || *p.name==0) p.name="plummer";
    pot=get_potential(p.name,p.pars,p.file);
    if (pot==NULL) error("Requested potential \"%s\" could not be loaded",
		getparam("potname"));
    omega = get_pattern();
    dprintf(1,"Pattern speed: %g\n",omega);
    omegasq = omega * omega;
    omegato = 2.0*omega;

    n = nemoinpd(getparam("step"),acc,2);   /* get steps; use acc[] as temp */
    if (n<1 || n>2) error("Parsing step=%s",getparam("step"));
    phasestep = acc[0];
    if (n==2)
        smallstep = acc[1];                 /* perturbing step for iteration */
    else
        smallstep = 0.1*phasestep;          /* take default of 10% of step */

    freqout = getiparam("freqout");

    ncross = getiparam("ncross");
    if (ncross > 0) {
      xsos = (double *) allocate(ncross * sizeof(double));
      ysos = (double *) allocate(ncross * sizeof(double));
      tsos = (double *) allocate(ncross * sizeof(double));
      esos = (double *) allocate(ncross * sizeof(double));
    } else
        error("Need ncross>0 (%d)\n",ncross);
    if (ncross>1)
        dprintf(0,"PERORB is now in interactive mode - no output\n");

    period = getiparam("period");
    if (period < 1 || period > 2) error("This version needs period to be 1 or 2");
    sign = (period==1 ? 1 : -1 );

    if (Qfix) {       /* get the true launch velocity if Energy was supplied */
        ndim = NDIM;
        time = 0.0;
        (*pot)(&ndim, phase, acc, &e_pot, &time); /* only to get e_pot */
        if (omega!=0)  /* NDIM=3 */
            e_pot -= 0.5*omegasq*( sqr(phase[0]) + sqr(phase[1]) );
        s = SGN(phase[NDIM+1-dirint]);         /* sign of velocity */
        vel = 2*(phase[2*NDIM]-e_pot);    /* new velocity squared */
        if (vel<0)                  /* check if vel not negative */
            error("V^2=%g < 0; try different position or energy",vel);     
        vel = sqrt(vel);
        if (s<0)
            vel = -vel;
        else if (s==0)
            warning("Sign of starting velocity assumed positive!!");
        phase[NDIM+1-dirint] = vel;                 /* and set it */
        dprintf(0,"Calculated launch velocity: %g\n",vel);
    }
    /* Now check phase space coordinates */
    if (phase[dirint]==0.0)   
        warning("launch position zero");
    if (phase[1-dirint]!=0.0) 
        warning("non-launch position non-zero",phase[1-dirint]);
    if (phase[NDIM+1-dirint]==0.0)   
        warning("launch velocity zero");
    if (phase[NDIM+dirint]!=0.0) 
        warning("non-launch velocity non-zero",phase[NDIM+1-dirint]);


    norb = getiparam("norbit");
    perorb = (double *) allocate(norb*NDIM*6);    /* store phase  */
        
    outfile = getparam("out");

    if (hasvalue("tab")) {
        cp = getparam("tab");
        dprintf(0,"Table with results in %s\n",cp);
        tab = stropen(cp,"w");
    } else
        tab = stdout;        

    match(getparam("mode"),integration_modes,&imode);
    if (imode==0x01)
        cycle = cycle_euler;
    else if (imode==0x02)
        cycle = cycle_leapfrog;
    else if (imode==0x04)
        cycle = cycle_rk2;
    else if (imode==0x08)
        cycle = cycle_rk4;
    else
      error("imode=0x%x; Illegal integration mode=%s, try one of: %s",
	    imode,getparam("mode"),integration_modes);
}


my_write_orbit(stream outstr, orbitptr o, int maxout, int freqout, int period)
{
    int i,j, nout;

    if (freqout > 1) {                    /* if only part of it to write .. */
        for (i=freqout, j=1; i<Nsteps(o); i+=freqout, j++)         /* shift */
            move_o(o,j,i);
        nout = j;        
        if ((Nsteps(o)-1) % freqout)      /* also force last one if odd end */
            move_o(o,nout++,Nsteps(o)-1);            
    } else
        nout = Nsteps(o);                       /* otherwise all is written */

    if (period==2) {        /* symmetrize if only half the orbit was stored */
        if (2*nout-1 <= maxout) {
          for (i=nout; i<2*nout; i++) {     /* symmetrize orbit in XY coord's */
            Posorb(o,i,1-dirint) = -Posorb(o,2*nout-i-2,1-dirint);
            Posorb(o,i, dirint ) =  Posorb(o,2*nout-i-2, dirint );
            Posorb(o,i,2)        =  Posorb(o,2*nout-i-2,2);
            Velorb(o,i,1-dirint) =  Velorb(o,2*nout-i-2,1-dirint);
            Velorb(o,i, dirint ) = -Velorb(o,2*nout-i-2, dirint );
            Velorb(o,i,2)        =  Velorb(o,2*nout-i-2,2);
          }
          nout = 2*nout-1;
       } else
	 warning("Not enough space to save symmetric part of orbit: nout=%d maxsteps=%d",
		 nout,maxout);
    }
    Nsteps(o) = nout;   /* make sure only all these are written */
    PotName(o) = p.name;    /* make sure any new potential is known */
    PotPars(o) = p.pars;
    PotFile(o) = p.file;
    write_orbit(outstr,o);
}

move_o (orbitptr o, int j, int i)        /* move all orbit elements from i to j */
{
    Torb(o,j) = Torb(o,i);
    Xorb(o,j) = Xorb(o,i);
    Yorb(o,j) = Yorb(o,i);
    Zorb(o,j) = Zorb(o,i);
    Uorb(o,j) = Uorb(o,i);
    Vorb(o,j) = Vorb(o,i);
    Worb(o,j) = Worb(o,i);
}

prepare()
{
    real acc[3], pos[3], e_pot, time, vnew;
    int ndim, dir0;

    if (iorb==0) {              /* first orbit: total initialize */
        Masso(o1) = 1.0;
        Xorb(o1,0) = phase[0];
        Yorb(o1,0) = phase[1];
        Zorb(o1,0) = phase[2];
        Uorb(o1,0) = phase[3];
        Vorb(o1,0) = phase[4];
        Worb(o1,0) = phase[5];
    } else {	         /* for third orbit a better extrapol. could be done */
        if (Qfix) {
            phase[2*NDIM] += phasestep;     /* tweak the energy */
            dir0 = SGN(Velorb(o1,0,1-dirint));
            time = Torb(o1,0);
            pos[0] = Xorb(o1,0);
            pos[1] = Yorb(o1,0);
            pos[2] = Zorb(o1,0);
            ndim=3;
            (*pot)(&ndim,pos,acc,&e_pot,&time);
            e_pot -= 0.5*omegasq*(sqr(pos[0])+sqr(pos[1]));
            vnew = 2*(phase[2*NDIM]-e_pot);
            if (vnew<0.0) {
                vnew = 0.0;
                warning("Resetting vnew=0 at xnew=%f  e_tot=%g\n",
                         Posorb(o1,0,dirint), e_pot);
            }
            Velorb(o1,0,1-dirint) = sqrt(vnew) * dir0;
        } else {            /* this has to be improved */
            Posorb(o1,0,dirint) += phasestep;
        }
    }
}

/*
 *  ITERATE:  iterates until the orbit matches itself reasonably well.
 *            on returning to its point of outset (or T/2 away from it)
 *              If no match is found, returns 0
 *              else it returns the number of steps for that period.
 *      If ncross>1 no iteration, but plain vanilla SOS coords are obtained
 */
iterate() 
{
  int i, iter, l1, l2, h, ndim, dir0;
  double eps, x1_s, x2_s, u1_s, u2_s, x1_e, x2_e, u1_e, u2_e, de;
  double a_s, b_s, a_e, b_e, xnew, vnew, time, e_pot, e_tot, lz_mean;
  double pos[NDIM], acc[NDIM];
  static bool first_out = TRUE;

  time = 0.0;   /* keep it fixed */
  ndim = Ndim(o1);

  dir0 = SGN(Velorb(o1,0,1-dirint));    /* initial direction of orbit */

  l1 = (*cycle)(o1);                  /* advance primary orbit (by T/period) */
  dprintf(1," cycle1 @ %d\n",l1);
  dprintf(1,"POS,VEL= %f %f %f %f\n         %f %f %f %f\n",
        Xorb(o1,0), Yorb(o1,0), Uorb(o1,0), Vorb(o1,0),
        Xorb(o1,l1-1), Yorb(o1,l1-1), Uorb(o1,l1-1), Vorb(o1,l1-1));

  pos[0] = Xorb(o1,0);
  pos[1] = Yorb(o1,0);
  pos[2] = Zorb(o1,0);
  (*pot)(&ndim, pos, acc, &e_pot, &time); /* only to get e_tot */
  e_pot -= 0.5*omegasq*(sqr(pos[0]) + sqr(pos[1]));
  e_tot = e_pot + 0.5*(sqr(Uorb(o1,0)) + sqr(Vorb(o1,0)) + sqr(Worb(o1,0)));
  I1(o1) = e_tot;			/* save first integral of motion */
  dprintf(1,"Etot for this orbit = %f\n",e_tot);

  if (l1<=0) return(l1);

  pos[0] = Xorb(o1,l1-1);
  pos[1] = Yorb(o1,l1-1);
  pos[2] = Zorb(o1,l1-1);
  (*pot)(&ndim, pos, acc, &e_pot, &time); /* e_pot; to get e_tot */
  e_pot -= 0.5*omegasq*(sqr(pos[0]) + sqr(pos[1]));
  e_pot += 0.5*(sqr(Uorb(o1,l1-1)) + sqr(Vorb(o1,l1-1)) + sqr(Worb(o1,l1-1)));
  de = (e_pot-I1(o1))/I1(o1);
  IE1(o1) = de;
  dprintf(1,"Etot                = %f after 1st cycle (dE/E=%g)\n",e_pot,de);

  if(ncross>1) return(l1);      /* if interactive SOS mode: don't iterate */

  copy_orbit(o1,o2);                  /* make a copy of this orbit */
  Posorb(o2,0,dirint) += smallstep;   /* now perturb the orbit a bit */
  pos[0] = Xorb(o2,0);
  pos[1] = Yorb(o2,0);
  pos[2] = Zorb(o2,0);
  (*pot)(&ndim,pos,acc,&e_pot,&time);        /* get new e_pot */
  e_pot -= 0.5*omegasq*(sqr(pos[0]) + sqr(pos[1]));
  vnew = 2*(e_tot-e_pot);           /* to get new launch velocity */
  if (vnew < 0.0) {                 /* in order to maintain e_tot */
     vnew =  0.0;
     warning("Resetting vnew=0 at xnew=%f New e_tot=%f\n",xnew,e_pot);
  }
  Velorb(o2,0,1-dirint) = sqrt(vnew) * dir0;

  for(iter=0;;iter++) {                  /* iterate until satisfied or exhausted */
    dprintf(2,"********************ITERATION %d *****************\n",iter+1);
    l2 = (*cycle)(o2);              /* cycle the new estimate (by T/period) */
    dprintf(2," cycle2 @ %d\n",l2);
    dprintf(2,"POS,VEL= %f %f %f %f\n         %f %f %f %f\n",
        Xorb(o2,0), Yorb(o2,0), Uorb(o2,0), Vorb(o2,0),
        Xorb(o2,l2-1), Yorb(o2,l2-1), Uorb(o2,l2-1), Vorb(o2,l2-1));

    x1_s = Posorb(o1,0,dirint);            /* surface of section coord's at t=0 */
    u1_s = Velorb(o1,0,dirint);
    x2_s = Posorb(o2,0,dirint);
    u2_s = Velorb(o2,0,dirint);
    a_s = (u2_s-u1_s)/(x2_s-x1_s);
    b_s = u1_s - x1_s*a_s;

    x1_e = sign*Posorb(o1,l1-1,dirint); /* surface of section coord's at t=T/2 */
    u1_e = sign*Velorb(o1,l1-1,dirint);
    x2_e = sign*Posorb(o2,l2-1,dirint);
    u2_e = sign*Velorb(o2,l2-1,dirint);
    a_e = (u2_e-u1_e)/(x2_e-x1_e);
    b_e = u1_e - x1_e*a_e;

    copy_orbit(o2,o1);                      /* put this one back */
    l1 = l2;

    xnew = -(b_e-b_s)/(a_e-a_s);            /* new estimate for pos */
    Posorb(o2,0,dirint) = xnew;             /* and set a new estimate for PO */
    pos[0] = Xorb(o2,0);
    pos[1] = Yorb(o2,0);
    pos[2] = Zorb(o2,0);
    (*pot)(&ndim,pos,acc,&e_pot,&time);        /* get new e_pot */
    e_pot -= 0.5*omegasq*(sqr(pos[0]) + sqr(pos[1]));
    vnew = 2*(e_tot-e_pot);
    if (vnew < 0.0) {
        vnew =  0.0;
        warning("Resetting vnew=0 at xnew=%g New e_tot=%g\n",xnew,e_pot);
    }
    Velorb(o2,0,1-dirint) = sqrt(vnew) * dir0;

    eps = ABS((xnew-x1_s)/(xnew));
#if 1
    pos[0] = Xorb(o2,l2-1);	/* look at energy conservation one-but-last */
    pos[1] = Yorb(o2,l2-1);
    pos[2] = Zorb(o2,l2-1);
    (*pot)(&ndim,pos,acc,&e_pot,&time);        /* get new e_pot */
    e_pot -= 0.5*omegasq*(sqr(pos[0]) + sqr(pos[1]));
    e_pot += 0.5*(sqr(Uorb(o2,l2-1)) + sqr(Vorb(o2,l2-1))); /* 2D */
#endif
    dprintf(1,"eps = %f xnew,vnew=%f %f e=%f\n",eps,xnew,vnew,e_pot);
    if (eps < eps_min) {                    /* check if good enough in phase space */
        l2 = (*cycle)(o2);                     /* compute best guess */
        dprintf(1,"iters = %d eps=%f\n",iter+1,eps);
        for (i=0, lz_mean=0.0; i<l2; i++)
            lz_mean += Xorb(o2,i)*Vorb(o2,i)-Yorb(o2,i)*Uorb(o2,i)
                       + omega*sqrt(sqr(Xorb(o2,i))+sqr(Yorb(o2,i)));
        lz_mean /= l2;

        h = (period==1 ? l2/4 : l2/2);       /* where to get the other */
	if (first_out) {
	  fprintf(tab,"# %s   %s   %s   %s   NPT   NITER   PERIOD   ETOT   LZ_MEAN   ETOT_ERR\n",
		  dirint==0 ? "x0" : "y0",
		  dirint==0 ? "v0" : "u0",
		  dirint==0 ? "y1" : "x1",
		  dirint==0 ? "u1" : "v1");
	  first_out = FALSE;
	}
        fprintf(tab,"%f %f %f %f %d %d %f %f %f %g\n",
            Posorb(o2,0,dirint), Velorb(o2,0,1-dirint),
            Posorb(o2,h,1-dirint),Velorb(o2,h,dirint),
	    l2,iter+1,period*Torb(o2,l2-1),e_tot,lz_mean, de);
        dprintf(1,"Final orbit #%d: %f %f %f %f %f %f %f %f\n", iorb+1,
           Xorb(o2,0), Yorb(o2,0), Uorb(o2,0), Vorb(o2,0),
           Xorb(o2,l2-1), Yorb(o2,l2-1), Uorb(o2,l2-1), Vorb(o2,l2-1));
        copy_orbit(o2,o1);
        return l2;                         /* and be done with iterations */
    } else if (eps > eps_max)
        break;
    if (iter==maxiter) 
        break;
  }
  return 0;
}


/*------------------------------------------------------------------------------
 *  CYCLE: integrate the orbit for (half an) orbit
 *
 *  It returns the number of steps it took to change sign in the [dirint] coordinate.
 *------------------------------------------------------------------------------
 */
int cycle_euler(orbitptr o_out)
{
        int i, ndim;
        double ax,ay,az,x,y,z,ekin,epot;
        double pos[3],vel[3],acc[3],time;

        double f,xnew,ynew,unew,vnew;
        int    s0,s;

        dprintf (1,"EULER integration\n");
        ndim=Ndim(o_out);               /* number of dimensions (2 or 3) */
        i=1;                            /* count the steps it took */
        s0 = SGN(Velorb(o_out,0,1-dirint));
        icross = 0;   /* counter to store SOS coordinates */
        for(;;) {   /* infinite loop until broken inside */
            pos[0] = Xorb(o_out,i-1);
            pos[1] = Yorb(o_out,i-1);
            pos[2] = Zorb(o_out,i-1);
            time = Torb(o_out,i-1);
            (*pot)(&ndim,pos,acc,&epot,&time);    /* only to get acc */
            /* note epot is not corrected for centrifugal term here */
            if (i>=Nsteps(o_out)) return 0;      /* no space to cycle */

                        /* simple forward euler integration */
                
            Torb(o_out,i) = Torb(o_out,i-1) + dt;

            Xorb(o_out,i) = Xorb(o_out,i-1) + dt*Uorb(o_out,i-1);       
            Yorb(o_out,i) = Yorb(o_out,i-1) + dt*Vorb(o_out,i-1);
            Zorb(o_out,i) = Zorb(o_out,i-1) + dt*Worb(o_out,i-1);

            Uorb(o_out,i) = Uorb(o_out,i-1) + 
                dt*(acc[0]+omegasq*Xorb(o_out,i-1)+omegato*Vorb(o_out,i-1));
            Vorb(o_out,i) = Vorb(o_out,i-1) + 
                dt*(acc[1]+omegasq*Yorb(o_out,i-1)-omegato*Uorb(o_out,i-1));
            Worb(o_out,i) = Worb(o_out,i-1) + dt*acc[2];


                /* the following only works for 2d orbits */
            f = Posorb(o_out,i,1-dirint) * Posorb(o_out,i-1,1-dirint);
            if (f < 0.0) {   /* found  a crossing */
                s = (Velorb(o_out,i,1-dirint) > 0 ? 1 : -1); /* get vel sign */
                if (period==1 && s!=s0) {
		    i++;
                    continue;           /* ai, not the right one yet */
		}
#if 0
                f = -Posorb(o_out,i-1,1-dirint)/
                         (dt*Velorb(o_out,i-1,1-dirint));
#else
                f = Posorb(o_out,i-1,1-dirint)/
                   (Posorb(o_out,i-1,1-dirint)-Posorb(o_out,i,1-dirint));
#endif 
                xnew = 0.0;     /* by def */
                unew = (1-f)*Velorb(o_out,i-1,1-dirint) + f*Velorb(o_out,i,1-dirint);
                ynew = (1-f)*Posorb(o_out,i-1,dirint) + f*Posorb(o_out,i,dirint);
                vnew = (1-f)*Velorb(o_out,i-1,dirint) + f*Velorb(o_out,i,dirint);

                if (unew<0) {
                  xsos[icross] = ynew;        /* save SOS coord's */
                  ysos[icross] = vnew;
                } else {
                  xsos[icross] = -ynew;        /* save SOS coord's */
                  ysos[icross] = -vnew;
                }
                tsos[icross] = Torb(o_out,i-1) + f*dt;
                esos[icross] = epot+(sqr(unew)+sqr(vnew)-omegasq*sqr(ynew))/2;
                if (++icross >= ncross) { /* and finished with cycling around */
                    Posorb(o_out,i,1-dirint) = xnew;
                    Posorb(o_out,i, dirint ) = ynew;
                    Velorb(o_out,i,1-dirint) = unew;
                    Velorb(o_out,i, dirint ) = vnew;
                    Torb(o_out,i) = Torb(o_out,i-1) + f*dt;
                    break;
                } else {
                    move_o(o_out,0,i);  /* copy slot i to 0 */
                    i=1;                /* reset local storage counter i */
                    continue;       /* and keep on cyclin' */
                }
            } /*if(f)*/
            i++;
        } 
        return i+1;
}

int cycle_leapfrog(orbitptr o_out)
{
    int i, ndim;
    double ax,ay,az,x,y,z,ekin,epot;
    double time,pos[3],vel[3],acc[3],acc1[3],vel0;
    double f,xnew,ynew,unew,vnew;
    int    s0, s;

    dprintf(1,"EXPERIMENTAL LEAPFROG integration\n");

            /* start at last step of input file */
    time = Torb(o_out,0);
    pos[0] = Xorb(o_out,0);
    pos[1] = Yorb(o_out,0);
    pos[2] = Zorb(o_out,0);
    vel[0] = Uorb(o_out,0);
    vel[1] = Vorb(o_out,0);
    vel[2] = Worb(o_out,0);
    ndim=Ndim(o_out);               /* number of dimensions (2 or 3) */

    i=1;                            /* counter of timesteps */
    s0 = SGN(Velorb(o_out,0,1-dirint)); /* sense of orbit */
    icross = 0;                     /* counter to store SOS coords */

    (*pot)(&ndim,pos,acc,&epot);
    I1(o_out) = epot + 0.5*(sqr(vel[0]) + sqr(vel[1])+ sqr(vel[2]));

    for(;;) {
        acc1[0] = acc[0];                       /* save old acc's */
        acc1[1] = acc[1];
        acc1[2] = acc[2];
        (*pot)(&ndim,pos,acc,&epot,&time);  /* get new acc's */
        if (i>=Nsteps(o_out)) return 0;     /* no space to cycle */
                                        /* second order correction to vel */
        vel[0] += dt2*(acc[0]-acc1[0]);         
        vel[1] += dt2*(acc[1]-acc1[1]);
        vel[2] += dt2*(acc[2]-acc1[2]);

        /* --> this is where check should be done <-- */
        
        time += dt;
        vel0=vel[0];
        vel[0] += dt2*(acc[0] + omegasq*pos[0] + omegato*vel[1]);
        vel[1] += dt2*(acc[1] + omegasq*pos[1] - omegato*vel0);
        vel[2] += dt2*acc[2];
        pos[0] += dt*vel[0];
        pos[1] += dt*vel[1];
        pos[2] += dt*vel[2];
        vel0 = vel[0];
        vel[0] += dt2*(acc[0] + omegasq*pos[0] + omegato*vel[1]);   
        vel[1] += dt2*(acc[1] + omegasq*pos[1] - omegato*vel0);
        vel[2] += dt2*acc[2];

        Torb(o_out,i) = time;
        Xorb(o_out,i) = pos[0];
        Yorb(o_out,i) = pos[1];
        Zorb(o_out,i) = pos[2];
        Uorb(o_out,i) = vel[0];
        Vorb(o_out,i) = vel[1];
        Worb(o_out,i) = vel[2];

        /* the following only works for 2d orbits */
        
        f = Posorb(o_out,i,1-dirint) * Posorb(o_out,i-1,1-dirint);
        if (f < 0.0) {     /* found a crossing, interpolate accordingly */
            s = Velorb(o_out,i,1-dirint) > 0 ? 1 : -1;  /* vel sign */
            if (period==1 && s!=s0) {
                i++;
                continue;               /* not the right one yet */
            }
#if 0                    
            f = -Posorb(o_out,i-1,1-dirint)/
                 (dt*Velorb(o_out,i-1,1-dirint));
#else
            f = Posorb(o_out,i-1,1-dirint)/
                 (Posorb(o_out,i-1,1-dirint)-Posorb(o_out,i,1-dirint));

#endif                                
            xnew = 0.0;     /* by def */
            unew = (1-f)*Velorb(o_out,i-1,1-dirint) + f*Velorb(o_out,i,1-dirint);
            ynew = (1-f)*Posorb(o_out,i-1,dirint) + f*Posorb(o_out,i,dirint);
            vnew = (1-f)*Velorb(o_out,i-1,dirint) + f*Velorb(o_out,i,dirint);

            if (unew<0) {
                xsos[icross] = ynew;        /* save SOS coord's */
                ysos[icross] = vnew;
            } else {
                xsos[icross] = -ynew;        /* save SOS coord's */
                ysos[icross] = -vnew;
            }
            tsos[icross] = Torb(o_out,i-1) + f*dt;
            esos[icross] = epot+(sqr(unew)+sqr(vnew)-omegasq*sqr(ynew))/2;
            if (++icross >= ncross) { /* and finished with cycling around */
                    Posorb(o_out,i,1-dirint) = xnew;
                    Posorb(o_out,i, dirint ) = ynew;
                    Velorb(o_out,i,1-dirint) = unew;
                    Velorb(o_out,i, dirint ) = vnew;
                    Torb(o_out,i) = Torb(o_out,i-1) + f*dt;
                    break;
            } else {
                    move_o(o_out,0,i);  /* copy slot i to 0 */
                    i=1;                /* reset local storage counter i */
                    continue;           /* and keep on cyclin' */
            }
        } /* if(f) */
        i++;
    }
    return i+1;
}

int cycle_rk2(orbitptr o_out)
{
        int i, ndim;
        double ax,ay,az,x,y,z,ekin,epot;
        double pos[3],vel[3],acc[3],time;
        double k1[6], k2[6];

        double f,xnew,ynew,unew,vnew;
        int    s0,s;

        dprintf (1,"RK2 integration\n");
        ndim=Ndim(o_out);               /* number of dimensions (2 or 3) */
        i=1;                            /* count the steps it took */
        s0 = SGN(Velorb(o_out,0,1-dirint));
        icross = 0;   /* counter to store SOS coordinates */
        for(;;) {   /* infinite loop until broken inside */

            if (i>=Nsteps(o_out)) return 0;      /* no space to cycle */

            pos[0] = Xorb(o_out,i-1);
            pos[1] = Yorb(o_out,i-1);
            pos[2] = Zorb(o_out,i-1);
            vel[0] = Uorb(o_out,i-1);
            vel[1] = Vorb(o_out,i-1);
            vel[2] = Worb(o_out,i-1);
            time   = Torb(o_out,i-1);

            (*pot)(&ndim,pos,acc,&epot,&time);    /* get acc */

                
            Torb(o_out,i) = Torb(o_out,i-1) + dt;

            set_rk(k1,pos,vel,acc,0,dt,k1);
            set_rk(k2,pos,vel,acc,1,dt,k1);

            Xorb(o_out,i) = Xorb(o_out,i-1) + k2[0];
            Yorb(o_out,i) = Yorb(o_out,i-1) + k2[1];
            Zorb(o_out,i) = Zorb(o_out,i-1) + k2[2];

            Uorb(o_out,i) = Uorb(o_out,i-1) + k2[3];
            Vorb(o_out,i) = Vorb(o_out,i-1) + k2[4];
            Worb(o_out,i) = Worb(o_out,i-1) + k2[5];

                /* the following only works for 2d orbits */
            f = Posorb(o_out,i,1-dirint) * Posorb(o_out,i-1,1-dirint);
            if (f < 0.0) {   /* found  a crossing */
                s = (Velorb(o_out,i,1-dirint) > 0 ? 1 : -1); /* get vel sign */
                if (period==1 && s!=s0) {
		    i++;
                    continue;           /* ai, not the right one yet */
		}
#if 0
                f = -Posorb(o_out,i-1,1-dirint)/
                         (dt*Velorb(o_out,i-1,1-dirint));
#else
                f = Posorb(o_out,i-1,1-dirint)/
                   (Posorb(o_out,i-1,1-dirint)-Posorb(o_out,i,1-dirint));
#endif 
                xnew = 0.0;     /* by def */
                unew = (1-f)*Velorb(o_out,i-1,1-dirint) + f*Velorb(o_out,i,1-dirint);
                ynew = (1-f)*Posorb(o_out,i-1,dirint) + f*Posorb(o_out,i,dirint);
                vnew = (1-f)*Velorb(o_out,i-1,dirint) + f*Velorb(o_out,i,dirint);

                if (unew<0) {
                  xsos[icross] = ynew;        /* save SOS coord's */
                  ysos[icross] = vnew;
                } else {
                  xsos[icross] = -ynew;        /* save SOS coord's */
                  ysos[icross] = -vnew;
                }
                tsos[icross] = Torb(o_out,i-1) + f*dt;
                esos[icross] = epot+(sqr(unew)+sqr(vnew)-omegasq*sqr(ynew))/2;
                if (++icross >= ncross) { /* and finished with cycling around */
                    Posorb(o_out,i,1-dirint) = xnew;
                    Posorb(o_out,i, dirint ) = ynew;
                    Velorb(o_out,i,1-dirint) = unew;
                    Velorb(o_out,i, dirint ) = vnew;
                    Torb(o_out,i) = Torb(o_out,i-1) + f*dt;
                    break;
                } else {
                    move_o(o_out,0,i);  /* copy slot i to 0 */
                    i=1;                /* reset local storage counter i */
                    continue;       /* and keep on cyclin' */
                }
            } /*if(f)*/
            i++;
        } 
        return i+1;
}

int cycle_rk4(orbitptr o_out)
{
        int i, ndim;
        double ax,ay,az,x,y,z,ekin,epot;
        double pos[3],vel[3],acc[3],time;
        double k1[6], k2[6], k3[6], k4[6];

        double f,xnew,ynew,unew,vnew;
        int    s0,s;

        dprintf (1,"RK4 integration\n");
        ndim=Ndim(o_out);               /* number of dimensions (2 or 3) */
        i=1;                            /* count the steps it took */
        s0 = SGN(Velorb(o_out,0,1-dirint));
        icross = 0;   /* counter to store SOS coordinates */
        for(;;) {   /* infinite loop until broken inside */

            if (i>=Nsteps(o_out)) return 0;      /* no space to cycle */

            pos[0] = Xorb(o_out,i-1);
            pos[1] = Yorb(o_out,i-1);
            pos[2] = Zorb(o_out,i-1);
            vel[0] = Uorb(o_out,i-1);
            vel[1] = Vorb(o_out,i-1);
            vel[2] = Worb(o_out,i-1);
            time   = Torb(o_out,i-1);

            (*pot)(&ndim,pos,acc,&epot,&time);    /* get acc */

                
            Torb(o_out,i) = Torb(o_out,i-1) + dt;

            set_rk(k1,pos,vel,acc,0,dt,k1);
            set_rk(k2,pos,vel,acc,1,dt,k1);
            set_rk(k3,pos,vel,acc,1,dt,k2);
            set_rk(k4,pos,vel,acc,2,dt,k3);

            Xorb(o_out,i) = Xorb(o_out,i-1) + 
                            (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0;
            Yorb(o_out,i) = Yorb(o_out,i-1) +
                            (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6.0;
            Zorb(o_out,i) = Zorb(o_out,i-1) + 
                            (k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6.0;

            Uorb(o_out,i) = Uorb(o_out,i-1) + 
                            (k1[3] + 2*k2[3] + 2*k3[3] + k4[3])/6.0;
            Vorb(o_out,i) = Vorb(o_out,i-1) + 
                            (k1[4] + 2*k2[4] + 2*k3[4] + k4[4])/6.0;
            Worb(o_out,i) = Worb(o_out,i-1) + 
                            (k1[5] + 2*k2[5] + 2*k3[5] + k4[5])/6.0;

                /* the following only works for 2d orbits */
            f = Posorb(o_out,i,1-dirint) * Posorb(o_out,i-1,1-dirint);
            if (f < 0.0) {   /* found  a crossing */
                s = (Velorb(o_out,i,1-dirint) > 0 ? 1 : -1); /* get vel sign */
                if (period==1 && s!=s0) {
		    i++;
                    continue;           /* ai, not the right one yet */
		}
#if 0
                f = -Posorb(o_out,i-1,1-dirint)/
                         (dt*Velorb(o_out,i-1,1-dirint));
#else
                f = Posorb(o_out,i-1,1-dirint)/
                   (Posorb(o_out,i-1,1-dirint)-Posorb(o_out,i,1-dirint));
#endif 
                xnew = 0.0;     /* by def */
                unew = (1-f)*Velorb(o_out,i-1,1-dirint) + f*Velorb(o_out,i,1-dirint);
                ynew = (1-f)*Posorb(o_out,i-1,dirint) + f*Posorb(o_out,i,dirint);
                vnew = (1-f)*Velorb(o_out,i-1,dirint) + f*Velorb(o_out,i,dirint);

                if (unew<0) {
                  xsos[icross] = ynew;        /* save SOS coord's */
                  ysos[icross] = vnew;
                } else {
                  xsos[icross] = -ynew;        /* save SOS coord's */
                  ysos[icross] = -vnew;
                }
                tsos[icross] = Torb(o_out,i-1) + f*dt;
                esos[icross] = epot+(sqr(unew)+sqr(vnew)-omegasq*sqr(ynew))/2;
                if (++icross >= ncross) { /* and finished with cycling around */
                    Posorb(o_out,i,1-dirint) = xnew;
                    Posorb(o_out,i, dirint ) = ynew;
                    Velorb(o_out,i,1-dirint) = unew;
                    Velorb(o_out,i, dirint ) = vnew;
                    Torb(o_out,i) = Torb(o_out,i-1) + f*dt;
                    break;
                } else {
                    move_o(o_out,0,i);  /* copy slot i to 0 */
                    i=1;                /* reset local storage counter i */
                    continue;       /* and keep on cyclin' */
                }
            } /*if(f)*/
            i++;
        } 
        return i+1;
}

set_rk(double *ko, double *pos, double *vel, double *acc, 
       int n, double dt, double *ki)
{
    double tpos[3], tvel[3], epot, tdum;
    int ndim=3; /* kludge */

    if (n>0) {
        tpos[0] = pos[0] + 0.5*n*ki[0];
        tpos[1] = pos[1] + 0.5*n*ki[1];
        tpos[2] = pos[2] + 0.5*n*ki[2];
	(*pot)(&ndim,tpos,acc,&epot,&tdum);
        tvel[0] = vel[0] + 0.5*n*ki[3];
        tvel[1] = vel[1] + 0.5*n*ki[4];
        tvel[2] = vel[2] + 0.5*n*ki[5];
    } else {
        tpos[0] = pos[0];
        tpos[1] = pos[1];
        tpos[2] = pos[2];
        tvel[0] = vel[0];
        tvel[1] = vel[1];
        tvel[2] = vel[2];
    }
    ko[0] = dt*tvel[0];
    ko[1] = dt*tvel[1];
    ko[2] = dt*tvel[2];
    ko[3] = dt*(acc[0] + omegasq*tpos[0] + omegato*tvel[1]);
    ko[4] = dt*(acc[1] + omegasq*tpos[1] - omegato*tvel[0]);
    ko[5] = dt*(acc[2] + omegasq*tpos[2]);

}

