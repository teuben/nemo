/*
 *  ORBINTV: integrate stellar orbits with variable timestep
 *           for the harder problems
 *
 *      15-may-2011    Cloned off orbint               Peter Teuben
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <orbit.h>

#include <dopri5.h>

string defv[] = {
    "in=???\n		  input filename (an orbit) ",
    "out=???\n		  output filename (an orbit) ",
    "dt=0.01\n            (initial) timestep",
    "tstop=10\n           Stop time",
    "ndiag=0\n		  frequency of diagnostics output (0=none)",
    "nsave=1\n		  frequency of storing ",
    "potname=\n	  	  potential name (default from orbit)",
    "potpars=\n	          parameters of potential ",
    "potfile=\n		  extra data-file for potential ",
    "mode=dopri5\n        integration method (dopri5 dop853)",
    "eta=\n               if used, stop if abs(de/e) > eta",
    "variable=f\n         Use variable timesteps (needs eta=)",
    "VERSION=1.0\n        15-may-11 PJT",
    NULL,
};

string usage = "integrate stellar orbits";

string cvsid="$Id$";


#ifndef HUGE
#  define HUGE 1.0e20
#endif

string	infile, outfile;		/* file names */
stream  instr, outstr;			/* file streams */

orbitptr o_in  = NULL;			/* pointer to input orbit */
orbitptr o_out = NULL;			/* pointer to output orbit */

int    nsteps, ndiag, nsave;		/* how often */
real   dt, dt2;				/* stepping */
real   tstop;                           /* stop time */
real   omega, omega2, tomega;  		/* pattern speed */
real   tdum=0.0;                        /* time used in potential() */
real   eta = -1.0;                      /* stop criterion parameter */
bool   Qstop = FALSE;                   /* global flag to stop intgr. */
bool   Qvar;



extern int match(string, string, int *);


proc pot;				/* pointer to the potential */
real print_diag();                      /* returns total energy/hamiltonian */
void setparams(), prepare();
void integrate_dopri5(),
     integrate_dop853();


/*----------------------------------------------------------------------------*/
void nemo_main ()
{
    int imode;

    setparams();
							       /* open files */
    instr = stropen (infile,"r");

    get_history(instr);			
    if (read_orbit(instr,&o_in)==0)   	 /* read input orbit */
		error ("error in reading input orbit");
    strclose(instr);
    /* replace default orbit with user supplied, if given */
    if (hasvalue("potname")) PotName(o_in) = getparam("potname");
    if (hasvalue("potpars")) PotPars(o_in) = getparam("potpars");
    if (hasvalue("potfile")) PotFile(o_in) = getparam("potfile");

    if (allocate_orbit (&o_out,Ndim(o_in),nsteps/nsave+1)==0)
		error ("Error allocating output orbit");
    pot=get_potential(PotName(o_in), PotPars(o_in), PotFile(o_in));
    if (pot==NULL) 
		error("Potential %s could not be loaded",PotName(o_in));
    PotName(o_out) =  PotName(o_in);
    PotPars(o_out) =  PotPars(o_in);
    PotFile(o_out) =  PotFile(o_in);

    omega = get_pattern();
    dprintf(0,"Pattern speed=%g\n",omega);
    omega2 = omega*omega;
    tomega = 2.0*omega;


    outstr = stropen (outfile,"w");
    prepare();
    match(getparam("mode"),"dopri5 dop853 end",&imode);
    if (imode==0x01)
        integrate_dopri5();
    else if (imode==0x02)
        integrate_dop853();
    else
        error("imode=0x%x; Illegal integration mode=",imode);

    put_history(outstr);
    write_orbit (outstr,o_out); 		/* write output file */
    strclose(outstr);
}

void setparams() 
{
    infile = getparam("in");
    outfile = getparam("out");

    nsteps = 1000;
    dt = getdparam("dt");
    dt2 = 0.5*dt;
    tstop = getdparam("tstop");
    ndiag=getiparam("ndiag");
    nsave=getiparam("nsave");
    Qvar = getbparam("variable");
    if (hasvalue("eta")) 
      eta = getdparam("eta");
    else if (Qvar)
      error("Variable timesteps choosen, it needs a control parameter eta=");
    
}

void prepare()
{
    Masso(o_out) = Masso(o_in);
}

/* helper functions for the integrator helpers */

void rhs(unsigned n, double x, double *y, double *f)
{

  double pos[3], vel[3], acc[3], epot, time;
  int ndim = 3;

  if (n!=6) error("not 3D");

  pos[0] = y[0];
  pos[1] = y[1];
  pos[2] = y[2];
  vel[0] = y[3];
  vel[1] = y[4];
  vel[2] = y[5];
  time = x;
  (*pot)(&ndim,pos,acc,&epot,&time);
  acc[0] += omega2*pos[0] + tomega*vel[1];    /* rotating frame */
  acc[1] += omega2*pos[1] - tomega*vel[0];    /* corrections    */

  f[0] = y[3];
  f[1] = y[4];
  f[2] = y[5];
  f[3] = acc[0];
  f[4] = acc[1];
  f[5] = acc[2];

}

void solout(long nr, double xold, double x, double *y, unsigned n, int *irtrn)
{
  static double xout;

  if (nr==1) {
    printf("%g   %g %g %g\n", x, y[0], y[1], y[2]);
    xout = x + 0.1;
  } else  {
    while (x >= xout) {
      printf("%g   %g %g %g\n", x, y[0], y[1], y[2]);
      xout += 0.1;
    }
  }
}


void integrate_dopri5()
{
  int i, ndim, kdiag, ksave, isave, res, iout, itoler;
  double time,epot,e_last, rtoler, atoler;
  double pos[3],vel[3],acc[3], xvel;
  double y[6];
  
  dprintf (1,"EULER integration\n");
  /* take last step of input file and set first step for outfile */
  time = Torb(o_out,0) = Torb(o_in,Nsteps(o_in)-1);
  y[0] = pos[0] = Xorb(o_out,0) = Xorb(o_in,Nsteps(o_in)-1);
  y[1] = pos[1] = Yorb(o_out,0) = Yorb(o_in,Nsteps(o_in)-1);
  y[2] = pos[2] = Zorb(o_out,0) = Zorb(o_in,Nsteps(o_in)-1);
  y[3] = vel[0] = Uorb(o_out,0) = Uorb(o_in,Nsteps(o_in)-1);
  y[4] = vel[1] = Vorb(o_out,0) = Vorb(o_in,Nsteps(o_in)-1);
  y[5] = vel[2] = Worb(o_out,0) = Worb(o_in,Nsteps(o_in)-1);
  
  ndim=Ndim(o_in);		/* number of dimensions (2 or 3) */
  
  iout = 1;     /* controlling solout */
  itoler = 0;
  rtoler = 1.0e-7;
  atoler = rtoler;
  res = dopri5(6, rhs, 0.0, y, tstop, &rtoler, &atoler, itoler, solout, iout,
	       stdout, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 6, NULL, 
	       3);

}


void integrate_dop853()
{
    int i, ndim, kdiag, ksave, isave;
    double time,epot,e_last;
    double pos[3],vel[3],acc[3], xvel;

    dprintf (0,"Modified EULER integration\n");
    /* take last step of input file and set first step for outfile */
    time = Torb(o_out,0) = Torb(o_in,Nsteps(o_in)-1);
    pos[0] = Xorb(o_out,0) = Xorb(o_in,Nsteps(o_in)-1);
    pos[1] = Yorb(o_out,0) = Yorb(o_in,Nsteps(o_in)-1);
    pos[2] = Zorb(o_out,0) = Zorb(o_in,Nsteps(o_in)-1);
    vel[0] = Uorb(o_out,0) = Uorb(o_in,Nsteps(o_in)-1);
    vel[1] = Vorb(o_out,0) = Vorb(o_in,Nsteps(o_in)-1);
    vel[2] = Worb(o_out,0) = Worb(o_in,Nsteps(o_in)-1);

    ndim=Ndim(o_in);		/* number of dimensions (2 or 3) */
    kdiag=0;			/* counter for diagnostics output */
    ksave=0;
    isave=0;
    i=0;				/* counter of timesteps */
    for(;;) {
        if (Qstop) break;

	/* first update the positions */

	pos[0] += dt*vel[0];
	pos[1] += dt*vel[1];
	pos[2] += dt*vel[2];

	/* get forces at new positions */

	(*pot)(&ndim,pos,acc,&epot,&time);

	time += dt;                     /* advance particle */
	if (omega != 0) {
	  acc[0] += omega2*pos[0] + tomega*vel[1];    /* rotating frame */
	  acc[1] += omega2*pos[1] - tomega*vel[0];    /* corrections    */
	}

	vel[0] += dt*acc[0];
	vel[1] += dt*acc[1];
	vel[2] += dt*acc[2];
	i++;

        if (i==0) I1(o_out) = print_diag(time,pos,vel,epot);
	if (ndiag && kdiag++ == ndiag) {	/* see if output needed */
	    e_last = print_diag(time,pos,vel,epot);
	    kdiag=0;
	}

	if (++ksave == nsave) {		/* see if need to store particle */
	    ksave=0;
	    isave++;
	    dprintf(2,"writing isave=%d for i=%d\n",isave,i);
            if (isave>=Nsteps(o_out)) error("Storage error modified EULER");
	    Torb(o_out,isave) = time;
	    Xorb(o_out,isave) = pos[0];
	    Yorb(o_out,isave) = pos[1];
	    Zorb(o_out,isave) = pos[2];
	    Uorb(o_out,isave) = vel[0];
	    Vorb(o_out,isave) = vel[1];
	    Worb(o_out,isave) = vel[2];
	}
	if (i>=nsteps) break;           /* see if need to quit looping */

    } /* for(;;) */
    if (ndiag)
    dprintf(0,"Energy conservation: %g\n", ABS((e_last-I1(o_out))/I1(o_out)));
}


/*
 *	PRINT_DIAG: print diagnostics, also add centrifugal term
 *		    returns the total energy in the rotating frame
 */

real print_diag(double time, double *pos, double *vel, double epot)
{
    double ekin;
    permanent bool first = TRUE;
    permanent double etot_0, etot_00;
    double err;
	
    ekin=sqr(vel[0]) + sqr(vel[1])+ sqr(vel[2]);
    ekin *= 0.5;
    epot -= 0.5*omega2*(sqr(pos[0]) + sqr(pos[1])+ sqr(pos[2]));
    if (first) {
        dprintf (1,"time   Etot =   ekin +  epot\n");
	first = FALSE;
	etot_0 = epot + ekin;
        etot_00 = (etot_0 == 0.0 ? 1.0 : etot_0);
    }
    err = (epot+ekin-etot_0)/etot_00;
    if (ndiag) dprintf(0,"%f %f %f %20.13g %g\n",time,ekin,epot,ekin+epot,err);
    if (eta > 0 && ABS(err) > eta) {
        warning("STOPPING: Time=%g E=%g E_0=%g Eta=%g",
                    time,epot+ekin,etot_0,eta);
        Qstop = TRUE;
    }
    return ekin+epot;
}
