/*
 *  ORBINTV: integrate stellar orbits with variable timestep
 *           for the harder problems
 *
 *      15-may-2011    Cloned off orbint               Peter Teuben
 *      11-jun-2016    fixed *** stack smashing detected ***
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <orbit.h>

#include <dopri5.h>
#include <dop853.h>

string defv[] = {
  "in=???\n		input filename (an orbit) ",
  "out=???\n		output filename (an orbit) ",
  "dtout=0.1\n          Output step time",
  "tstop=10\n           Stop time",
  "dt=0.01\n            (initial) timestep (not used)",
  "ndiag=0\n		** frequency of diagnostics output (0=none)",
  "potname=\n	  	potential name (default from orbit)",
  "potpars=\n	        parameters of potential ",
  "potfile=\n		extra data-file for potential ",
  "mode=dopri5\n        integration method (dopri5 dop853)",
  "tol=-7\n             tolerance of integration",
  "VERSION=1.1a\n       11-jun-2016 PJT",
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

int    nsteps, ndiag;          		/* how often */
real   dt, dtout;			/* stepping */
real   tstop;                           /* stop time */
real   omega, omega2, tomega;  		/* pattern speed */
real   tdum=0.0;                        /* time used in potential() */
real   eta ;                            /* tolerance */



extern int match(string, string, int *);


proc pot;				/* pointer to the potential */

void integrate_dopri5(void),
     integrate_dop853(void),
     setparams(void), 
     prepare(void);
real print_diag(double time, double *posvel);


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

    if (allocate_orbit (&o_out,Ndim(o_in),nsteps)==0)
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

void setparams(void) 
{
    infile = getparam("in");
    outfile = getparam("out");

    nsteps = 1000;
    dt = getdparam("dt");
    dtout = getdparam("dtout");
    tstop = getdparam("tstop");
    nsteps = (int)(tstop/dtout + 1.001);
    dprintf(0,"nsteps = %d\n",nsteps);

    ndiag=getiparam("ndiag");
    eta = getdparam("tol");
    if (eta < 0) 
      eta = pow(10.0,eta);
}

void prepare(void)
{
    Masso(o_out) = Masso(o_in);
}

/* helper functions for the integrator */

static double last_pot;

void rhs(unsigned n, double x, double *y, double *f)
{

  double pos[3], vel[3], acc[3], time;
  int ndim = 3;

  if (n!=6) error("not 3D");

  pos[0] = y[0];
  pos[1] = y[1];
  pos[2] = y[2];
  vel[0] = y[3];
  vel[1] = y[4];
  vel[2] = y[5];
  time = x;
  (*pot)(&ndim,pos,acc,&last_pot,&time);
  acc[0] += omega2*pos[0] + tomega*vel[1];    /* rotating frame */
  acc[1] += omega2*pos[1] - tomega*vel[0];    /* corrections    */

  f[0] = y[3];
  f[1] = y[4];
  f[2] = y[5];
  f[3] = acc[0];
  f[4] = acc[1];
  f[5] = acc[2];

}

void solout5(long nr, double xold, double x, double *y, unsigned n, int *irtrn)
{
  static double xout;
  static int isave;
  double pv[6];

  if (nr==1) {
    xout = x + dtout;
    isave = 0;
    Torb(o_out,isave) = x;
    Xorb(o_out,isave) = y[0]; Yorb(o_out,isave) = y[1]; Zorb(o_out,isave) = y[2];
    Uorb(o_out,isave) = y[3]; Vorb(o_out,isave) = y[4]; Worb(o_out,isave) = y[5];
    (void) print_diag(x, y);
  } else  {
    while (x >= xout) {
      isave += 1;
      Torb(o_out,isave) = xout;
      pv[0] = Xorb(o_out,isave) = contd5(0,xout);
      pv[1] = Yorb(o_out,isave) = contd5(1,xout);
      pv[2] = Zorb(o_out,isave) = contd5(2,xout);
      pv[3] = Uorb(o_out,isave) = contd5(3,xout);
      pv[4] = Vorb(o_out,isave) = contd5(4,xout);
      pv[5] = Worb(o_out,isave) = contd5(5,xout);
      (void) print_diag(xout, pv);
      xout += dtout;
    }
  }
}


void integrate_dopri5(void)
{
  int i, ndim, kdiag, ksave, isave, res, iout, itoler;
  double time,epot,e_last, rtoler, atoler;
  double y[6];
  
  dprintf (1,"DOPRI5 integration\n");
  /* take last step of input file and set first step for outfile */
  time = Torb(o_out,0) = Torb(o_in,Nsteps(o_in)-1);
  y[0] = Xorb(o_out,0) = Xorb(o_in,Nsteps(o_in)-1);
  y[1] = Xorb(o_out,0) = Yorb(o_in,Nsteps(o_in)-1);
  y[2] = Zorb(o_out,0) = Zorb(o_in,Nsteps(o_in)-1);
  y[3] = Uorb(o_out,0) = Uorb(o_in,Nsteps(o_in)-1);
  y[4] = Vorb(o_out,0) = Vorb(o_in,Nsteps(o_in)-1);
  y[5] = Worb(o_out,0) = Worb(o_in,Nsteps(o_in)-1);
  
  ndim=Ndim(o_in);		/* number of dimensions (2 or 3) */
  
  iout = 2;     /* controlling solout */
  itoler = 0;
  rtoler = eta;
  atoler = rtoler;
  res = dopri5(6, rhs, 0.0, y, tstop, &rtoler, &atoler, itoler, solout5, iout,
	       stdout, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 6, NULL, 
	       3);
}

void solout8(long nr, double xold, double x, double *y, unsigned n, int *irtrn)
{
  static double xout;
  static int isave;
  double pv[6];

  if (nr==1) {
    xout = x + dtout;
    isave = 0;
    Torb(o_out,isave) = x;
    Xorb(o_out,isave) = y[0]; Yorb(o_out,isave) = y[1];  Zorb(o_out,isave) = y[2];
    Uorb(o_out,isave) = y[3]; Vorb(o_out,isave) = y[4];  Worb(o_out,isave) = y[5];
    (void) print_diag(x, y);
  } else  {
    while (x >= xout) {
      isave += 1;
      Torb(o_out,isave) = xout;
      pv[0] = Xorb(o_out,isave) = contd8(0,xout);
      pv[1] = Yorb(o_out,isave) = contd8(1,xout);
      pv[2] = Zorb(o_out,isave) = contd8(2,xout);
      pv[3] = Uorb(o_out,isave) = contd8(3,xout);
      pv[4] = Vorb(o_out,isave) = contd8(4,xout);
      pv[5] = Worb(o_out,isave) = contd8(5,xout);
      (void) print_diag(xout, pv);
      xout += dtout;
    }
  }
}

void integrate_dop853()
{
  int i, ndim, kdiag, ksave, isave, res, iout, itoler;
  double time,epot,e_last, rtoler, atoler;
  double y[6];
  
  dprintf (1,"DOP853 integration\n");
  /* take last step of input file and set first step for outfile */
  time = Torb(o_out,0) = Torb(o_in,Nsteps(o_in)-1);
  y[0] = Xorb(o_out,0) = Xorb(o_in,Nsteps(o_in)-1);
  y[1] = Xorb(o_out,0) = Yorb(o_in,Nsteps(o_in)-1);
  y[2] = Zorb(o_out,0) = Zorb(o_in,Nsteps(o_in)-1);
  y[3] = Uorb(o_out,0) = Uorb(o_in,Nsteps(o_in)-1);
  y[4] = Vorb(o_out,0) = Vorb(o_in,Nsteps(o_in)-1);
  y[5] = Worb(o_out,0) = Worb(o_in,Nsteps(o_in)-1);
  
  ndim=Ndim(o_in);		/* number of dimensions (2 or 3) */
  
  iout = 2;     /* controlling solout */
  itoler = 0;
  rtoler = eta;
  atoler = rtoler;
  res = dop853(6, rhs, 0.0, y, tstop, &rtoler, &atoler, itoler, solout8, iout,
	       stdout, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 6, NULL, 
	       3);
}


/*
 *	PRINT_DIAG: print diagnostics, also add centrifugal term
 *		    returns the total energy in the rotating frame
 */

real print_diag(double time, double *posvel)
{
    double ekin;
    double epot;
    permanent bool first = TRUE;
    permanent double etot_0, etot_00;
    permanent int idiag = 0;
    double err;
    double f[6];
    double *pos = posvel;
    double *vel = pos+3;

    rhs(6,time,posvel,f);
    epot = last_pot;
	
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
    if (ndiag) {
      if (idiag%ndiag == 0)
	dprintf(0,"%f %f %f %20.13g %g\n",time,ekin,epot,ekin+epot,err);
    }
    idiag++;
    return ekin+epot;
}
