/*
 *  ORBINT: testroutine to integrate stellar orbits
 *
 *	xx-jul-87	Princeton - first version	P.J. Teuben
 *	16-jul-87	V1.1 small mod in orbit-struct	PJT
 *	20-jul-87	V1.2 potential-interface defined PJT
 *      28-jul-87       V2.0 mod's for new orbit(5) naming  PJT
 *
 *	20-jul-87:
 *      newton0 implementation is preferred now
 *	the two produce the same answer	so from here on we'll 
 *	abandon this program for a while and use newton0 as integrator.
 *
 *	 8-apr-88:	V2.1 updated for new potential(5)	PJT
 *	 2-May-88	V2.2 have a LEAPFROG implementation	PJT
 *				experimenting with improved LEAPFROG
 *	 2-Jun-88	V2.3 new filestruct			PJT
 *      21-nov-90       V2.4 quick fix for Nemo 2.x             PJT
 *	12-jun-91	V2.5 added rotating frame of reference  PJT
 *                           also needed to add time parameter
 *	 7-mar-92	V2.5a happy gcc2.0			pjt
 *	24-may-92	V2.6 default potential from in= orbit itself	PJT
 *       9-jun-92       V2.7 fixed bug in rotating potentials   PJT
 *	18-oct-93       V3.0 new method of get_pattern()	PJT
 *	 3-dec-93 	V3.1 allow ndiag=0 with no output	pjt
 *      27-mar-95       V3.2 fixed leapfrog1 bug for omega != 0 PJT
 *                           introduced rk2 and rk4
 *	19-apr-95	V3.3 all diagnostics output now with dprintf()
 *			     leapfrog2 is now called test
 *       3-feb-98       V3.4 stop criterion if energy not conserved
 *				well enough			PJT
 *      29-oct-00       a    don't report of ndiag=0 given      PJT
 *      10-feb-04       V4.0 variable timestepping              PJT
 *      14-jul-09       V4.1 Bastille Day @ PiTP - added some extra integration modes PJT
 *                           after Tremaines nice lecture
 *      10-dec-2019     V4.2 Add optional Phi/Acc to output     PJT
 *                           but not implemented for all cases - also fixed pattern speed bug
 *                           
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>	/* careful: dangerous with potentials */
#include <orbit.h>

string defv[] = {
    "in=???\n		  input filename (an orbit) ",
    "out=???\n		  output filename (an orbit) ",
    "nsteps=10\n          number of steps",
    "dt=0.1\n             (initial) timestep",
    "ndiag=0\n		  frequency of diagnostics output (0=none)",
    "nsave=1\n		  frequency of storing ",
    "potname=\n	  	  potential name (default from orbit)",
    "potpars=\n	          parameters of potential ",
    "potfile=\n		  extra data-file for potential ",
    "mode=rk4\n           integration method (euler,leapfrog,rk2,rk4)",
    "eta=\n               if used, stop if abs(de/e) > eta",
    "variable=f\n         Use variable timesteps (needs eta=)",
    "VERSION=4.2\n        10-dec-2019 PJT",
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
real   omega, omega2, tomega;  		/* pattern speed */
real   tdum=0.0;                        /* time used in potential() */
real   eta = -1.0;                      /* stop criterion parameter */
bool   Qstop = FALSE;                   /* global flag to stop intgr. */
bool   Qvar;



extern int match(string, string, int *);


proc pot;				/* pointer to the potential */
real print_diag();                      /* returns total energy/hamiltonian */
void setparams(), prepare();
void integrate_euler1(), integrate_euler2(), 
     integrate_leapfrog1(), integrate_leapfrog2(),
     integrate_rk2(), integrate_rk4();
void set_rk(double *ko, double *pos, double *vel, double *acc, 
	    int n, double dt, double *ki);


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
    match(getparam("mode"),"euler leapfrog test rk2 rk4 me end",&imode);
    if (imode==0x01)           // euler
        integrate_euler1();
    else if (imode==0x02)      // leapfrog
        integrate_leapfrog1();
    else if (imode==0x04)      // test
        integrate_leapfrog2();
    else if (imode==0x08)      // rk2
        integrate_rk2();
    else if (imode==0x10)      // rk4
        integrate_rk4();
    else if (imode==0x20)      // me
        integrate_euler2();
    else
        error("imode=0x%x; Illegal integration mode=",imode);

    put_history(outstr);
    write_orbit (outstr,o_out); 
    strclose(outstr);
}

void setparams() 
{
    infile = getparam("in");
    outfile = getparam("out");

    nsteps = getiparam("nsteps");
    dt = getdparam("dt");
    dt2 = 0.5*dt;
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

/* Standard Euler integration */
void integrate_euler1()
{
    int i, ndim, kdiag, ksave, isave;
    double time,epot,e_last;
    double pos[3],vel[3],acc[3], xvel;

    dprintf (1,"EULER integration\n");
    /* take last step of input file and set first step for outfile */
    time = Torb(o_out,0) = Torb(o_in,Nsteps(o_in)-1);
    pos[0] = Xorb(o_out,0) = Xorb(o_in,Nsteps(o_in)-1);
    pos[1] = Yorb(o_out,0) = Yorb(o_in,Nsteps(o_in)-1);
    pos[2] = Zorb(o_out,0) = Zorb(o_in,Nsteps(o_in)-1);
    vel[0] = Uorb(o_out,0) = Uorb(o_in,Nsteps(o_in)-1);
    vel[1] = Vorb(o_out,0) = Vorb(o_in,Nsteps(o_in)-1);
    vel[2] = Worb(o_out,0) = Worb(o_in,Nsteps(o_in)-1);
#ifdef ORBIT_PHI
    Porb(o_out,0)  = Porb(o_in,Nsteps(o_in)-1);
    AXorb(o_out,0) = AXorb(o_in,Nsteps(o_in)-1);
    AYorb(o_out,0) = AYorb(o_in,Nsteps(o_in)-1);
    AZorb(o_out,0) = AZorb(o_in,Nsteps(o_in)-1);
#endif    

    ndim=Ndim(o_in);		/* number of dimensions (2 or 3) */
    kdiag=0;			/* counter for diagnostics output */
    ksave=0;
    isave=0;
    i=0;				/* counter of timesteps */
    for(;;) {
        if (Qstop) break;
	(*pot)(&ndim,pos,acc,&epot,&time);
        if (i==0) I1(o_out) = print_diag(time,pos,vel,epot);
	if (ndiag && kdiag++ == ndiag) {	/* see if output needed */
	    e_last = print_diag(time,pos,vel,epot);
	    kdiag=0;
	}
	if (i>=nsteps) break;           /* see if need to quit looping */

	time   += dt;                     /* advance particle */
        acc[0] += omega2*pos[0] + tomega*vel[1];    /* rotating frame */
        acc[1] += omega2*pos[1] - tomega*vel[0];    /* corrections    */
	

	pos[0] += dt*vel[0];
	pos[1] += dt*vel[1];
	pos[2] += dt*vel[2];

	vel[0] += dt*acc[0];
	vel[1] += dt*acc[1];
	vel[2] += dt*acc[2];
	i++;

	if (++ksave == nsave) {		/* see if need to store particle */
	    ksave=0;
	    isave++;
	    dprintf(2,"writing isave=%d for i=%d\n",isave,i);
            if (isave>=Nsteps(o_out)) error("Storage error EULER");
	    Torb(o_out,isave) = time;
	    Xorb(o_out,isave) = pos[0];
	    Yorb(o_out,isave) = pos[1];
	    Zorb(o_out,isave) = pos[2];
	    Uorb(o_out,isave) = vel[0];
	    Vorb(o_out,isave) = vel[1];
	    Worb(o_out,isave) = vel[2];
#ifdef ORBIT_PHI
	    (*pot)(&ndim,pos,acc,&epot,&time);	    
	    Porb(o_out,isave) = epot   - 0.5*omega2*(pos[0]*pos[0]+pos[1]*pos[1]);
	    AXorb(o_out,isave)= acc[0] + omega2*pos[0] + tomega*vel[1];
	    AYorb(o_out,isave)= acc[1] + omega2*pos[1] - tomega*vel[0];
	    AZorb(o_out,isave)= acc[2];
#endif	    
	}
    } /* for(;;) */
    if (ndiag)          // @todo   make this also print if ndiag=0
      dprintf(0,"Energy conservation: %g\n", ABS((e_last-I1(o_out))/I1(o_out)));
}

/* Modified Euler integration (drift/kick) */
void integrate_euler2()
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
#ifdef ORBIT_PHI
	    Porb(o_out,isave) = epot   - 0.5*omega2*(pos[0]*pos[0]+pos[1]*pos[1]);
	    AXorb(o_out,isave)= acc[0];
	    AYorb(o_out,isave)= acc[1];
	    AZorb(o_out,isave)= acc[2];
#endif	    
	}
	if (i>=nsteps) break;           /* see if need to quit looping */

    } /* for(;;) */
    if (ndiag)
    dprintf(0,"Energy conservation: %g\n", ABS((e_last-I1(o_out))/I1(o_out)));
}

/* Standard Leapfrog */

/* with an attempt to solve the velocities in an implicit scheme */

void integrate_leapfrog1_implicit_test()
{
	int i, kdiag, ksave, isave, ndim;
	double epot,e_last;
	double time,pos[3],vel[3],acc[3];
	double det, det0, det1, vel0, vel1, dvel0, dvel1;

        dprintf(1,"LEAPFROG integration (implicit)\n");

            /* start at last step of input file */
	time = Torb(o_out,0) = Torb(o_in,Nsteps(o_in)-1);
	pos[0] = Xorb(o_out,0) = Xorb(o_in,Nsteps(o_in)-1);
	pos[1] = Yorb(o_out,0) = Yorb(o_in,Nsteps(o_in)-1);
	pos[2] = Zorb(o_out,0) = Zorb(o_in,Nsteps(o_in)-1);
	vel[0] = Uorb(o_out,0) = Uorb(o_in,Nsteps(o_in)-1);
	vel[1] = Vorb(o_out,0) = Vorb(o_in,Nsteps(o_in)-1);
	vel[2] = Worb(o_out,0) = Worb(o_in,Nsteps(o_in)-1);
	ndim=Ndim(o_in);		/* number of dimensions (2 or 3) */
	kdiag=0;			/* reset-counter for diagnostics */
	ksave=0;			/* reset-counter for saving */
	isave=0;			/* save counter */
	i=0;				/* counter of timesteps */
	(*pot)(&ndim,pos,acc,&epot,&tdum);
	I1(o_out) = print_diag(time,pos,vel,epot);
        /* prepare half step for LEAPFROG to get VEL and POS out of sync */
#if 0
	vel[0] += dt2*(acc[0]+omega2*pos[0]+tomega*vel[1]);
	vel[1] += dt2*(acc[1]+omega2*pos[1]-tomega*vel[0]);
#else
	det = 1 + sqr(dt*omega);
	vel0 = vel[0] + dt2*(acc[0]+omega2*pos[0]+omega*vel[1]);
	vel1 = vel[1] + dt2*(acc[1]+omega2*pos[1]-omega*vel[0]);
	vel[0] = (vel0 + dt*omega*vel1)/det;
	vel[1] = (vel1 - dt*omega*vel0)/det;
#endif
	vel[2] += dt2*acc[2];
	while (i<nsteps) {
		i++;
		time += dt;
		pos[0] += dt*vel[0];
		pos[1] += dt*vel[1];
		pos[2] += dt*vel[2];
		(*pot)(&ndim,pos,acc,&epot,&time);
                /* bring back to sync for possible output */
#if 0
        	vel[0] += dt2*(acc[0]+omega2*pos[0]+tomega*vel[1]);
	        vel[1] += dt2*(acc[1]+omega2*pos[1]-tomega*vel[0]);
#else
		vel0 = vel[0] + dt*(acc[0]+omega2*pos[0]+omega*vel[1]);
		vel1 = vel[1] + dt*(acc[1]+omega2*pos[1]-omega*vel[0]);
		dvel0 = ((vel0 + 2*dt*omega*vel1)/det - vel[0])*0.5;
		dvel1 = ((vel1 - 2*dt*omega*vel0)/det - vel[1])*0.5;
		vel[0] += dvel0;	/* half an integration */
		vel[1] += dvel1;
#endif
	        vel[2] += dt2*acc[2];

		if (ndiag && ++kdiag == ndiag) {
		    e_last=print_diag(time,pos,vel,epot);
		    kdiag = 0;
		}
		if (++ksave == nsave) {
		    ksave=0;
		    isave++;
		    dprintf(2,"writing isave=%d for i=%d\n",isave,i);
		    Torb(o_out,isave) = time;
		    Xorb(o_out,isave) = pos[0];
		    Yorb(o_out,isave) = pos[1];
		    Zorb(o_out,isave) = pos[2];
		    Uorb(o_out,isave) = vel[0];
		    Vorb(o_out,isave) = vel[1];
		    Worb(o_out,isave) = vel[2];
#ifdef ORBIT_PHI
		    Porb(o_out,isave) = epot;    // plus centrifugal term?
		    AXorb(o_out,isave)= acc[0];  // not correction for omega ?
		    AYorb(o_out,isave)= acc[1];
		    AZorb(o_out,isave)= acc[2];
#endif	    
		}
                /* put back out of sync */
#if 0
        	vel[0] += dt2*(acc[0]+omega2*pos[0]+tomega*vel[1]);
	        vel[1] += dt2*(acc[1]+omega2*pos[1]-tomega*vel[0]);
#else
		vel[0] += dvel0;	/* fix up the remaining half */
		vel[1] += dvel1;
#endif
	        vel[2] += dt2*acc[2];
	}
    if (ndiag)
    dprintf(0,"Energy conservation: %g\n", ABS((e_last-I1(o_out))/I1(o_out)));
}
/* Standard Leapfrog */
void integrate_leapfrog1()
{
	int i, kdiag, ksave, isave, ndim;
	double epot,e_last;
	double time,pos[3],vel[3],acc[3];

        dprintf(1,"LEAPFROG integration\n");

            /* start at last step of input file */
	time = Torb(o_out,0) = Torb(o_in,Nsteps(o_in)-1);
	pos[0] = Xorb(o_out,0) = Xorb(o_in,Nsteps(o_in)-1);
	pos[1] = Yorb(o_out,0) = Yorb(o_in,Nsteps(o_in)-1);
	pos[2] = Zorb(o_out,0) = Zorb(o_in,Nsteps(o_in)-1);
	vel[0] = Uorb(o_out,0) = Uorb(o_in,Nsteps(o_in)-1);
	vel[1] = Vorb(o_out,0) = Vorb(o_in,Nsteps(o_in)-1);
	vel[2] = Worb(o_out,0) = Worb(o_in,Nsteps(o_in)-1);
#ifdef ORBIT_PHI
	Porb(o_out,0)  = Porb(o_in,Nsteps(o_in)-1);
	AXorb(o_out,0) = AXorb(o_in,Nsteps(o_in)-1);
	AYorb(o_out,0) = AYorb(o_in,Nsteps(o_in)-1);
	AZorb(o_out,0) = AZorb(o_in,Nsteps(o_in)-1);
#endif    
	
	ndim=Ndim(o_in);		/* number of dimensions (2 or 3) */
	kdiag=0;			/* reset-counter for diagnostics */
	ksave=0;			/* reset-counter for saving */
	isave=0;			/* save counter */
	i=0;				/* counter of timesteps */
	(*pot)(&ndim,pos,acc,&epot,&tdum);
	I1(o_out) = print_diag(time,pos,vel,epot);
        /* prepare half step for LEAPFROG to get VEL and POS out of sync */
	vel[0] += dt2*(acc[0]+omega2*pos[0]+tomega*vel[1]);
	vel[1] += dt2*(acc[1]+omega2*pos[1]-tomega*vel[0]);
	vel[2] += dt2*acc[2];
	while (i<nsteps) {
		i++;
		time += dt;
		pos[0] += dt*vel[0];
		pos[1] += dt*vel[1];
		pos[2] += dt*vel[2];
		(*pot)(&ndim,pos,acc,&epot,&time);
                /* bring back to sync for possible output */
        	vel[0] += dt2*(acc[0]+omega2*pos[0]+tomega*vel[1]);
	        vel[1] += dt2*(acc[1]+omega2*pos[1]-tomega*vel[0]);
	        vel[2] += dt2*acc[2];

		if (ndiag && ++kdiag == ndiag) {
		    e_last=print_diag(time,pos,vel,epot);
		    kdiag = 0;
		}
		if (++ksave == nsave) {
		    ksave=0;
		    isave++;
		    dprintf(2,"writing isave=%d for i=%d\n",isave,i);
		    Torb(o_out,isave) = time;
		    Xorb(o_out,isave) = pos[0];
		    Yorb(o_out,isave) = pos[1];
		    Zorb(o_out,isave) = pos[2];
		    Uorb(o_out,isave) = vel[0];
		    Vorb(o_out,isave) = vel[1];
		    Worb(o_out,isave) = vel[2];
		}
                /* put back out of sync */
        	vel[0] += dt2*(acc[0]+omega2*pos[0]+tomega*vel[1]);
	        vel[1] += dt2*(acc[1]+omega2*pos[1]-tomega*vel[0]);
	        vel[2] += dt2*acc[2];

	}
    if (ndiag)
      dprintf(0,"Energy conservation: %g\n", ABS((e_last-I1(o_out))/I1(o_out)));
}

/* Modified Leapfrog, as used in hackcode1 */
/* it is really the same as leapfrog1, except 2nd order corrections */
/* to the velocities are applied to correct for the fact that the */
/* order of things is odd */

void integrate_leapfrog2()
{
	int i, kdiag, ksave, isave, ndim;
	double epot,e_last;
	double time,pos[3],vel[3],acc[3];
	double acc1[3];

        dprintf(1,"EXPERIMENTAL LEAPFROG integration\n");
	if (omega!=0.0) warning("LEAPFROG2; cannot do omega!=0.0 (%g)",omega);

            /* start at last step of input file */
	time = Torb(o_out,0) = Torb(o_in,Nsteps(o_in)-1);
	pos[0] = Xorb(o_out,0) = Xorb(o_in,Nsteps(o_in)-1);
	pos[1] = Yorb(o_out,0) = Yorb(o_in,Nsteps(o_in)-1);
	pos[2] = Zorb(o_out,0) = Zorb(o_in,Nsteps(o_in)-1);
	vel[0] = Uorb(o_out,0) = Uorb(o_in,Nsteps(o_in)-1);
	vel[1] = Vorb(o_out,0) = Vorb(o_in,Nsteps(o_in)-1);
	vel[2] = Worb(o_out,0) = Worb(o_in,Nsteps(o_in)-1);
	ndim=Ndim(o_in);		/* number of dimensions (2 or 3) */
	kdiag=0;			/* reset-counter for diagnostics */
	ksave=0;			/* reset-counter for saving */
	isave=0;			/* save counter */
	i=0;				/* counter of timesteps */
	(*pot)(&ndim,pos,acc,&epot,&tdum);
	I1(o_out) = print_diag(time,pos,vel,epot);
	while (i<nsteps) {
		acc1[0] = acc[0];	/* save old acc's for */
		acc1[1] = acc[1];	/* a much needed 2nd order */
		acc1[2] = acc[2];	/* correction - see below */

		(*pot)(&ndim,pos,acc,&epot,&tdum);  /* get new acc's */

		vel[0] += dt2*(acc[0]-acc1[0]);	   /* second order */
		vel[1] += dt2*(acc[1]-acc1[1]);	   /* correction */
		vel[2] += dt2*(acc[2]-acc1[2]);
#if 1
		if (ndiag && ++kdiag == ndiag) {
			e_last = print_diag(time,pos,vel,epot);
			kdiag = 0;
		}
#endif

		time += dt;
		i++;
		vel[0] += dt2*acc[0];
		vel[1] += dt2*acc[1];
		vel[2] += dt2*acc[2];
		pos[0] += dt*vel[0];
		pos[1] += dt*vel[1];
		pos[2] += dt*vel[2];
		vel[0] += dt2*acc[0];	/* should be new acc's, but */
		vel[1] += dt2*acc[1];	/* we don't have those yet */
		vel[2] += dt2*acc[2];	/* and need to correct next */

		/* these are not good to save, the 2nd order correction */
		/* has not been applied to Vel */

		if (++ksave == nsave) {
			ksave=0;
			isave++;
			dprintf(2,"writing isave=%d for i=%d\n",isave,i);
			Torb(o_out,isave) = time;
			Xorb(o_out,isave) = pos[0];
			Yorb(o_out,isave) = pos[1];
			Zorb(o_out,isave) = pos[2];
			Uorb(o_out,isave) = vel[0];
			Vorb(o_out,isave) = vel[1];
			Worb(o_out,isave) = vel[2];
		}
	}
    if (ndiag)
    dprintf(0,"Energy conservation: %g\n", ABS((e_last-I1(o_out))/I1(o_out)));
}

/* Standard 2nd order RK2 integration */
/* wow: I did this right on the first trial !! 27-mar-95 !! */

void integrate_rk2()
{
    int i, ndim, kdiag, ksave, isave;
    double time,epot,e_last;
    double pos[3],vel[3],acc[3], xvel;
    double k1[6], k2[6];

    dprintf (1,"RK2 integration\n");
    /* take last step of input file and set first step for outfile */
    time = Torb(o_out,0) = Torb(o_in,Nsteps(o_in)-1);
    pos[0] = Xorb(o_out,0) = Xorb(o_in,Nsteps(o_in)-1);
    pos[1] = Yorb(o_out,0) = Yorb(o_in,Nsteps(o_in)-1);
    pos[2] = Zorb(o_out,0) = Zorb(o_in,Nsteps(o_in)-1);
    vel[0] = Uorb(o_out,0) = Uorb(o_in,Nsteps(o_in)-1);
    vel[1] = Vorb(o_out,0) = Vorb(o_in,Nsteps(o_in)-1);
    vel[2] = Worb(o_out,0) = Worb(o_in,Nsteps(o_in)-1);
#ifdef ORBIT_PHI
    Porb(o_out,0)  = Porb(o_in,Nsteps(o_in)-1);
    AXorb(o_out,0) = AXorb(o_in,Nsteps(o_in)-1);
    AYorb(o_out,0) = AYorb(o_in,Nsteps(o_in)-1);
    AZorb(o_out,0) = AZorb(o_in,Nsteps(o_in)-1);
#endif    

    ndim=Ndim(o_in);		/* number of dimensions (2 or 3) */
    kdiag=0;			/* counter for diagnostics output */
    ksave=0;
    isave=0;
    i=0;				/* counter of timesteps */
    for(;;) {
        if (Qstop) break;        
	(*pot)(&ndim,pos,acc,&epot,&time);
        if (i==0) I1(o_out) = print_diag(time,pos,vel,epot);
	if (ndiag && kdiag++ == ndiag) {	/* see if output needed */
	    e_last = print_diag(time,pos,vel,epot);
	    kdiag=0;
	}
	if (i>=nsteps) break;           /* see if need to quit looping */

	time += dt;                     /* advance particle */

        set_rk(k1,pos,vel,acc,0,dt,k1);
        set_rk(k2,pos,vel,acc,1,dt,k1);

	pos[0] += k2[0];
	pos[1] += k2[1];
	pos[2] += k2[2];

	vel[0] += k2[3];
	vel[1] += k2[4];
	vel[2] += k2[5];

	i++;

	if (++ksave == nsave) {		/* see if need to store particle */
	    ksave=0;
	    isave++;
	    dprintf(2,"writing isave=%d for i=%d\n",isave,i);
            if (isave>=Nsteps(o_out)) error("Storage error RK2");
	    Torb(o_out,isave) = time;
	    Xorb(o_out,isave) = pos[0];
	    Yorb(o_out,isave) = pos[1];
	    Zorb(o_out,isave) = pos[2];
	    Uorb(o_out,isave) = vel[0];
	    Vorb(o_out,isave) = vel[1];
	    Worb(o_out,isave) = vel[2];
#ifdef ORBIT_PHI
	    (*pot)(&ndim,pos,acc,&epot,&time);
	    Porb(o_out,isave) = epot   - 0.5*omega2*(pos[0]*pos[0]+pos[1]*pos[1]);
	    AXorb(o_out,isave)= acc[0] + omega2*pos[0] + tomega*vel[1];
	    AYorb(o_out,isave)= acc[1] + omega2*pos[1] - tomega*vel[0];
	    AZorb(o_out,isave)= acc[2];
#endif	    
	}
    } /* for(;;) */
    if (ndiag)
    dprintf(0,"Energy conservation: %g\n", ABS((e_last-I1(o_out))/I1(o_out)));
}

/* 4th order RK4 integration */

void integrate_rk4()
{
    int i, ndim, kdiag, ksave, isave;
    double time,epot,e_last;
    double pos[3],vel[3],acc[3], xvel;
    double k1[6], k2[6], k3[6], k4[6];

    dprintf (1,"RK4 integration\n");
    /* take last step of input file and set first step for outfile */
    time = Torb(o_out,0) = Torb(o_in,Nsteps(o_in)-1);
    pos[0] = Xorb(o_out,0) = Xorb(o_in,Nsteps(o_in)-1);
    pos[1] = Yorb(o_out,0) = Yorb(o_in,Nsteps(o_in)-1);
    pos[2] = Zorb(o_out,0) = Zorb(o_in,Nsteps(o_in)-1);
    vel[0] = Uorb(o_out,0) = Uorb(o_in,Nsteps(o_in)-1);
    vel[1] = Vorb(o_out,0) = Vorb(o_in,Nsteps(o_in)-1);
    vel[2] = Worb(o_out,0) = Worb(o_in,Nsteps(o_in)-1);
#ifdef ORBIT_PHI
    Porb(o_out,0)  = Porb(o_in,Nsteps(o_in)-1);
    AXorb(o_out,0) = AXorb(o_in,Nsteps(o_in)-1);
    AYorb(o_out,0) = AYorb(o_in,Nsteps(o_in)-1);
    AZorb(o_out,0) = AZorb(o_in,Nsteps(o_in)-1);
#endif    

    ndim=Ndim(o_in);		/* number of dimensions (2 or 3) */
    kdiag=0;			/* counter for diagnostics output */
    ksave=0;
    isave=0;
    i=0;				/* counter of timesteps */
    for(;;) {
        if (Qstop) break;
	(*pot)(&ndim,pos,acc,&epot,&time);
        if (i==0) I1(o_out) = print_diag(time,pos,vel,epot);
	if (ndiag && kdiag++ == ndiag) {	/* see if output needed */
	    e_last = print_diag(time,pos,vel,epot);
	    kdiag=0;
	}
	if (i>=nsteps) break;           /* see if need to quit looping */

	time += dt;                     /* advance particle */

        set_rk(k1,pos,vel,acc,0,dt,k1);
        set_rk(k2,pos,vel,acc,1,dt,k1);
        set_rk(k3,pos,vel,acc,1,dt,k2);
        set_rk(k4,pos,vel,acc,2,dt,k3);

	pos[0] += (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0;
	pos[1] += (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6.0;
	pos[2] += (k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6.0;

	vel[0] += (k1[3] + 2*k2[3] + 2*k3[3] + k4[3])/6.0;
	vel[1] += (k1[4] + 2*k2[4] + 2*k3[4] + k4[4])/6.0;
	vel[2] += (k1[5] + 2*k2[5] + 2*k3[5] + k4[5])/6.0;

	i++;

	if (++ksave == nsave) {		/* see if need to store particle */
	    ksave=0;
	    isave++;
	    dprintf(2,"writing isave=%d for i=%d\n",isave,i);
            if (isave>=Nsteps(o_out)) error("Storage error RK4");
	    Torb(o_out,isave) = time;
	    Xorb(o_out,isave) = pos[0];
	    Yorb(o_out,isave) = pos[1];
	    Zorb(o_out,isave) = pos[2];
	    Uorb(o_out,isave) = vel[0];
	    Vorb(o_out,isave) = vel[1];
	    Worb(o_out,isave) = vel[2];
#ifdef ORBIT_PHI
	    (*pot)(&ndim,pos,acc,&epot,&time);
	    Porb(o_out,isave) = epot - 0.5*omega2*(pos[0]*pos[0]+pos[1]*pos[1]);
	    AXorb(o_out,isave)= acc[0] + omega2*pos[0] + tomega*vel[1];
	    AYorb(o_out,isave)= acc[1] + omega2*pos[1] - tomega*vel[0];
	    AZorb(o_out,isave)= acc[2];
#endif	    
	}
    } /* for(;;) */
    if (ndiag)
    dprintf(0,"Energy conservation: %g\n", ABS((e_last-I1(o_out))/I1(o_out)));
}

void set_rk(double *ko, double *pos, double *vel, double *acc, 
	    int n, double dt, double *ki)
{
    double tpos[3], tvel[3], epot, tdum;
    int ndim=3;

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
    ko[3] = dt*(acc[0] + omega2*tpos[0] + tomega*tvel[1]);
    ko[4] = dt*(acc[1] + omega2*tpos[1] - tomega*tvel[0]);
    ko[5] = dt*acc[2];
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
	
    ekin = 0.5*(sqr(vel[0]) + sqr(vel[1])+ sqr(vel[2]));
    epot -= 0.5*omega2*(sqr(pos[0]) + sqr(pos[1]));   
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
