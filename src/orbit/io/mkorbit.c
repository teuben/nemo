/*
 * MKORBIT:	generate an orbit from various initial conditions
 *		can also read an N-body snapshot(5) file and take a
 *		a particle from that one. It needs a potential
 *		fitter or tree-builder etc. of that kind.
 *
 *	xx-jul-87	first version			P.J. Teuben
 *	15-jul-87	V1.1 small mod in orbit-struct		PJT
 *	28-jul-87	V2.0 new orbit(5)			PJT
 *	 4-May-88	V2.1 allow energy also to be given	PJT
 *	 2-Jun-88	V2.2 new filestruct, no code change	PJT
 *	 2-mar-92	V2.3 NEMO 2.x helpvec usage etc.	PJT
 *			changed names of keywords, added lz=
 *	 7-mar-92       gcc happy				pjt
 *	24-may-92	new orbit(5) with potential		pjt
 *      20-may-93       V2.5 new z= and vz= keys to make 3D     pjt
 *       6-oct-93       V3.0 added rotating pattern support     pjt
 *      19-oct           3.1 using get_pattern() now            pjt
 *       3-dec           3.2 added headline=                    pjt
 *	17-dec		 3.2a fixed bug when vz != 0		pjt
 *	22-feb-95           b ansi headers
 *	17-apr-95           c compacted header file 		pjt
 *      14-sep-01	    d using potproc_ types		pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <orbit.h>

string defv[] = {
    "out=???\n		  output filename (an orbit)",
    "x=\n		  x-position",
    "y=\n                 y-",
    "z=\n                 z-",
    "vx=\n		  x-velocity",
    "vy=\n                y-",
    "vz=\n                z-",
    "etot=\n		  total energy (ekin+epot)",
    "lz=\n		  angular momentum in z (**ignored**)",
    "time=0.0\n           time",
    "potname=\n		  optional potential(5NEMO)",
    "potpars=\n		  .. with optional parameters",
    "potfile=\n		  .. and optional datafile name",
    "headline=\n          random verbiage",
    "VERSION=3.2b\n       22-feb-94 PJT",
    NULL,
};

string usage = "Make an orbit with from given initial conditions";

string	infile,outfile;			/* file names */
stream  instr,outstr;			/* file streams */

orbitptr optr;
a_potential p;

double x,y,z,u,v,w;
double etot, lz;
double tnow, omega;

void setparams();

void nemo_main ()
{
	setparams();	/* get cmdline stuff and compute x,y,u,v,etot,lz */

	optr = NULL;				/* make an orbit */
	allocate_orbit (&optr, 3, 1);
	Masso(optr)  = 1.0;                     /* and set Mass */
	Torb(optr,0) = tnow;
	Xorb(optr,0) = x;			/* .. positions */
	Yorb(optr,0) = y;
	Zorb(optr,0) = z;
	Uorb(optr,0) = u;			/* .. velocities */
	Vorb(optr,0) = v;
	Worb(optr,0) = w;
	I1(optr) = etot;			/*  energy (zero if not used) */
	I2(optr) = lz;                          /* angular momentum */

	dprintf(0,"pos: %f %f %f  \nvel: %f %f %f  \netot: %f\nlz=%f\n",
		  x,y,z,u,v,w,etot,lz);

	outstr = stropen (outfile,"w");		/* write to file */
	put_history(outstr);
        PotName(optr) = p.name;
        PotPars(optr) = p.pars;
        PotFile(optr) = p.file;
	write_orbit(outstr,optr);
	strclose(outstr);
}

void setparams()
{
    potproc_double pot;
    double pos[3],acc[3],epot;
    int ndim=3;

    p.name = getparam("potname");
    p.pars = getparam("potpars");
    p.file = getparam("potfile");

    outfile = getparam("out");
    if(hasvalue("etot")) {	  /* energy given; calculate missing vx or vy */
 	   etot = getdparam("etot");
           if(!hasvalue("potname"))
              error("No potential given (potname=)");
	   pot = get_potential_double(p.name, p.pars, p.file);
	   if (pot==NULL) 
		error("potential %s cannot be loaded",p.name);
	   omega = get_pattern();
	   dprintf(0,"Using pattern speed = %g\n",omega);

	   /* allow only u or v missing for now */
	   /* if x or y not given, 0.0 will be taken */

	   pos[0] = x = getdparam("x");
	   pos[1] = y = getdparam("y");
	   pos[2] = z = 0.0;			/* force z=0 if E_tot given */
	   (*pot)(&ndim,pos,acc,&epot,&tnow);
	   epot -= 0.5*omega*omega*(x*x+y*y);

           if(!hasvalue("vx")) {                /* vx missing */
		v = getdparam("vy");
		u = 2*(etot-epot) - v*v;
		if (u<0.0)
			error("vy too large for this energy: epot=%g vy(max)=%g\n",
					epot,2*sqrt(etot-epot));
		else
			u = sqrt(u);
		lz = x*v-y*u;	
		if(lz<0) u = -u;     /* force counter clockwise orbit */
	   } else {
                if(hasvalue("vy"))
			error("cannot give e=,  vx= and vy= at the same time");
		u = getdparam("vx");
		v = 2*(etot-epot) - u*u;
		if (v<0.0)
			error("vx too large for this energy: epot=%g vx(max)=%g\n",
					epot,2*sqrt(etot-epot));
		else
			v = sqrt(v);
		lz = x*v-y*u;	
		if(lz<0) v = -v;     /* force counter clockwise orbit */
	   }
	   w = 0.0;
           lz = x*v-y*u;
    } else {
	   if(hasvalue("etot")) warning("Initial value for etot= not used");
	   if(hasvalue("lz")) warning("Initial value for lz= not used");
	   x = getdparam("x");			
	   y = getdparam("y");
	   z = getdparam("z");
	   u = getdparam("vx");
	   v = getdparam("vy");
	   w = getdparam("vz");
	   
           if(hasvalue("potname")) {
    	       pot = get_potential_double(p.name, p.pars, p.file);
	       if (pot==NULL) {
		  warning("potential %s cannot be loaded",p.name);
                  etot = 0.0;
               } else {
	          omega = get_pattern();
                  dprintf(0,"Using pattern speed = %g\n",omega);
	          pos[0] = x;
	          pos[1] = y;
	          pos[2] = z;
	          (*pot)(&ndim,pos,acc,&epot,&tnow);
                  etot = epot + 0.5*(u*u+v*v+w*w) - 0.5*omega*omega*(x*x+y*y);
               }
           } else {
               warning("Potential potname=%s not used; set etot=0.0",p.name);
               etot = 0.0;
           }
           lz = x*v-y*u;
    } 
    tnow = getdparam("time");
    if (hasvalue("headline")) set_headline(getparam("headline"));
}
