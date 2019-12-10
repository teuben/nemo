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
 *      11-feb-02       V4.0 start of the new "+/-" notation    pjt
 *      27-aug-04           a   messing with signs
 *      14-oct-10       V4.1a   finished off the +/- sign       pjt
 *
 * TODO:   allow 'circular' orbit, i.e. the one that locally balances centrifugally
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <orbit.h>

string defv[] = {
    "out=???\n		  output filename (an orbit)",
    "x=\n		  x-position (use + or - to designate sign)",
    "y=\n                 y-",
    "z=\n                 z-",
    "vx=\n		  x-velocity (use + or - to designate sign)",
    "vy=\n                y-",
    "vz=\n                z-",
    "etot=\n		  total energy (ekin+epot)",
    "lz=\n		  sign of angular momentum in z",
    "time=0.0\n           time",
    "potname=\n		  optional potential(5NEMO)",
    "potpars=\n		  .. with optional parameters",
    "potfile=\n		  .. and optional datafile name",
    "headline=\n          random verbiage",
    "VERSION=4.4\n        10-dec-2019 PJT",
    NULL,
};

string usage = "Make an orbit with from given initial conditions";
string cvsid="$Id$";

string	infile,outfile;			/* file names */
stream  instr,outstr;			/* file streams */

orbitptr optr;
a_potential p;

double x,y,z,u,v,w;
double etot, lz;
double tnow, omega;
double acc[3],epot;

int  Dpos, Dvel;

void setparams(void);

void nemo_main (void)
{
  setparams();	/* get cmdline stuff and compute x,y,u,v,etot,lz */

  optr = NULL;				/* make an orbit */
  allocate_orbit (&optr, 3, 1);
  Masso(optr)  = 1.0;                   /* and set Mass */
  Key(optr) = 0;
  Torb(optr,0) = tnow;
  Xorb(optr,0) = x;			/* .. positions */
  Yorb(optr,0) = y;
  Zorb(optr,0) = z;
  Uorb(optr,0) = u;			/* .. velocities */
  Vorb(optr,0) = v;
  Worb(optr,0) = w;
#ifdef ORBIT_PHI  
  Porb(optr,0)  = epot;
  AXorb(optr,0) = acc[0];
  AYorb(optr,0) = acc[1];
  AZorb(optr,0) = acc[2];
#endif    
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

void setparams(void)
{
    potproc_double pot;
    double pos[3],vel2,per;
    int ndim=3;
    int signcount;
    int lzsign;

    p.name = getparam("potname");
    p.pars = getparam("potpars");
    p.file = getparam("potfile");

    signcount = 0;
    Dpos = -1;
    if (streq(getparam("x"),"+") || streq(getparam("x"),"-")) {
      Dpos = 0;
      signcount++;
    }
    if (streq(getparam("y"),"+") || streq(getparam("y"),"-")) {
      Dpos = 1;
      signcount++;
    }
    if (streq(getparam("z"),"+") || streq(getparam("z"),"-")) {
      Dpos = 2;
      signcount++;
    }
    if (signcount > 1) error("Can only set one of x,y,z to + or -");

    signcount = 0;
    Dvel = -1;
    if (streq(getparam("vx"),"+") || streq(getparam("vx"),"-")) {
      Dvel = 0;
      signcount++;
    }
    if (streq(getparam("vy"),"+") || streq(getparam("vy"),"-")) {
      Dvel = 1;
      signcount++;
    }
    if (streq(getparam("vz"),"+") || streq(getparam("vz"),"-")) {
      Dvel = 2;
      signcount++;
    }
    if (signcount > 1) error("Can only set one of vx,vy,vz to + or -");
    dprintf(1,"Dpos=%d Dvel=%d\n",Dpos,Dvel);
    
    tnow = getdparam("time");

    outfile = getparam("out");
    if(hasvalue("etot")) {	  /* energy given; calculate missing vx (for y) or vy (for x) */      
 	   etot = getdparam("etot");
	   dprintf(1,"etot=%g\n", etot);
           if(!hasvalue("potname"))
              error("No potential given (potname=)");
	   pot = get_potential_double(p.name, p.pars, p.file);
	   if (pot==NULL) 
		error("potential %s cannot be loaded",p.name);
	   omega = get_pattern();
	   dprintf(0,"Using pattern speed = %g\n",omega);
	   if (hasvalue("lz")) {
	     lz = getdparam("lz");
	     lzsign = SGN(lz);
	   } else
	     lz = 0;

	   /* launching orbit is perpendicular to the axis for which X or Y is given */

	   pos[0] = x = getdparam("x");
	   pos[1] = y = getdparam("y");
	   pos[2] = z = 0.0;			/* force z=0 if E_tot given */
	   (*pot)(&ndim,pos,acc,&epot,&tnow);
	   dprintf(1,"acc: %g %g %g \n", acc[0],acc[1],acc[2]);
	   epot -= 0.5*omega*omega*(x*x+y*y);

           if(hasvalue("y")) {                /* vx missing */
		v = getdparam("vy");
		u = 2*(etot-epot) - v*v;
		if (u<0.0)
		  error("vy too large for this energy: epot=%g vy(max)=%g\n",
			epot,2*sqrt(etot-epot));
		else
		  u = sqrt(u);
		lz = x*v-y*u;	
		if(lzsign != SGN(lz)) u = -u;     /* switch rotation */
	   } else if(hasvalue("x")) {
		u = getdparam("vx");
		v = 2*(etot-epot) - u*u;
		if (v<0.0)
		  error("vx too large for this energy: epot=%g vx(max)=%g\n",
					epot,2*sqrt(etot-epot));
		else
		  v = sqrt(v);
		lz = x*v-y*u;	
		if(lzsign != SGN(lz)) v = -v;     /* switch rotation */
	   }
	   w = 0.0;
           lz = x*v-y*u;
    } else {
	   if(hasvalue("lz")) warning("Initial value for lz= not used");
	   x = getdparam("x");			
	   y = getdparam("y");
	   z = getdparam("z");
	   dprintf(0,"Dvel=%d\n",Dvel);
	   u = Dvel==0 ? 0 : getdparam("vx");
	   v = Dvel==1 ? 0 : getdparam("vy");
	   w = Dvel==2 ? 0 : getdparam("vz");
	   
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
		  dprintf(1,"acc: %g %g %g \n", acc[0],acc[1],acc[2]);
		  if (Dvel==0) {
		    error("Dvel=0 not implemented");
		  } else if (Dvel==1) {
		    vel2 = -acc[0]*x;
		    v = sqrt(vel2);
		    per = TWO_PI*x/v;
		    dprintf(0,"vel2=%g  p=%g\n",vel2,per);
		  } else if (Dvel==2) {
		    error("Dvel=2 not implemented");
		  } 
		  etot = epot + 0.5*(u*u+v*v+w*w) - 0.5*omega*omega*(x*x+y*y);
               }
           } else {
               warning("Potential potname=%s not used; set etot=0.0",p.name);
               etot = 0.0;
           }
           lz = x*v-y*u;
    } 
    if (hasvalue("headline")) set_headline(getparam("headline"));
}
