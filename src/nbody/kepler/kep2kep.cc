/*  kep2kep.c - ... , main<TOOL>, .... */

/*
 *  kep2kep.c:  Kepler orbit transformations, including TOOL version
 *
 *      Dec 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 *	Mar 1991  -  adapted for NEMO V2 - no local error() and warning() 
 */

#include  "kepler.h"

/*-----------------------------------------------------------------------------
 *  transkepler  --  performs a transformation to the six parameters which
 *                   specify a Kepler orbit (circle, ellipse, parabole or
 *                   hyperbole).
 *                   accepts: total_mass: sum of masses of the two bodies;
 *                            orbit: a pointer to a 6-element array which
 *                                   contains the parameters describing the
 *                                   given orbit, in the representation
 *                                   'oldtype';
 *                            oldtype: one of the allowed types of
 *                                     representation of Kepler orbits,
 *                                     indicated by a capitalized macro which
 *                                     stands for an integer (see  kepler.h );
 *                            newtype: another allowed type of representation
 *                                     of Kepler orbits, indicating the desired
 *                                     new type of orbit to transform to.
 *                   note: the new orbit parameters appear in the same array
 *                          orbit[] , in which the old parameters were given.
 *-----------------------------------------------------------------------------
 */
void  transkepler(total_mass, orbit, oldtype, newtype)
real  total_mass;
realptr  orbit;
int  oldtype;
int  newtype;
    {
    register int  k;
    KEPLER  two_body;

    KTime( &two_body ) = 0.0;                           /* useful convention */
/*
 * for historical reasons, concerning the form of data representation, the
 * procedure  transkep()  below requires a data structure which separately
 * contains both masses. Since really only the sum of the two masses is 
 * required for performing Kepler transformations, we use the simple trick 
 * of assigning half the total mass to mass1, and the other half to mass2
 * (any other division would have given the same result):
 */
    if (total_mass < 0.0)
	error("transkepler: total_mass = %lg is negative\n", total_mass);
    else
        KMass1( &two_body ) = KMass2( &two_body ) = 0.5 * total_mass;

    if (iskeptype(oldtype) && iskeptype(newtype))
	KType( &two_body ) = oldtype;
    else
	{
	if (! iskeptype(oldtype))
	    error("transkepler: oldtype = %d is not a valid Kepler type\n",
                  oldtype);
	else
	    error("transkepler: newtype = %d is not a valid Kepler type\n",
                  newtype);
	}

    for (k = 0; k < 2*NDIM; k++)
	KSpace( &two_body )[k] = orbit[k];

    transkep( &two_body, newtype);

    for (k = 0; k < 2*NDIM; k++)
	orbit[k] = KSpace( &two_body )[k];
    }

/*===========================================================================*/

#ifdef  TOOLBOX

/*-----------------------------------------------------------------------------
 *  main  --  interactive transformations between kepler orbit coordinates,
 *	      with a double purpose:
 *            -- to debug the function    transkep(,)   interactively;
 *            -- to form the first part of a Celestial Mechanics Calculator.
 *            You are asked to provide the type of initial coordinates of a
 *            kepler orbit; its masses and the values of the six phase space
 *            coordinates; and then you can transform between any of the
 *            given types of coordinate systems.
 *-----------------------------------------------------------------------------
 */
main()
	{
	KEPLER  comet;
	int  i;

	printf("type in an integer to indicate the type of initial orbit.\n");
	printf("Your choices are:\n");
	printf("   11: CARTESIAN  ");
	printf("     12: SPHERICAL ");
	printf("     31: SCATTERING ");
	printf("     32: TSCATTERING\n");
	printf("   21: TRUEANOMALY");
	printf("     22: ECCANOMALY");
	printf("     23: MEANANOMALY");
	printf("     24: PERIPASSAGE\n");
	scanf("%d", &KType(&comet) );
	if (iskeptype(KType(&comet)))
	        kepinit( &comet );	   /*  initializes the kepler orbit  */
	else
		error("main(): type %d is an unknown type kepler orbit\n", 
		      KType(&comet));

	for (;;)		/*  endless loop  */
		{
		kepprint( &comet );
		printf("type in an integer to get the corresponding");
		printf(" transformation. Your choices:\n");
	        printf("   11: CARTESIAN  ");
	        printf("     12: SPHERICAL ");
	        printf("     31: SCATTERING ");
	        printf("     32: TSCATTERING\n");
		printf("   21: TRUEANOMALY");
		printf("     22: ECCANOMALY");
		printf("     23: MEANANOMALY");
		printf("     24: PERIPASSAGE\n");
		scanf("%d", &i);
		switch (i)
			{
			case 	11 :  transkep( &comet , CARTESIAN );
				     break;
			case 	12 :  transkep( &comet , SPHERICAL );
				     break;
			case 	21 :  transkep( &comet , TRUEANOMALY );
				     break;
			case 	22 :  transkep( &comet , ECCANOMALY );
				     break;
			case 	23 :  transkep( &comet , MEANANOMALY );
				     break;
			case 	24 :  transkep( &comet , PERIPASSAGE );
				     break;
			case 	31 :  transkep( &comet , SCATTERING );
				     break;
			case 	32 :  transkep( &comet , TSCATTERING );
				     break;
			default:
			    error("main(): type %d unknown kepler orbit\n", i);
			}
		}
	}
	
/*-----------------------------------------------------------------------------
 *  kepinit  --  initializes the six phase space coordinates
 *		 for the kepler orbit and sets time t = 0.0
 *-----------------------------------------------------------------------------
 */
local void  kepinit(body)
KEPLERPTR  body;
	{
	long  seed;		/* seed for the random number generator	*/
	long  time();
        real  randunit();

	if ( KType(body) == TRUEANOMALY || KType(body) == ECCANOMALY ||
	     KType(body) == MEANANOMALY || KType(body) == PERIPASSAGE)
		{
		printf("seed = ?  ( 0 < int < 2e9;");
		printf("  0 invokes a random seed )\n");
		scanf("%ld", &seed);
		if (seed == 0L)
			{
/*
 * uncomment the following line on a UNIX machine:
 */
			seed = time(0);
/*
 * uncomment the following line on a non-UNIX machine:
 *			error("no UNIX clock available\n");
 */
			}
		srandunit(seed);
		}
	KTime(body) = 0.0;
	printf("KMass1 = ?\n");
	scanf("%lf", &KMass1(body) );
	if ( KMass1(body) < 0.0 )
		error("mass1 is negative\n");
	printf("KMass2 = ?\n");
	scanf("%lf", &KMass2(body) );
	if ( KMass2(body) < 0.0 )
		error("mass2 is negative\n");

	switch (KType(body))
		{
		case	 CARTESIAN :
				printf("z = ?\t( z component of position )");
				printf("  [-infinity, +infinity] )\n");
				scanf("%lf", &KRvect(body)[0] );
				break;
		case	  SPHERICAL :
				printf("r = ?\t( radius vector: a scalar )");
				printf("  [0, +infinity] )\n");
				scanf("%lf", &KRmod(body) );
				if (KRmod(body) < 0.0)
					error("no negative radial distance\n");
				break;
		case	TRUEANOMALY :
		case	 ECCANOMALY :
		case	MEANANOMALY :
		case	PERIPASSAGE :
				printf("a = ?\t( semimajor axis )");
				printf("  [0, +infinity] )\n");
				scanf("%lf", &KAsemi(body) );
				if (KAsemi(body) < 0.0)
					error("no negative semimajor axis\n");
				break;
		case     SCATTERING :
				printf("rho = ?\t( impact parameter )");
				printf("  (0, +infinity) \n");
				scanf("%lf", &KImp_param(body) );
				if (KImp_param(body) < 0.0)
				       error("no negative impact parameter\n");
				break;
		case     TSCATTERING :
				printf("rho = ?\t( impact parameter )");
				printf("  (0, +infinity) \n");
				scanf("%lf", &KImp_param(body) );
				if (KImp_param(body) < 0.0)
				       error("no negative impact parameter\n");
				break;
		default	:  error("kepinit(): unknown type of kepler orbit\n");
		}

	if ( KType(body) == TRUEANOMALY || KType(body) == ECCANOMALY ||
						KType(body) == MEANANOMALY )
		{
		printf("A random number will be substituted");
		printf(" whenever you type in a negative number.\n");
		printf("This random number is drawn from a thermal");
		printf(" distribution of orbital elements,\n");
		printf("i.e. with equal probabilities for equal");
		printf(" phase space volumes\n");
		}


	switch (KType(body))
		{
		case	 CARTESIAN :
				printf("x = ?\t( x component of position )");
				printf("  [-infinity, +infinity] )\n");
				scanf("%lf", &KSpace(body)[1] );
				break;
		case	  SPHERICAL :
				printf("theta = ?\t( angle with z axis )");
				printf("  [0, pi] )\n");
				scanf("%lf", &KTheta(body) );
				if (KTheta(body) > PI)
					error("too large a theta\n");
				if (KTheta(body) < 0.0)
					error("no negative theta\n");
				break;
		case	TRUEANOMALY :
		case	 ECCANOMALY :
		case	MEANANOMALY :
		case	PERIPASSAGE :
				printf("e = ?\t( eccentricity;");
				printf("  [0, infinity] )\n");
				scanf("%lf", &KEcc(body) );
				if (KEcc(body) < 0.0)
					KEcc(body) = sqrt( randunit() );
				break;
		case     SCATTERING :
				printf("psi = ?\t( impact angle )");
				printf("  [0, 2pi] \n");
				scanf("%lf", &KImp_angle(body) );
				if (KImp_angle(body) > TWO_PI)
					error("too large an impact angle\n");
				if (KImp_angle(body) < 0.0)
					KImp_angle(body) = TWO_PI * randunit();
				break;
		case     TSCATTERING :
				printf("psi = ?\t( impact angle )");
				printf("  [0, 2pi] \n");
				scanf("%lf", &KImp_angle(body) );
				if (KImp_angle(body) > TWO_PI)
					error("too large an impact angle\n");
				if (KImp_angle(body) < 0.0)
					KImp_angle(body) = TWO_PI * randunit();
				break;
		default	:  error("kepinit(): unknown type of kepler orbit\n");
		}

	switch (KType(body))
		{
		case	 CARTESIAN :
				printf("y = ?\t( y component of position )");
				printf("  [-infinity, +infinity] )\n");
				scanf("%lf", &KSpace(body)[2] );
				break;
		case	  SPHERICAL :
				printf("phi = ?\t( angle with x axis in the");
				printf(" projection onto the x,y plane )");
				printf("  [0, 2pi] )\n");
				scanf("%lf", &KPhi(body) );
				if (KPhi(body) > TWO_PI)
					error("too large a phi\n");
				if (KPhi(body) < 0.0)
					error("no negative phi\n");
				break;
		case	TRUEANOMALY :
		case	 ECCANOMALY :
		case	MEANANOMALY :
		case	PERIPASSAGE :
				printf("i = ?\t( inclination;  [0, pi] )\n");
				scanf("%lf", &KIncl(body) );
				if (KIncl(body) > PI)
					error("too large an inclination\n");
				if (KIncl(body) < 0.0)
					KIncl(body) = acos( randunit() );
				break;
		case     SCATTERING :
				printf("d = ?\t( initial distance )");
				printf("  (0, +infinity) \n");
				scanf("%lf", &KImp_separ(body) );
				if (KImp_separ(body) < 0.0)
				       error("no negative initial distance\n");
				break;
		case     TSCATTERING :
				printf("t = ?\t( encounter time (< 0.0)");
				printf("  (-infinity,0) \n");
				scanf("%lf", &KImp_time(body) );
				break;
		default	:  error("kepinit(): unknown type of kepler orbit\n");
		}

	switch (KType(body))
		{
		case	 CARTESIAN :
				printf("vz = ?\t( z component of velocity )");
				printf("  [-infinity, +infinity] )\n");
				scanf("%lf", &KSpace(body)[3] );
				break;
		case	  SPHERICAL :
				printf("v_rad = ?\t( radial velocity )");
				printf("  [-infinity, +infinity] )\n");
				scanf("%lf", &KSpace(body)[3] );
				break;
		case	TRUEANOMALY :
		case	 ECCANOMALY :
		case	MEANANOMALY :
		case	PERIPASSAGE :
				printf("OMEGA = ?\t( longitude of the");
				printf(" ascending node;  [0, 2pi] )\n");
				scanf("%lf", &KLonasc(body) );
				if (KLonasc(body) > TWO_PI)
				      error("too large a long. of asc. node\n");
				if (KLonasc(body) < 0.0)
					KLonasc(body) = TWO_PI * randunit();
				break;
		case     SCATTERING :
				printf("Vinf = ?\t( velocity at infinity )");
				printf("  (0, +infinity) \n");
				scanf("%lf", &KVinf_abs(body) );
				if (KVinf_abs(body) < 0.0)
					warning("outgoing orbit\n");
				break;
		case     TSCATTERING :
				printf("Vinf = ?\t( velocity at infinity )");
				printf("  (0, +infinity) \n");
				scanf("%lf", &KVinf_abs(body) );
				if (KVinf_abs(body) < 0.0)
					error("no negative vel. at inf.\n");
				break;
		default	:  error("kepinit(): unknown type of kepler orbit\n");
		}

	switch (KType(body))
		{
		case	 CARTESIAN :
				printf("vx = ?\t( x component of velocity )");
				printf("  [-infinity, +infinity] )\n");
				scanf("%lf", &KSpace(body)[4] );
				break;
		case	  SPHERICAL :
				printf("v_theta = ?\t( velocity in the");
				printf(" direction of increasing theta )");
				printf("  [-infinity, +infinity] )\n");
				scanf("%lf", &KSpace(body)[4] );
				break;
		case	TRUEANOMALY :
		case	 ECCANOMALY :
		case	MEANANOMALY :
		case	PERIPASSAGE :
				printf("omega = ?\t( argument of pericenter;");
				printf("  [0, 2pi] )\n");
				scanf("%lf", &KArgper(body) );
				if (KArgper(body) > TWO_PI)
				       error("too large an arg. of perictr\n");
				if (KArgper(body) < 0.0)
					KArgper(body) = TWO_PI * randunit();
				break;
		case     SCATTERING :
				printf("theta = ?\t( incoming angle )");
				printf("  [0, pi] \n");
				scanf("%lf", &KVinf_theta(body) );
				if (KVinf_theta(body) > PI)
					error("too large a theta\n");
				if (KVinf_theta(body) < 0.0)
					KVinf_theta(body) = acos( -1.0
							   + 2.0 * randunit());
				break;
		case     TSCATTERING :
				printf("theta = ?\t( incoming angle )");
				printf("  [0, pi] \n");
				scanf("%lf", &KVinf_theta(body) );
				if (KVinf_theta(body) > PI)
					error("too large a theta\n");
				if (KVinf_theta(body) < 0.0)
					KVinf_theta(body) = acos( -1.0
							   + 2.0 * randunit());
				break;
		default	:  error("kepinit(): unknown type of kepler orbit\n");
		}

	switch (KType(body))
		{
		case	 CARTESIAN :
				printf("vy = ?\t( y component of velocity )");
				printf("  [-infinity, +infinity] )\n");
				scanf("%lf", &KSpace(body)[5] );
				break;
		case	  SPHERICAL :
				printf("v_phi = ?\t( velocity in the");
				printf(" direction of increasing phi )");
				printf("  [-infinity, +infinity] )\n");
				scanf("%lf", &KSpace(body)[5] );
				break;
		case	TRUEANOMALY :
				printf("f = ?\t( true anomaly;");
				printf("  [0, 2pi] )\n");
				scanf("%lf", &KTruean(body) );
				if (KTruean(body) >  /* hyperbola constraint */
				   (PI - KAsemi(body)*acos(1.0/KEcc(body)))
				   && KEcc(body) > 1.0)
					error("too large an true anomaly\n");
				if (KTruean(body) > TWO_PI)
					error("too large a true anomaly\n");
				if (KTruean(body) < 0.0 && KEcc(body) < 1.0)
					KTruean(body) = TWO_PI * randunit();
				if (KTruean(body) <  /* hyperbola constraint */
				    (PI - KAsemi(body)*acos(1.0/KEcc(body)))
				    && KEcc(body) > 1.0)
					error("too small a true anomaly\n");
				break;
		case	 ECCANOMALY :
				printf("E = ?\t( eccentric anomaly;");
				printf("  [0, 2pi] )\n");
				scanf("%lf", &KEccan(body) );
				if (KEccan(body) > TWO_PI && KEcc(body) < 1.0)
					error("too large an ecc. anomaly\n");
				if (KEccan(body) < 0.0 && KEcc(body) < 1.0)
					KEccan(body) = TWO_PI * randunit();
				break;
		case	MEANANOMALY :
				printf("M = ?\t( mean anomaly;");
				printf("  [0, 2pi] )\n");
				scanf("%lf", &KMeanan(body) );
				if (KMeanan(body) > TWO_PI && KEcc(body) < 1.0)
					error("too large a mean anomaly\n");
				if (KMeanan(body) < 0.0 && KEcc(body) < 1.0)
					KMeanan(body) = TWO_PI * randunit();
				break;
		case	PERIPASSAGE :
				printf("T = ?\t( time of pericenter passage");
				printf("  [-infinity, +infinity] )\n");
				printf("\n");
				scanf("%lf", &KPerpas(body) );
				break;

		case     SCATTERING :
				printf("phi = ?\t( incoming angle )");
				printf("  [0, 2pi] \n");
				scanf("%lf", &KVinf_phi(body) );
				if (KVinf_phi(body) > TWO_PI)
					error("too large a phi\n");
				if (KVinf_phi(body) < 0.0)
					KVinf_phi(body) = TWO_PI * randunit();
				break;
		case     TSCATTERING :
				printf("phi = ?\t( incoming angle )");
				printf("  [0, 2pi] \n");
				scanf("%lf", &KVinf_phi(body) );
				if (KVinf_phi(body) > TWO_PI)
					error("too large a phi\n");
				if (KVinf_phi(body) < 0.0)
					KVinf_phi(body) = TWO_PI * randunit();
				break;
		default	:  error("kepinit(): unknown type of kepler orbit\n");
		}
	}

/*-----------------------------------------------------------------------------
 *  randunit  --  returns a random real number within the unit interval
 *		  note: based on      @(#)rand.c   4.1 (Berkeley) 12/21/80,
 *			but returning a positive number smaller than unity.
 *  srandunit  --  accepts a seed to start up randunit
 *-----------------------------------------------------------------------------
 */
#define  MAXNUM	2147483647.0 	/* the maximum value which rand() can return */

static	long	randx = 1;

local void  srandunit(x)
long x;
	{
	randx = x;
	}

local real  randunit()
    {
    return((real)((randx= randx * 1103515245 + 12345) & 0x7fffffff)/MAXNUM);
    }

#undef  MAXNUM

/*-----------------------------------------------------------------------------
 *  kepprint  --  prints the six components of a kepler orbit   *body
 *	          note: different type of representations are invoked by
 *			the same function call, the proper coordinate
 *			being chosen using the information from KType(body)
 *-----------------------------------------------------------------------------
 */
local void  kepprint( body )
KEPLERPTR  body;
	{
        if (NDIM != 3)
	        error("transkep(): NDIM != 3\n");

	switch (KType(body))
		{
		case	CARTESIAN    :	print_cartesian(body); break;
		case	SPHERICAL    :	print_spherical(body); break;
		case	TRUEANOMALY  :	print_trueanomaly(body); break;
		case	ECCANOMALY   :	print_eccanomaly(body); break;
		case	MEANANOMALY  :	print_meananomaly(body); break;
		case	PERIPASSAGE  :	print_peripassage(body); break;
		case	SCATTERING   :	print_scattering(body); break;
		case	TSCATTERING  :	print_tscattering(body); break;
		default	:  error("kepprint(): kepler orbit type unknown\n");
		}
	}

/*-----------------------------------------------------------------------------
 *  print_cartesian  --  print a kepler orbit with coordinate type CARTESIAN
 *-----------------------------------------------------------------------------
 */
local void  print_cartesian(body)
KEPLERPTR  body;
	{
	printf(" x = %f    y = %f    z = %f\n",
			KRvect(body)[1], KRvect(body)[2], KRvect(body)[0]);
	printf("vx = %f   vy = %f   vz = %f\n",
			KVvect(body)[1], KVvect(body)[2], KVvect(body)[0]);
	}

/*-----------------------------------------------------------------------------
 *  print_spherical  --  print a kepler orbit of coordinate type SPHERICAL
 *-----------------------------------------------------------------------------
 */
local void  print_spherical(body)
KEPLERPTR  body;
	{
	printf("rmod = %f    theta = %f    phi = %f\n",
					KRmod(body), KTheta(body), KPhi(body));
	printf("vrad = %f   vtheta = %f   vphi = %f\n",
				      KVrad(body), KVtheta(body), KVphi(body));
	}

/*-----------------------------------------------------------------------------
 *  print_trueanomaly  --  print a kepler orbit of coordinate type TRUEANOMALY
 *-----------------------------------------------------------------------------
 */
local void  print_trueanomaly(body)
KEPLERPTR  body;
	{
	printf("      a = %f         e = %f   i = %f\n",
					KAsemi(body), KEcc(body), KIncl(body));
	printf("OMEGA = %f     omega = %f     f = %f\n",
				  KLonasc(body), KArgper(body), KTruean(body));
	}

/*-----------------------------------------------------------------------------
 *  print_eccanomaly  --  print a kepler orbit of coordinate type ECCANOMALY
 *-----------------------------------------------------------------------------
 */
local void  print_eccanomaly(body)
KEPLERPTR  body;
	{
	printf("      a = %f         e = %f        i = %f\n",
					KAsemi(body), KEcc(body), KIncl(body));
	printf("OMEGA = %f     omega = %f     E = %f\n",
				   KLonasc(body), KArgper(body), KEccan(body));
	}

/*-----------------------------------------------------------------------------
 *  print_meananomaly  --  print a kepler orbit of coordinate type MEANANOMALY
 *-----------------------------------------------------------------------------
 */
local void  print_meananomaly(body)	
KEPLERPTR  body;
	{
	printf("      a = %f         e = %f        i = %f\n",
					KAsemi(body), KEcc(body), KIncl(body));
	printf("OMEGA = %f     omega = %f     M = %f\n",
				  KLonasc(body), KArgper(body), KMeanan(body));
	}

/*-----------------------------------------------------------------------------
 *  print_peripassage  --  print a kepler orbit of coordinate type PERIPASSAGE
 *-----------------------------------------------------------------------------
 */
local void  print_peripassage(body)	
KEPLERPTR  body;
	{
	printf("      a = %f         e = %f        i = %f\n",
					KAsemi(body), KEcc(body), KIncl(body));
	printf("OMEGA = %f     omega = %f     T = %f\n",
				  KLonasc(body), KArgper(body), KPerpas(body));
	}

/*-----------------------------------------------------------------------------
 *  print_scattering  --  print a kepler orbit of coordinate type SCATTERING
 *-----------------------------------------------------------------------------
 */
local void  print_scattering(body)
KEPLERPTR  body;
	{
	printf("imp.param. = %f    imp.angle = %f   init.separ. = %f\n",
			 KImp_param(body), KImp_angle(body), KImp_separ(body));
	printf(" v at inf. = %f        theta = %f           phi = %f\n",
			  KVinf_abs(body), KVinf_theta(body), KVinf_phi(body));
	}

/*-----------------------------------------------------------------------------
 *  print_tscattering  --  print a kepler orbit of coordinate type TSCATTERING
 *-----------------------------------------------------------------------------
 */
local void  print_tscattering(body)
KEPLERPTR  body;
	{
	printf("imp.param. = %f    imp.angle = %f   encount.time = %f\n",
			 KImp_param(body), KImp_angle(body), KImp_time(body));
	printf(" v at inf. = %f        theta = %f           phi = %f\n",
			  KVinf_abs(body), KVinf_theta(body), KVinf_phi(body));
	}

#endif

/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*...........................................................................*/
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*===========================================================================*/
/*
 *    here is where the actual kepler transformations are done :
 */
/*===========================================================================*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/
/*...........................................................................*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/

/*  keptrans.c - ITERAC, MAXITER, asinh, acosh, can, can_eccanomaly,
		 can_meananomaly, can_peripassage, can_spherical,
                 can_tscattering, cast, cast_eccanomaly, 
		 cast_meananomaly, cast_peripassage, cast_spherical, 
		 cast_tscattering, classify, keplerseq, keplershypeq,
		 tran_cs_cartesian_to_scattering,
		 tran_ct_cartesian_to_trueanomaly,
		 tran_sc_scattering_to_cartesian,
		 tran_tc_trueanomaly_to_cartesian,
		 trans_to_cartesian, trans_to_scattering, trans_to_trueanomaly,
		 transclass, transkep, */

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 *
 *  keptrans.c:			(c) 1986  Piet Hut  Princeton, NJ, USA
 *
 *  This block of functions performs transformations between different types
 *  of phase space coordinate choices for the description of Kepler orbits.
 *  All transformations are invoked by a single call to the function
 *	transkep(body, type)
 *  where the information of the old Kepler orbit is contained in the
 *  structure   *body  , and the the new type is given by   type  .
 *
 *           $!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$
 *	     $								 $
 *	     $  BEWARE: 3-dim. Cartesian coordinates are numbered as	 $
 *	     $								 $
 *	     $	            { z, x, y }  =  { 0, 1, 2 }  .		 $
 *	     $								 $
 *	     $		Using other conventions is likely to lead to	 $
 *	     $		errors in invoking coordinate transformations!!  $
 *	     $								 $
 *	     $   BEWARE!!       BEWARE!!       BEWARE!!       BEWARE!!   $
 *	     $								 $
 *           $!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$
 *
 *:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 *
 *  Implementation:
 *
 *  The final transformation is a combination of only three types
 *  of transformation. Within each class of coordinate systems there
 *  is one prototype which has the same name as the class. To get from
 *  A to B, the first transformation is from A to prototype(A), then from
 *  prototype(A) to prototype(B), and finally from prototype(B) to B.
 *  these three transformations are called can, trans, cast, respectively.
 *
 *   ___________     trans     __________________________________________
 *  |           |------------>|                                          |
 *  | CARTESIAN |<------------|               TRUEANOMALY                |
 *  |___________|             |__________________________________________|
 *       |  A                  c |  A             |  A             |  A
 *       |  |                  a |  | c           |  |             |  |
 *       |  |                  s |  | a           |  |             |  |
 *       |  |                  t |  | n           |  |             |  |
 *   ____V__|___            _____V__|___     _____V__|____     ____V__|_____
 *  |           |          |            |   |             |   |             |
 *  | SPHERICAL |          | ECCANOMALY |   | MEANANOMALY |   | PERIPASSAGE |
 *  |___________|          |____________|   |_____________|   |_____________|
 *
 *:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 *
 *  BUGS:	  Singularities in the coordinate systems should be handled
 *		gracefully, and lead to reasonable transformed values,
 *		perhaps with a slight arbitrary error to circumvent having
 *		to divide by zero.
 *		  Near-parabolic orbits should really be treated separately,
 *		with an algorithm that does not produce singularities
 *		around eccentricity e = 1.
 *		  Home work!
 *
 *:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 */

/*
 * commented out here, since we have combined everything in one file:
 *
 *  #include "kepler.h"
 *
 */

/*-----------------------------------------------------------------------------
 *  transkep  --  transforms a kepler orbit to a different coordinate system
 *		  accepts: body: a pointer to a kepler orbit,
 *			   type: the type of the new coordinate system.
 *		  effect: *body will have its phasespace coordinates recast
 *			  into the new coordinate form after the function call.
 *-----------------------------------------------------------------------------
 */
local void  transkep(body, type)
KEPLERPTR  body;
int  type;
	{
	int  class;

        if (NDIM != 3)
	        error("transkep(): NDIM != 3\n");

	if (KType(body) == type)		/* no change 		*/
		return;
	can(body);			/* into a canonical type */
	classify( &class, type );
	if (KType(body) != class)	/* into a different class */
		transclass(body, class);
	cast(body, type);		/* into a particular type */
	}

/*-----------------------------------------------------------------------------
 *  classify -- determines the class of a type, given the type
 *-----------------------------------------------------------------------------
 */
local void  classify( class, type )
int  *class, type;
	{
	switch (type)
		{
		case	CARTESIAN   :
		case	SPHERICAL    :  *class = CARTESIAN; break;
		case	TRUEANOMALY  :
		case	ECCANOMALY   :
		case	MEANANOMALY  :
		case	PERIPASSAGE  :  *class = TRUEANOMALY; break;
                case    SCATTERING   :
                case    TSCATTERING  :  *class = SCATTERING; break;
		default	: error("classify(,): unknown type of kepler orbit\n");
		}
	}

/*-----------------------------------------------------------------------------
 *  can  --  transform a kepler orbit into the canonical type within its class
 *-----------------------------------------------------------------------------
 */
local void  can(body)
KEPLERPTR  body;
	{
	switch (KType(body))
		{
		case	CARTESIAN   :  break;
		case	SPHERICAL   :  can_spherical(body); break;
		case	TRUEANOMALY :  break;
		case	ECCANOMALY  :  can_eccanomaly(body); break;
		case	MEANANOMALY :  can_meananomaly(body); break;
		case	PERIPASSAGE :  can_peripassage(body); break;
		case	SCATTERING  :  break;
		case	TSCATTERING :  can_tscattering(body); break;
		default	:  error("can(): type of kepler orbit not valid\n\n");
		}
	}

/*-----------------------------------------------------------------------------
 *  can_spherical  --  transforms from spherical to cartesian coordinates
 *-----------------------------------------------------------------------------
 */
local void  can_spherical(body)
KEPLERPTR  body;
	{
	KEPLER  temp, *p;

	if (KType(body) == SPHERICAL)
		KType(body) = CARTESIAN;	
	else
		error("can_spherical(): wrong type of kepler orbit\n\n");

	p = &temp;
	scopy(p, body);

	if (KRmod(p) < 0.0)		/* inconsistent */
		error("can_spherical(): radius vector has negative length\n");

	if (KTheta(p) < 0.0)	       /*  unusual theta values are suspect: */
		error("can_spherical(): spherical coordinate theta < 0\n");
	if (KTheta(p) > PI)
		error("can_spherical(): spherical coordinate theta > pi\n");

		    /*  unusual phi values might arize from wrapping around: */
	if (KPhi(p) < 0.0)
		warning("can_spherical(): spherical coordinate phi < 0\n");
	if (KPhi(p) > TWO_PI)
		warning("can_spherical(): spherical coordinate phi > 2pi\n");
	KPhi(p) = pos_angle( KPhi(p) );

	KRvect(body)[0] = KRmod(p) * cos(KTheta(p));
	KRvect(body)[1] = KRmod(p) * sin(KTheta(p)) * cos(KPhi(p));
	KRvect(body)[2] = KRmod(p) * sin(KTheta(p)) * sin(KPhi(p));
	KVvect(body)[0] = KVrad(p) * cos(KTheta(p))
		        - KVtheta(p) * sin(KTheta(p));
	KVvect(body)[1] = KVrad(p) * sin(KTheta(p)) * cos(KPhi(p)) +
			  KVtheta(p) * cos(KTheta(p)) * cos(KPhi(p))
		        - KVphi(p) * sin(KPhi(p));
	KVvect(body)[2] = KVrad(p) * sin(KTheta(p)) * sin(KPhi(p)) +
			  KVtheta(p) * cos(KTheta(p)) * sin(KPhi(p))
			+ KVphi(p) * cos(KPhi(p));
	}

/*-----------------------------------------------------------------------------
 *  can_eccanomaly  --  transform from eccentric anomaly to true anomaly coord.
 *-----------------------------------------------------------------------------
 */
local void  can_eccanomaly(body)
KEPLERPTR  body;
	{
	real  ecc_an;		/*  eccentric anomaly (ellipse), or
				 *  the equivalent (hyperbola)
				 */
	if (KType(body) == ECCANOMALY)
		KType(body) = TRUEANOMALY;	
	else
		error("can_eccanomaly(): wrong type of kepler orbit\n");

	if (KEcc(body) < 0.0)		/* inconsistent */
		error("can_eccanomaly(): negative eccentricity\n");

	ecc_an = KEccan(body);
	if (KEcc(body) < 1.0)		/* ellipse */
		{
		if (ecc_an < 0.0)
			warning("can_eccanomaly(): eccentric anomaly < 0\n");
		if (ecc_an > TWO_PI)
			warning("can_eccanomaly(): eccentric anomaly > 2pi\n");
		ecc_an = pos_angle( ecc_an );
		KTruean(body) = 2.0 * atan( sqrt( (1.0+KEcc(body)) /
				      (1.0-KEcc(body)) ) * tan( ecc_an/2.0 ) );
		if ( ecc_an > PI )
		    KTruean(body) += TWO_PI;    /* to guarantee *true_an > 0 */
		}
	else			/* hyperbola: -infinity < ecc_an < +infinity */
		KTruean(body) = 2.0 * atan( sqrt( (KEcc(body)+1.0) /
				      (KEcc(body)-1.0) ) * tanh( ecc_an/2.0 ));
	}

/*-----------------------------------------------------------------------------
 *  can_meananomaly  --  transforms from mean anomaly to true anomaly coord.
 *-----------------------------------------------------------------------------
 */
local void  can_meananomaly(body)
KEPLERPTR  body;
	{
	real  mean_an;			/*  mean anomaly  */

	if (KType(body) == MEANANOMALY)
		KType(body) = TRUEANOMALY;	
	else
		error("can_meananomaly(): wrong type of kepler orbit\n");

	if (KEcc(body) < 0.0)		/* inconsistent */
		error("can_meananomaly(): negative eccentricity\n");

	mean_an = KMeanan(body);
	if (KEcc(body) < 1.0)		/* ellipse */
		{
		if (mean_an < 0.0)
			warning("can_meananomaly(): mean anomaly < 0\n");
		if (mean_an > TWO_PI)
			warning("can_meananomaly(): mean anomaly > 2pi\n");
		mean_an = pos_angle( mean_an );
		keplerseq( mean_an, KEcc(body), &KTruean(body) );	
		}
	else		       /* hyperbola: -infinity < mean_an < +infinity */
		keplershypeq( mean_an, KEcc(body), &KTruean(body) );
	}

/*-----------------------------------------------------------------------------
 *  can_peripassage  --  transforms from time-of-pericenter-passage to
 *			 true anomaly coordinates
 *-----------------------------------------------------------------------------
 */
local void  can_peripassage(body)
KEPLERPTR  body;
	{
	real  mtot;		/*  sum of the two masses		    */
	real  mean_an;		/*  mean anomaly			    */
	real  omega;		/*  orbit average of the angular velocity   */
			/*  these definitions are for an elliptic orbit;
			 *  for a hyperbola, the orbital period is infinite,
			 *  and Kepler's third law is used instead.
			 */

	if (KType(body) == PERIPASSAGE)
		KType(body) = TRUEANOMALY;	
	else
		error("can_peripassage(): wrong type of kepler orbit\n");

	if (KAsemi(body) < 0.0)		/* inconsistent */
		error("can_peripassage(): negative semimajor axis\n");

	if (KEcc(body) < 0.0)		/* inconsistent */
		error("can_peripassage(): negative eccentricity\n");

	mtot = KMass1(body) + KMass2(body);
		/*
		 *  Kepler's third law can yield a definition of
		 *  mean anomaly also in the hyperbolic case:
		 */
	omega = sqrt( mtot / (KAsemi(body)*KAsemi(body)*KAsemi(body)) );
	mean_an = omega * (KTime(body) - KPerpas(body));
	if (KEcc(body) < 1.0)		/* ellipse */
		{
		mean_an = pos_angle( mean_an );	/* 0 < mean_an < 2pi */
		keplerseq( mean_an, KEcc(body), &KTruean(body) );	
		}
	else				/* hyperbola */
		keplershypeq( mean_an, KEcc(body), &KTruean(body) );
	}

/*-----------------------------------------------------------------------------
 *  keplerseq  --  solves Kepler's equation by iteration, for elliptic orbits
 *		   litt: H. Goldstein (1980), Classical Mechanics, eq. (3-76);
 *  			 S.W. McCuskey (1963), Introduction to Celestial
 *					       Mechanics, Sec. 3-7.
 *		   method: true_an follows from ecc_an, which is determined
 *			   by inverting Kepler's equation:
 *			   mean_an = ecc_an - ecc * sin(ecc_an)
 *-----------------------------------------------------------------------------
 */
#define  ITERAC	   1.0e-10	/* iteration accuracy for Kepler's equation */
#define  MAXITER   100		/* maximum number of iterations		    */

local int  keplerseq( mean_an, ecc, true_an )
real  mean_an;			/*  mean anomaly			*/
real  ecc;			/*  eccentricity			*/
real  *true_an;			/*  true anomaly			*/
	{
	real  ecc_an;		/*  eccentric anomaly			*/
	real  delta_ecc_an;	/*  iterative increment in ecc_an	*/
	real  function;		/*  function = 0  solves kepler's eq.	*/
	real  derivative;	/*  d function / d ecc_an		*/ 
	int  i;

	if (ecc < 0.0)			/* inconsistent */
		error("keplerseq(,,): negative eccentricity\n");
	if (ecc > 1.0)			/* inconsistent */
		error("keplerseq(,,): eccentricity >1 in an elliptic orbit\n");

	if (mean_an < 0.0)
		warning("keplerseq(,,): mean anomaly < 0\n");
	if (mean_an > TWO_PI)
		warning("keplerseq(,,): mean anomaly > 2pi\n");
	mean_an = pos_angle( mean_an );	   /* 0 < mean_an < 2pi */
	ecc_an = mean_an;		/* first guess for ecc_an	*/
	i = 0;
	delta_ecc_an = 1.0;		/* just to start the while loop	*/
	while ( ABS( delta_ecc_an ) > ITERAC )
		{
		if (++i > MAXITER)		/* convergence too slow	*/
			error("keplerseq(,,): convergence too slow\n");
		function = -mean_an + ecc_an - ecc * sin(ecc_an);
				/* function = 0 solves Kepler's equation */
		derivative = 1.0 - ecc * cos(ecc_an);
				/*  d(function) / d(ecc_an)  */
		delta_ecc_an = -function / derivative;
				/* Newton's method to find roots of f()	*/
		if ( delta_ecc_an > 1.0 )
			delta_ecc_an = 1.0;
		else if ( delta_ecc_an < -1.0 )
			delta_ecc_an = -1.0;	/* avoids large jumps	*/
		ecc_an += delta_ecc_an;
		}
	*true_an = 2.0 * atan( sqrt((1.0+ecc)/(1.0-ecc)) * tan(ecc_an/2.0) );
	if ( ecc_an > PI )
		*true_an += TWO_PI;	/* to guarantee *true_an > 0	*/
	return(i);
	}
#undef   ITERAC
#undef   MAXITER

/*-----------------------------------------------------------------------------
 *  keplershypeq  --  solves Kepler's equation by iteration: hyperbolic orbits
 *		      litt: S.W. McCuskey (1963), Introduction to Celestial
 *					          Mechanics, Sec. 3-10.
 *		      method: true_an follows from ecc_an, which is determined
 *			      by inverting the hyperbolic analogy of
 *			      Kepler's equation:
 *			      mean_an = -ecc_an + ecc * sinh(ecc_an)
 *-----------------------------------------------------------------------------
 */
#define  ITERAC	   1.0e-10	/* iteration accuracy for Kepler's equation */
#define  MAXITER   100		/* maximum number of iterations		    */

#define	  asinh(x)   log( (x) + sqrt( (x)*(x) + 1 ) )  /* inverse of sinh() */
		/* 
		 * only correctly defined this way for positive arguments
		 * for negative arguments (not needed here) simply change
		 * sign of argument and result.
		 */

local int  keplershypeq( mean_an, ecc, true_an )
real  mean_an;			/*  mean anomaly			   */
real  ecc;			/*  eccentricity			   */
real  *true_an;			/*  true anomaly			   */
	{
	real  ecc_an;		/*  eccentric anomaly (hyperbolic analogy) */
	real  delta_ecc_an;	/*  iterative increment in ecc_an	   */
	real  function;		/*  function = 0  solves kepler's eq.	   */
	real  derivative;	/*  d function / d ecc_an		   */ 
 	int  i;

	if (ecc < 1.0)		/* inconsistent */
		if (ecc < 0.0)
			error("keplershypeq(,,): negative eccentricity\n");
		else
			error("keplershypeq(,,):e <1 in a hyperbolic orbit\n");
	if (mean_an > 0.0)
		ecc_an = asinh( mean_an / ecc );   /* first guess for ecc_an */
	else
		ecc_an = -asinh( -mean_an / ecc ); /* asinh > 0  at present  */

	i = 0;
	delta_ecc_an = 1.0;		/* just to start the while loop	*/
	while ( ABS( delta_ecc_an ) > ITERAC )
		{
		if (++i > MAXITER)		/* convergence too slow	*/
			error("keplershypeq(,,): convergence too slow\n");
		function = -mean_an - ecc_an + ecc * sinh(ecc_an);
				/* function = 0 solves Kepler's equation */
		derivative = -1.0 + ecc * cosh(ecc_an);
				/*  d(function) / d(ecc_an)  */
		delta_ecc_an = -function / derivative;
				/* Newton's method to find roots of f()	*/
		if  ( ABS(derivative) < 1.0 )	       /* avoids large jumps */
			{
			if ( delta_ecc_an > 1.0 )
				delta_ecc_an = 1.0;
			else if ( delta_ecc_an < -1.0 )
				delta_ecc_an = -1.0;
			}
		ecc_an += delta_ecc_an;
		}
	*true_an = 2.0 * atan( sqrt((ecc+1.0)/(ecc-1.0)) * tanh(ecc_an/2.0) );
	return(i);
	}

#undef   ITERAC
#undef   MAXITER
#undef	 asinh

/*-----------------------------------------------------------------------------
 *  can_tscattering --  transforms from tscattering to scattering coordinates,
 *			i.e.: changing from time-before-pericenter-passage to
 *                            initial-separation.
 *                      note: for outgoing orbits, this time is negative.
 *-----------------------------------------------------------------------------
 */
local void  can_tscattering(body)
KEPLERPTR  body;
	{
	real  mtot;		/*  sum of the two masses		    */
	real  true_an;		/*  true anomaly			    */
	real  mean_an;		/*  mean anomaly			    */
	real  omega;		/*  hyperbolic generalization,              */
                                /*    cf. can_peripassage()                 */
        real  semi_major_axis;
	real  eccentricity;
	real  a_cubed;
        real  ecc_sq;
	real  r_cross_v_sq;

	if (KType(body) == TSCATTERING)
		KType(body) = SCATTERING;
	else
		error("can_tscattering(): wrong type of kepler orbit\n");

	mtot = KMass1(body) + KMass2(body);
        semi_major_axis = mtot / (KVinf_abs(body) * KVinf_abs(body));
	a_cubed = semi_major_axis * semi_major_axis * semi_major_axis;
		/*
		 *  Kepler's third law can yield a definition of
		 *  mean anomaly also in the hyperbolic case:
		 */
	omega = sqrt( mtot / a_cubed );
	mean_an = omega * (KTime(body) - KImp_time(body));
        r_cross_v_sq = KImp_param(body) * KImp_param(body)
                      * KVinf_abs(body) * KVinf_abs(body);
        ecc_sq = 1.0 + r_cross_v_sq / (mtot * semi_major_axis);
        eccentricity = sqrt(ecc_sq);
	keplershypeq( mean_an, eccentricity, &true_an );

	KImp_separ(body) = semi_major_axis * (ecc_sq - 1.0)
	                  / (1.0 + eccentricity * cos(true_an));
	}

/*-----------------------------------------------------------------------------
 *  transclass  --  transforms from one canonical type to another in a
 *		    different class
 *-----------------------------------------------------------------------------
 */
local void  transclass(body, class)
KEPLERPTR  body;
int  class;
	{
	switch (class)
		{
		case	CARTESIAN   :  trans_to_cartesian(body); break;
		case	TRUEANOMALY :  trans_to_trueanomaly(body); break;
		case    SCATTERING  :  trans_to_scattering(body); break;
		default: error("transclass(,): kepler orbit type not valid\n");
		}
	}

/*-----------------------------------------------------------------------------
 *  trans_to_cartesian  --  transforms a kepler orbit to the canonical type of
 *			    cartesian coordinates, starting from another
 *			    canonical form in a different class
 *-----------------------------------------------------------------------------
 */
local void  trans_to_cartesian(body)
KEPLERPTR  body;
	{
	switch (KType(body))
		{
		case	TRUEANOMALY :  tran_tc_trueanomaly_to_cartesian(body);
				       return;
		case	SCATTERING  :  tran_sc_scattering_to_cartesian(body);
				       return;
		}
	error("trans_to_cartesian(): wrong type of kepler orbit\n");
	}

/*-----------------------------------------------------------------------------
 *  trans_to_trueanomaly  --  transforms a kepler orbit to the canonical type
 *			      of orbital elements (i.e. true anomaly), starting
 *			      from another canonical form in a different class
 *-----------------------------------------------------------------------------
 */
local void  trans_to_trueanomaly(body)
KEPLERPTR  body;
	{
	switch (KType(body))
		{
		case	SCATTERING :  tran_sc_scattering_to_cartesian(body);
						/* drop through to :   */
		case	CARTESIAN  :  tran_ct_cartesian_to_trueanomaly(body);
				      return;
		}
	error("trans_to_trueanomaly(): wrong type of kepler orbit\n");
	}

/*-----------------------------------------------------------------------------
 *  trans_to_scattering  --  transforms a kepler orbit to the canonical type of
 *			     scattering coordinates, starting from another
 *			     canonical form in a different class
 *-----------------------------------------------------------------------------
 */
local void  trans_to_scattering(body)
KEPLERPTR  body;
	{
	switch (KType(body))
		{
		case	TRUEANOMALY :  tran_tc_trueanomaly_to_cartesian(body);
						/* drop through to :   */
		case	CARTESIAN   :  tran_cs_cartesian_to_scattering(body);
				       return;
		}
	error("trans_to_scattering(): wrong type of kepler orbit\n");
	}

/*-----------------------------------------------------------------------------
 *  tran_ct_cartesian_to_trueanomaly  --  transforms from cartesian to
 *					  true anomaly coordinates
 *					  litt: S.W. McCuskey (1963),
 *						Introduction to Celestial
 *						Mechanics, Sec. 3-4.
 *					  note: the abbreviation _ct_ is added 
 *						to guarantee uniqueness of
 *						first eight characters.
 *-----------------------------------------------------------------------------
 */
local void  tran_ct_cartesian_to_trueanomaly( body )
KEPLERPTR  body;
	{
	KEPLER  temp, *p;
	real  mtot;		/*  sum of the two masses		*/
	real  rmod;		/*  absolute value of the radius vector	*/
	real  vsquare;		/*  square of the velocity		*/
	real  h, hh, c0, c1, c2;
			/*
			 *  h is the absolute value of r cross v, the
			 *  outer product of r and v, the components of
			 *  which are  c0, c1, c2;   hh = h*h.
			 */
	real  p0;		/*  semi latus rectum, i.e. the distance from
				 *  a focus to the conic in the direction
				 *  perpendicular to the major axis
				 */
	real  u0;		/*  angle from the ascending node to the 
				 *  position in the orbit
				 */

	if (KType(body) == CARTESIAN)
		KType(body) = TRUEANOMALY;	
	else
		error("tran_ct_cartesian_...(): wrong type of kepler orbit\n");

	p = &temp;
	scopy(p, body);

	mtot = KMass1(p) + KMass2(p);
	rmod = sqrt( KRvect(p)[0]*KRvect(p)[0] + KRvect(p)[1]*KRvect(p)[1] +
						 KRvect(p)[2]*KRvect(p)[2] );
	vsquare = KVvect(p)[0]*KVvect(p)[0] + KVvect(p)[1]*KVvect(p)[1] +
						 KVvect(p)[2]*KVvect(p)[2];
	KAsemi(body) = 1.0 / ( 2.0/rmod - vsquare/(mtot) );
			/*
			 *  to allow identical formula in both cases,
			 *  a trick is employed:
			 *  KAsemi(p) > 0 for an ellipse
			 *  KAsemi(p) < 0 for a hyperbola (restored at bottom)
			 *  parabolas are not yet considered
			 */
	c0 = (KRvect(p)[1]*KVvect(p)[2] - KRvect(p)[2]*KVvect(p)[1]);
	c1 = (KRvect(p)[2]*KVvect(p)[0] - KRvect(p)[0]*KVvect(p)[2]);
	c2 = (KRvect(p)[0]*KVvect(p)[1] - KRvect(p)[1]*KVvect(p)[0]);
	hh = c0*c0 + c1*c1 + c2*c2;
	h = sqrt( hh );
	KEcc(body) = sqrt( 1.0 - hh/(mtot*KAsemi(body)) );
	KIncl(body) = acos( c0/h );				/* [0,pi]   */
	KLonasc(body) = acos( -c2 / (h*sin(KIncl(body))) );	/* [0,2*pi) */
	if ( h*sin(KIncl(body))*sin(KLonasc(body))/c1 < 0.0 )
		KLonasc(body) = TWO_PI - KLonasc(body);
	p0 = hh / (mtot);
	KTruean(body) = acos( (p0/rmod - 1.0)/KEcc(body) );	/* [0,2*pi) */
	if ( KRvect(p)[0]*KVvect(p)[0] + KRvect(p)[1]*KVvect(p)[1]
					+ KRvect(p)[2]*KVvect(p)[2] < 0.0 )
		KTruean(body) = TWO_PI - KTruean(body);
	if (KEcc(body) > 1.0)		/* hyperbola */
		KTruean(body) = sym_angle( KTruean(body) );
	u0 = acos( (KRvect(p)[1]*cos(KLonasc(body))
		    + KRvect(p)[2]*sin(KLonasc(body)))/rmod );	/* [0,2*pi) */
	if ( (-KRvect(p)[1]*sin(KLonasc(body))
			  + KRvect(p)[2]*cos(KLonasc(body))) * cos(KIncl(body))
			  + KRvect(p)[0]*sin(KIncl(body)) < 0.0 )
		u0 = TWO_PI - u0;
	KArgper(body) = u0 - KTruean(body);			/* [0,2*pi) */
	if ( KArgper(body) < 0.0 )
		KArgper(body) += TWO_PI;
	KAsemi(body) = ABS( KAsemi(body) );	/*  also > 0 for hyperbolae */
	}

/*-----------------------------------------------------------------------------
 *  tran_tc_trueanomaly_to_cartesian  --  transforms from true anomaly to
 *					  cartesian coordinates.
 *					  litt: S.W. McCuskey (1963),
 *						Introduction to Celestial
 *						Mechanics, Sec. 3-12.
 *					  note: the abbreviation _tc_ is added 
 *						to guarantee uniqueness of
 *						first eight characters.
 *-----------------------------------------------------------------------------
 */
local void  tran_tc_trueanomaly_to_cartesian( body )
KEPLERPTR  body;
	{
	KEPLER  temp, *p;
	real  maxtruean;	/*  max. value for hyperbolic true anomaly */
	real  mtot;		/*  sum of the two masses		   */
	real  rmod;		/*  absolute value of the radius vector	   */
	real  h;		/*  the absolute value of r cross v,
				 *  the outer product of r and v.
				 */
	real  u0;		/*  angle from the ascending node to the 
				 *  position in the orbit
				 */
	real  drmoddf, dfdtoverr, dfdttimesr;		/*  see below  */

	if (KType(body) == TRUEANOMALY)
		KType(body) = CARTESIAN;	
	else
		error("tran_tc_trueanomaly..(): wrong type of kepler orbit\n");

	p = &temp;
	scopy(p, body);

	if (KAsemi(p) < 0.0)		/* inconsistent */
		error("tran_tc_trueanomaly...(): negative semimajor axis\n");

	if (KIncl(p) < 0.0)		/* inconsistent */
		error("tran_tc_trueanomaly...(): negative inclination\n");

	if (KIncl(p) > PI)		/* inconsistent */
		error("tran_tc_trueanomaly...(): inclination > pi\n");

	if (KLonasc(p) < 0.0)
		warning("tran_tc_trueanomaly...(): long. of asc. node < 0\n");
	if (KLonasc(p) > TWO_PI)
		warning("tran_tc_trueanomaly..(): long. of asc. node > 2pi\n");
	KLonasc(p) = pos_angle( KLonasc(p) );   /* 0 < KLonasc(p) < 2pi */

	if (KArgper(p) < 0.0)
		warning("tran_tc_trueanomaly...(): arg. of pericenter < 0\n");
	if (KArgper(p) > TWO_PI)
		warning("tran_tc_trueanomaly..(): arg. of pericenter > 2pi\n");
	KArgper(p) = pos_angle( KArgper(p) );   /* 0 < KArgper(p) < 2pi */

	if (KEcc(p) < 0.0)		/* inconsistent */
		error("tran_tc_trueanomaly...(): negative eccentricity\n");

	if (KEcc(p) < 1.0)		/* ellipse */
		{
		if (KTruean(p) < 0.0)
			warning("tran_tc_trueanomaly..(): true anomaly < 0\n");
		if (KTruean(p) > TWO_PI)
		       warning("tran_tc_trueanomaly..(): true anomaly >2pi\n");
		KTruean(p) = pos_angle( KTruean(p) );
		}

	if (KEcc(p) > 1.0)		/* hyperbola */
		{
		maxtruean = PI - acos( 1.0 / KEcc(p) );
				/*
				 *  the asymptotic value of the true anomaly
				 *  for very large separation has absolute
				 *  value   maxtruean  .
				 */
		if (KTruean(p) < -maxtruean)
		      error("tran_tc_trueanomaly..(): hyp. true an. < min.\n");
		if (KTruean(p) > maxtruean)
		      error("tran_tc_trueanomaly..(): hyp. true an. > max.\n");
		}

	mtot = KMass1(p) + KMass2(p);
	u0 = KTruean(p) + KArgper(p);
	u0 = pos_angle( u0 );	/*   0 < u0 < 2pi ;  not really      */
					/*  necessary, but more consistent.  */
			/*
			 *  to allow identical formula in both cases,
			 *  hereafter a trick is employed:
			 *  KAsemi(p) > 0 for an ellipse
			 *  KAsemi(p) < 0 for a hyperbola (restored at bottom)
			 *  parabolas are not yet considered
			 */
	if (KEcc(p) > 1.0)
		KAsemi(p) *= -1.0;
	rmod = KAsemi(p)*(1.0 - KEcc(p)*KEcc(p))
					     / (1.0 + KEcc(p)*cos(KTruean(p)));

			/*
			 *  calculating the radius vector is straightforward;
			 *  cf. McCuskey, eq. (3-93):
			 */
	KRvect(body)[0] = rmod * sin(u0) * sin(KIncl(p));
	KRvect(body)[1] = rmod * ( cos(u0)*cos(KLonasc(p))
			    - sin(u0)*sin(KLonasc(p))*cos(KIncl(p)) );
	KRvect(body)[2] = rmod * ( cos(u0)*sin(KLonasc(p))
			    + sin(u0)*cos(KLonasc(p))*cos(KIncl(p)) );
			/*
			 *  calculating the velocity vector involves two steps:
			 *
			 *  _v = d_r/dt = (d_r/df) * (df/dt) =
			 *  (d_r/df) * [ h/(rmod*rmod) ].
			 *
			 *  d_r/df = (d rmod/df) * _r / rmod +
			 *	     + rmod * (d [_r / rmod] /df).
			 *  The first term at the right hand side is
			 *  the radial velocity component:
			 */
	drmoddf = KAsemi(p) * KEcc(p) * ( 1.0 - KEcc(p)*KEcc(p) )
			    * sin (KTruean(p))
			    /( (1.0 + KEcc(p)*cos(KTruean(p)))
			    * (1.0 + KEcc(p)*cos(KTruean(p))) );
	h = sqrt( mtot*KAsemi(p)*(1.0 - KEcc(p)*KEcc(p)) );
	dfdtoverr = h /(rmod*rmod*rmod);
	KVvect(body)[0] = drmoddf * KRvect(body)[0] * dfdtoverr;
	KVvect(body)[1] = drmoddf * KRvect(body)[1] * dfdtoverr;
	KVvect(body)[2] = drmoddf * KRvect(body)[2] * dfdtoverr;
			/*
			 *  The second term is the velocity component
			 *  perpendicular to the radius vector:
			 */
	dfdttimesr = h/rmod;
	KVvect(body)[0] += dfdttimesr * cos(u0) * sin(KIncl(p));
	KVvect(body)[1] += dfdttimesr * ( -sin(u0)*cos(KLonasc(p))
			     - cos(u0)*sin(KLonasc(p))*cos(KIncl(p)) );
	KVvect(body)[2] += dfdttimesr * ( -sin(u0)*sin(KLonasc(p))
			     + cos(u0)*cos(KLonasc(p))*cos(KIncl(p)) );
	}

/*-----------------------------------------------------------------------------
 *  tran_cs_cartesian_to_scattering  --  transforms from cartesian to
 *					 true anomaly coordinates
 *					 litt: P. Hut and J.N. Bahcall (1983),
 *					       Astrophys. J. 268, p. 319.
 *					 restriction: only valid for
 *						      hyperbolic orbits!
 *					 note: the abbreviation _ct_ is added 
 *					       to guarantee uniqueness of
 *					       first eight characters.
 *-----------------------------------------------------------------------------
 */
local void  tran_cs_cartesian_to_scattering( body )
KEPLERPTR  body;
	{
	KEPLER  temp1, *p1;
	KEPLER  temp2, *p2;
	real  mtot;		/*  sum of the two masses		   */
	int  inorout;		/*  +1 for ingoing orbits; -1 for outgoing */
	real  rho_aux[3];	/*  auxiliary vector parallel to rho,      */
				/*    the impact vector			   */
	real  vinf_aux[3];	/*  auxiliary vector parallel to vinf,     */
				/*    the velocity vector at infinity      */
	real  semi_deflect;	/*  half of the total angle of deflection, */
				/*    i.e. half the angle between the      */
				/*    asymptotes of the hyperbolic orbit   */
	real  normalization;	/*  temporary variable used for 	   */
				/*    normalization of vectors		   */
	real  unitv[3];		/*  unit vector in the plane spanned by    */
				/*    the z-axis and the velocity vector   */
				/*    at infinity, perpendicular to the    */
				/*    latter; the impact angle is the      */
				/*    angle from this unit vector to the   */
				/*    impact vector, in counterclockwise   */
				/*    direction as seen from the origin    */

	if (KType(body) == CARTESIAN)
		KType(body) = SCATTERING;
	else
		error("tran_cs_cartesian_...(): wrong type of kepler orbit\n");

	mtot = KMass1(body) + KMass2(body);

	p1 = &temp1;
	p2 = &temp2;
	scopy(p1, body);
	KImp_separ(body) = sqrt( KRvect(p1)[0] * KRvect(p1)[0] +
			 	KRvect(p1)[1] * KRvect(p1)[1] +
			 	KRvect(p1)[2] * KRvect(p1)[2] );
	KType(p1) = CARTESIAN;   /* in order to allow the next function call */
	tran_ct_cartesian_to_trueanomaly( p1 );
	if (KEcc(p1) < 1.0)		/* ellipse */
		error("tran_cs_cartesian_...(): elliptic orbit\n");
	if (KTruean(p1) < 0.0)
		inorout = 1.0;		/* ingoing orbit */
	else
		inorout = -1.0;		/* outgoing orbit */

	semi_deflect = asin( 1.0 / KEcc(p1) );
	KVinf_abs(body) = sqrt( mtot/KAsemi(p1) ) * inorout;
	scopy(p2, p1);
	KTruean(p2) = (HALF_PI - semi_deflect) * inorout;
	tran_tc_trueanomaly_to_cartesian( p2 );
	vinf_aux[0] = KRvect(p2)[0];
	vinf_aux[1] = KRvect(p2)[1];
	vinf_aux[2] = KRvect(p2)[2];
	normalization = sqrt(vinf_aux[0] * vinf_aux[0] +
			     vinf_aux[1] * vinf_aux[1] +
			     vinf_aux[2] * vinf_aux[2]);
	vinf_aux[0] /= normalization;
	vinf_aux[1] /= normalization;
	vinf_aux[2] /= normalization;
	          /*
                   * in the following four lines minus signs have
                   * been added to  vinf_aux  because the theta and phi
                   * angles point in the direction of the incoming star,
                   * and thus in the opposite direction as  vinf_aux
                   */
	KVinf_theta(body) = acos( -vinf_aux[0] );
	KVinf_phi(body) = atan( vinf_aux[2] / vinf_aux[1] );
	if (-vinf_aux[1] < 0.0)
		KVinf_phi(body) += PI;
	scopy(p2, p1);
	KTruean(p2) = - semi_deflect * inorout;
	tran_tc_trueanomaly_to_cartesian( p2 );
	rho_aux[0] = KRvect(p2)[0];
	rho_aux[1] = KRvect(p2)[1];
	rho_aux[2] = KRvect(p2)[2];
	normalization = sqrt(rho_aux[0] * rho_aux[0] +
			     rho_aux[1] * rho_aux[1] +
			     rho_aux[2] * rho_aux[2] );
	rho_aux[0] /= normalization;
	rho_aux[1] /= normalization;
	rho_aux[2] /= normalization;
	unitv[0] = sin( KVinf_theta( body ) );
	unitv[1] = - cos( KVinf_theta( body ) ) * cos( KVinf_phi( body ) );
	unitv[2] = - cos( KVinf_theta( body ) ) * sin( KVinf_phi( body ) );
	KImp_param( body ) = KAsemi(p1) * sqrt( KEcc(p1)*KEcc(p1) - 1.0 );
	KImp_angle( body ) = acos( unitv[0] * rho_aux[0] +
				  unitv[1] * rho_aux[1] +
				  unitv[2] * rho_aux[2] );
	if (vinf_aux[0] * (unitv[1]*rho_aux[2] - unitv[2]*rho_aux[1]) +
	    vinf_aux[1] * (unitv[2]*rho_aux[0] - unitv[0]*rho_aux[2]) +
	    vinf_aux[2] * (unitv[0]*rho_aux[1] - unitv[1]*rho_aux[0]) < 0.0)
		KImp_angle( body ) *= -1.0;
	}

/*-----------------------------------------------------------------------------
 *  tran_sc_scattering_to_cartesian  --  transforms from cartesian to
 *					 true anomaly coordinates
 *					 litt: P. Hut and J.N. Bahcall (1983),
 *					       Astrophys. J. 268, p. 319.
 *					 restriction: only valid for
 *						      hyperbolic orbits!
 *					 note: the abbreviation _ct_ is added 
 *					       to guarantee uniqueness of
 *					       first eight characters.
 *-----------------------------------------------------------------------------
 */
local void  tran_sc_scattering_to_cartesian( body )
KEPLERPTR  body;
	{
	KEPLER  temp, *p;
	real  mtot;		/*  sum of the two masses		 */
	real  rin;		/*  initial separation of the two bodies */
	real  rho;		/*  impact parameter			 */
	real  psi;		/*  impact angle, measured in		 */
				/*  counterclockwise direction 		 */
				/*  as seen from the origin		 */
	real  theta, phi;	/*  direction of asymptotic motion of    */
				/*  hyperbolic orbit (ingoing if	 */
				/*  KVinf_abs(body) > 0.0; outgoing if	 */
				/*  KVinf_abs(body) < 0.0)		 */
	real  ain;		/*  semimajor axis of hyperbolic orbit   */
	real  ein;		/*  eccentricity of hyperbolic orbit     */
	real  semi_deflect;	/*  half the total angle of deflection,  */
				/*  i.e. half the angle between the      */
				/*  asymptotes of the hyperbolic orbit   */
	real  meananin;		/*  initial mean anomaly 		 */
	real  phiin;		/*  initial angle in orbital plane	 */
				/*  between ingoing (outgoing) body	 */
				/*  and normal to ingoing (outgoing)	 */
				/*  asymptote				 */
	real  vin;		/*  initial velocity in absolute value	 */
				/*    (at initial position, therefore	 */
				/*     larger than KVinf_abs(body) ! )	 */
	real  xin, yin;		/*  initial positions projected on the   */
				/*  x and y axes in the orbital plane;	 */
				/*    the y axis is parallel to the	 */
				/*    incoming (outgoing) asymptote	 */
				/*    if the body is ingoing (outgoing)	 */
	real  vxin, vyin;	/*  initial velocities projected on the	 */
				/*  x and y axes in the orbital plane	 */
	int  inorout;		/*  +1 for ingoing orbits;		 */
				/*  -1 for outgoing orbits		 */
	real  sp, cp;		/*  dummy variables			 */


	if (KType(body) == SCATTERING)
		KType(body) = CARTESIAN;	
	else
	       error("tran_sc_scattering_...(): wrong type of kepler orbit\n");

	mtot = KMass1(body) + KMass2(body);
	rin = KImp_separ(body);
	rho = KImp_param(body);
	psi = KImp_angle(body);
	psi = PI - psi;  /* to conform to 1981 notation */
	theta = KVinf_theta(body);
	phi = KVinf_phi(body);
	inorout = (KVinf_abs(body) > 0.0)  ?  1  :  -1;
	ain = mtot / (KVinf_abs(body)*KVinf_abs(body));
	ein = sqrt( 1.0 + (rho*rho) / (ain*ain) );
	semi_deflect = asin( 1.0 / ein );
	meananin = acos( ( (rho*rho) / (ain*rin) - 1.0 ) / ein );
	phiin = meananin - semi_deflect;
	sp = sin( phiin );
	cp = cos( phiin );
	xin = cp * rin;
	yin = sp * rin;
	vin = sqrt( (KVinf_abs(body)*KVinf_abs(body)) + 2.0*mtot/rin );
	vxin = inorout * vin * (1.0 - sp) /
			 sqrt( (1.0-sp)*(1.0-sp) + (rho/ain+cp)*(rho/ain+cp) );
	vyin = inorout * sqrt( vin*vin - vxin*vxin );
	KRvect(body)[0] = yin * cos(theta) - xin * sin(theta) * cos(psi);
	KRvect(body)[1] = yin * sin(theta) * cos(phi) + xin *
		      (cos(theta) * cos(phi) * cos(psi) - sin(phi) * sin(psi));
	KRvect(body)[2] = yin * sin(theta) * sin(phi) + xin *
		      (cos(theta) * sin(phi) * cos(psi) + cos(phi) * sin(psi));
	KVvect(body)[0] = -vyin * cos(theta) + vxin * sin(theta) * cos(psi);
	KVvect(body)[1] = -vyin * sin(theta) * cos(phi) - vxin *
		      (cos(theta) * cos(phi) * cos(psi) - sin(phi) * sin(psi));
	KVvect(body)[2] = -vyin * sin(theta) * sin(phi) - vxin *
		      (cos(theta) * sin(phi) * cos(psi) + cos(phi) * sin(psi));
	}

/*-----------------------------------------------------------------------------
 *  cast  --  transform from the canonical type to a particular type
 *	      within the same class
 *-----------------------------------------------------------------------------
 */
local void  cast(body, type)
KEPLERPTR  body;
int  type;
	{
	switch (type)
		{
		case	CARTESIAN  :  break;
		case	SPHERICAL   :  cast_spherical(body); break;
		case	TRUEANOMALY :  break;
		case	ECCANOMALY  :  cast_eccanomaly(body); break;
		case	MEANANOMALY :  cast_meananomaly(body); break;
		case	PERIPASSAGE :  cast_peripassage(body); break;
		case	TSCATTERING :  cast_tscattering(body); break;
		}
	if (KType(body) != type)
		error("cast(,) :  something seriously wrong\n");
	}

/*-----------------------------------------------------------------------------
 *  cast_spherical  --  transform from cartesian to spherical coordinates
 *-----------------------------------------------------------------------------
 */
local void  cast_spherical(body)
KEPLERPTR  body;
	{
	KEPLER  temp, *p;

	if (KType(body) == CARTESIAN)
		KType(body) = SPHERICAL;
	else
		error("cast_spherical(): wrong type of kepler orbit\n");

	p = &temp;
	scopy(p, body);

	KRmod(body) = sqrt(KRvect(p)[0]*KRvect(p)[0]
				+ KRvect(p)[1]*KRvect(p)[1]
				+ KRvect(p)[2]*KRvect(p)[2]);
	KTheta(body) = atan( sqrt( KRvect(p)[1]*KRvect(p)[1]
				+ KRvect(p)[2]*KRvect(p)[2] ) / KRvect(p)[0] );
	if (KTheta(body) < 0.0)
		KTheta(body) += PI;
	KPhi(body) = atan( KRvect(p)[2] / KRvect(p)[1] );
	if (KPhi(body) < 0.0) 
		KPhi(body) += PI;
	if (KRvect(p)[2] < 0.0)
		KPhi(body) += PI;
	KVrad(body) = KVvect(p)[0]*cos(KTheta(body)) + sin(KTheta(body)) *
	       ( KVvect(p)[1]*cos(KPhi(body)) + KVvect(p)[2]*sin(KPhi(body)) );
	KVtheta(body) = -KVvect(p)[0]*sin(KTheta(body)) + cos(KTheta(body)) *
	       ( KVvect(p)[1]*cos(KPhi(body)) + KVvect(p)[2]*sin(KPhi(body)) );
	KVphi(body) = -KVvect(p)[1]*sin(KPhi(body)) +
		 KVvect(p)[2]*cos(KPhi(body));
	}

/*-----------------------------------------------------------------------------
 *  cast_eccanomaly  --  transform from true anomaly to eccentric anomaly 
 *-----------------------------------------------------------------------------
 */
#define	  acosh(x)   log( (x) + sqrt( (x)*(x) - 1 ) )  /* inverse of cosh() */

local void  cast_eccanomaly( body )
KEPLERPTR  body;
	{
	real  maxtruean;	/*  max. value for hyperbolic true anomaly */
	real  ecc_an;		/*  eccentric anomaly (ellipse), or
				 *  the equivalent (hyperbola)
				 */

	if (KType(body) == TRUEANOMALY)
		KType(body) = ECCANOMALY;
	else
		error("cast_eccanomaly(): wrong type of kepler orbit\n");

	if (KEcc(body) < 0.0)		/* inconsistent */
		error("cast_eccanomaly(): negative eccentricity\n");

	if (KEcc(body) < 1.0)		/* ellipse */
		{
		if (KTruean(body) < 0.0)
			warning("cast_eccanomaly(): true anomaly < 0\n");
		if (KTruean(body) > TWO_PI)
		       warning("cast_eccanomaly(): true anomaly > 2pi\n");
		KTruean(body) = pos_angle( KTruean(body) );
		ecc_an = acos( (KEcc(body) + cos(KTruean(body))) /
				       (1.0 + KEcc(body)*cos(KTruean(body))) );
		if ( KTruean(body) < PI )
			KEccan(body) = ecc_an;
		else
			KEccan(body) = TWO_PI - ecc_an;
		}
	else				/* hyperbola */
		{
		maxtruean = PI - acos( 1.0 / KEcc(body) );
				/*
				 *  the asymptotic value of the true anomaly
				 *  for very large separation has absolute
				 *  value   maxtruean  .
				 */
		if (KTruean(body) < -maxtruean)
			error("cast_eccanomaly(): hyp. true anomaly < min.\n");
		if (KTruean(body) > maxtruean)
			error("cast_eccanomaly(): hyp. true anomaly > max.\n");
		ecc_an = acosh( (KEcc(body) + cos(KTruean(body))) /
				       (1.0 + KEcc(body)*cos(KTruean(body))) );
		if ( KTruean(body) > 0.0 )
			KEccan(body) = ecc_an;
		else
			KEccan(body) = -ecc_an;
		}
	}

#undef	  acosh

/*-----------------------------------------------------------------------------
 *  cast_meananomaly  --  transform from true anomaly to mean anomaly
 *-----------------------------------------------------------------------------
 */
#define	  acosh(x)   log( (x) + sqrt( (x)*(x) - 1 ) )  /* inverse of cosh() */

local void  cast_meananomaly( body )
KEPLERPTR  body;
	{
	real  maxtruean;	/*  max. value for hyperbolic true anomaly */
	real  ecc_an;		/*  eccentric anomaly (ellipse), or
				 *  the equivalent (hyperbola)
				 */

	if (KType(body) == TRUEANOMALY)
		KType(body) = MEANANOMALY;
	else
		error("cast_meananomaly(): wrong type of kepler orbit\n");

	if (KEcc(body) < 0.0)		/* inconsistent */
		error("cast_meananomaly(): negative eccentricity\n");

	if (KEcc(body) < 1.0)		/* ellipse */
		{
		if (KTruean(body) < 0.0)
			warning("cast_meananomaly(): true anomaly < 0\n");
		if (KTruean(body) > TWO_PI)
		       warning("cast_meananomaly(): true anomaly > 2pi\n");
		KTruean(body) = pos_angle( KTruean(body) );
		ecc_an = acos( (KEcc(body) + cos(KTruean(body))) /
				       (1.0 + KEcc(body)*cos(KTruean(body))) );
		if ( KTruean(body) > PI )
			ecc_an = TWO_PI - ecc_an;
		}
	else				/* hyperbola */
		{
		maxtruean = PI - acos( 1.0 / KEcc(body) );
				/*
				 *  the asymptotic value of the true anomaly
				 *  for very large separation has absolute
				 *  value   maxtruean  .
				 */
		if (KTruean(body) < -maxtruean)
			error("cast_meananomaly(): hyp. true anomaly <min.\n");
		if (KTruean(body) > maxtruean)
			error("cast_meananomaly(): hyp. true anomaly >max.\n");
		ecc_an = acosh( (KEcc(body) + cos(KTruean(body))) /
				       (1.0 + KEcc(body)*cos(KTruean(body))) );
		if ( KTruean(body) < 0.0 )
			ecc_an *= -1.0;
		}
	if (KEcc(body) < 1.0)		/* ellipse */
		KMeanan(body) = ecc_an - KEcc(body) * sin(ecc_an);
	else				/* hyperbola */
		KMeanan(body) = -ecc_an + KEcc(body) * sinh(ecc_an);
	}

#undef	  acosh

/*-----------------------------------------------------------------------------
 *  cast_peripassage  --  transform from true anomaly to
 *			  time-of-pericenter-passage coordinates
 *-----------------------------------------------------------------------------
 */
#define	  acosh(x)   log( (x) + sqrt( (x)*(x) - 1 ) )  /* inverse of cosh() */

local void  cast_peripassage( body )
KEPLERPTR  body;
	{
	real  maxtruean;	/*  max. value for hyperbolic true anomaly */
	real  mtot;		/*  sum of the two masses		    */
	real  mean_an;		/*  mean anomaly			    */
	real  ecc_an;		/*  eccentric anomaly (ellipse), or
				 *  the equivalent (hyperbola)
				 */
	real  omega;		/*  orbit average of the angular velocity   */
			/*  these definitions are for an elliptic orbit;
			 *  for a hyperbola, the orbital period is infinite,
			 *  and Kepler's third law is used instead.
			 */

	if (KType(body) == TRUEANOMALY)
		KType(body) = PERIPASSAGE;
	else
		error("cast_peripassage(): wrong type of kepler orbit\n");

	if (KAsemi(body) < 0.0)		/* inconsistent */
		error("cast_peripassage(): negative semimajor axis\n");

	if (KEcc(body) < 0.0)		/* inconsistent */
		error("cast_peripassage(): negative eccentricity\n");

	if (KEcc(body) < 1.0)		/* ellipse */
		{
		if (KTruean(body) < 0.0)
			warning("cast_peripassage(): true anomaly < 0\n");
		if (KTruean(body) > TWO_PI)
		       warning("cast_peripassage(): true anomaly > 2pi\n");
		KTruean(body) = pos_angle( KTruean(body) );
		ecc_an = acos( (KEcc(body) + cos(KTruean(body))) /
				       (1.0 + KEcc(body)*cos(KTruean(body))) );
		if ( KTruean(body) > PI )
			ecc_an = TWO_PI - ecc_an;
		}
	else				/* hyperbola */
		{
		maxtruean = PI - acos( 1.0 / KEcc(body) );
				/*
				 *  the asymptotic value of the true anomaly
				 *  for very large separation has absolute
				 *  value   maxtruean  .
				 */
		if (KTruean(body) < -maxtruean)
			error("cast_peripassage(): hyp. true anomaly <min.\n");
		if (KTruean(body) > maxtruean)
			error("cast_peripassage(): hyp. true anomaly >max.\n");
		ecc_an = acosh( (KEcc(body) + cos(KTruean(body))) /
				       (1.0 + KEcc(body)*cos(KTruean(body))) );
		if ( KTruean(body) < 0.0 )
			ecc_an *= -1.0;
		}
	if (KEcc(body) < 1.0)		/* ellipse */
		mean_an = ecc_an - KEcc(body) * sin(ecc_an);
	else				/* hyperbola */
		mean_an = -ecc_an + KEcc(body) * sinh(ecc_an);
	mtot = KMass1(body) + KMass2(body);
		/*
		 *  Kepler's third law can yield a definition of
		 *  mean anomaly also in the hyperbolic case:
		 */
	omega = sqrt( mtot / (KAsemi(body)*KAsemi(body)*KAsemi(body)) );
	KPerpas(body) = KTime(body) - mean_an / omega;
	}

#undef	  acosh

/*-----------------------------------------------------------------------------
 *  cast_tscattering  --  transforms from tscattering to scattering
 *                        coordinates,
 * 			  i.e.: changing from initial-separation to
 *                              time-before-pericenter-passage.
 *                        note: for outgoing orbits, this time is negative.
 *-----------------------------------------------------------------------------
 */
local void  cast_tscattering( body )
KEPLERPTR  body;
	{
	KEPLER  temp, *p;

	p = &temp;
	scopy(p, body);

	if (KType(body) == SCATTERING)
		KType(body) = TSCATTERING;
	else
		error("cast_tscattering(): wrong type of kepler orbit\n");

        transkep(p, PERIPASSAGE);

	KImp_time(body) = KPerpas(p);
	}


/*-----------------------------------------------------------------------------
 *  scopyf  --  copies an array of characters, which can be the elementary
 *		components of structures if invoked by the scopy() macro
 *		in "kepler.h"
 *-----------------------------------------------------------------------------
 */
local void  scopyf(a, b, n)
register char  *a, *b;
register int  n;
	{
	while (--n)
		*a++ = *b++;
	}

#if defined(NONEMO)
/*			This code should now be deleted for NEMO V2
typedef  char  *charptr;


/*-----------------------------------------------------------------------------
 *  warning  --  similar to error(), but continues execution
 *		 usage: as a warning when suspicious values occur which
 *			may or may not be intentional.
 *-----------------------------------------------------------------------------
 */
local void  warning(str, p1, p2, p3, p4, p5, p6, p7, p8)
string  str;                       		   /* string to be printed */
charptr  p1, p2, p3, p4, p5, p6, p7, p8;           /* cf. error()          */
	{
	fprintf(stderr,"WARNING: ");
	fprintf(stderr, str, p1, p2, p3, p4, p5, p6, p7, p8);
	}

/* end of: keptrans.c */

/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/

/* end of: kep2kep.c */

/*
 * Cannot use the NEMO error() call, since it assumes
 *      the user interface initparam()   
 */ 
error(msg,p1,p2)
char *msg, *p1,*p2;
{                   
    fprintf(stderr,msg,p1,p2);
    exit(-1);   
}
#endif
