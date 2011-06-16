/* 
 * for NEMO, mostly taken from $TIPSY/code/defs.h
 *
 *	2-sep-94	my last mod			PJT
 *	17-aug-00	fixed padding byte problem	PJT
 *			-DNEEDPAD  for now .....
 *      28-aug-01       no NEEDPAD but warn if sizeof() not 32.....????
 *      16-jun-11       trying the new Balin padding scheme         PJT
 *                      but enabling TIPSY_NEEDPAD again
 */

#define MAXDIM 3
#define TIPSY_NEEDPAD   
#define forever for(;;)

typedef float Real;

struct gas_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real rho;
    Real temp;
    Real hsmooth;
    Real metals ;
    Real phi ;
} ;

struct gas_particle *gas_particles;

struct dark_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real eps;
    Real phi ;
} ;

struct dark_particle *dark_particles;

struct star_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real metals ;
    Real tform ;
    Real eps;
    Real phi ;
} ;

struct star_particle *star_particles;

struct dump {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
#ifdef TIPSY_NEEDPAD
    int junk ;			/* padding byte  !!! ??? !!! */
#endif
#if 0
        /* Jeremy Balin addition */
    char align[ (32 - sizeof(double) - 5 * sizeof(int)) / sizeof(char) ];
    /* total size should be 32 bytes to make alignment okay with
     * 64-bit architectures (ie. alphas) */
#endif
} ;

struct dump header ;
