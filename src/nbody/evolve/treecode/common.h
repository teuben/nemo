/*
 *  common.h:   definition of the common blocks used in the fortran code
 *		assumed NMAX and NDIM have been defined before 
 */
#define NBODCELL  (NMAX+NMAX)



struct fcb1 {                           /* PARAMCOM */
    int     nbodies;
    float   tol, eps, eps2, epsinv;
    int     usequad; /* logical !! */
};

struct fcb2 {                           /* MSGCOM */
    char    headline[20];
};

struct fcb3 {                           /* CELLCOM */
    real    rsize, rmin[NDIM];
    int     incells, imax, imax2, nindex[NDIM];
};
#if 0
struct fcb4 {                           /* POINTERS */
    int     toot, pskip, subp[NSUBCELL*NJUNK];
};
#endif

struct fcb5 {
    float   mass[NBODCELL], phi[NMAX];
    float   pos[NDIM*NBODCELL], sizetol2[NBODCELL];
    float   vel[NDIM*NMAX], acc[NDIM*NMAX];
};

#if 0
struct fcb6 {
};
#endif

struct fcb7 {                               /* FORCECOM */
    int     nterms, nttot, ntmin, ntmax, ntavg;
};

struct fcb8 {                               /* TIMECOM */
    int     nsteps, nout;
    float   tnow, tpos, tvel, dtime, dtime2;
};

struct fcb9 {                               /* CPUCOM */
    float   cputime0, cputime1, cputime;
};

#if 0
struct fcb10 {                              /* QUADCOM */
};
#endif

#if 0
struct fcb11 {                              /* TEMPQUAD */
};
#endif
