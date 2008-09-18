/*
 *  POTENTIAL: 
 *	interfaces to aid in dynamically loading potentials
 *
 *	jul 1987:	original implementation
 *	sep 2001:	added C++ support, including const'ing 
 */

#ifndef _potential_h
#define _potential_h


typedef struct a_potential {
    string name;
    string pars;
    string file;
    struct a_potential *next;
} a_potential;

/* should we make args 1,2 and 5 const ?? */

typedef void (*potproc_double)(const int *, const double *, double *, double *, const double *);
typedef void (*potproc_float) (const int *, const float *,  float *,  float *,  const float *);
#ifdef SINGLEPREC
typedef potproc_float potproc_real;
#else
typedef potproc_double potproc_real;
#endif

#if defined(__cplusplus)
extern "C" {
#endif


potproc_real   get_potential        (const string, const string, const string);
potproc_float  get_potential_float  (const string, const string, const string);
potproc_double get_potential_double (const string, const string, const string);
proc           get_inipotential     (void);
real           get_pattern          (void);

#if defined(__cplusplus)
}
#endif


#endif
