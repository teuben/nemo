/*
 *  POTENTIAL: 
 *	interfaces to aid in dynamically loading potentials
 *
 *	jul 1987:	original implementation
 *	sep 2001:	added C++ support
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

typedef void (*potproc)       (int *, double *, double *, double *, double *);
typedef void (*potproc_double)(int *, double *, double *, double *, double *);
typedef void (*potproc_float) (int *,  float *,  float *,  float *,  float *);


#if defined(__cplusplus)
extern "C" {
#endif


potproc        get_potential        (string, string, string);
potproc_float  get_potential_float  (string, string, string);
potproc_double get_potential_double (string, string, string);
potproc        get_inipotential     (void);
real           get_pattern          (void);

#if defined(__cplusplus)
}
#endif


#endif
