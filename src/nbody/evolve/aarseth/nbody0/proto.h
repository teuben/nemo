/* Mangled names for C defined functions that need to be used in Fortran */

#  define nbody0  F77_FUNC(nbody0,NBODY0)
#  define inbods  F77_FUNC(inbods,INBODS)
#  define inpars  F77_FUNC(inpars,INPARS)
#  define outbods F77_FUNC(outbods,OUTBODS)
#  define outene  F77_FUNC(outene,OUTENE)

extern void nbody0(void);
extern void inpars (int *, int *, double *, double *, double *, double *, int *, int *);
extern void inbods (int *, double *, double *, double *);
extern void outbods (double *, double *, double *, double *, int *);
extern void outene (double *, int *, double *);


#if 0

		/* here's the old stuff */

/*============================================================================*/
/*                  Fortran-to-C interface stuff:
 *  Since we only handle simple parameters here (i.e. no character
 *  variables), the only thing to be done is to properly define the name
 *  of the C routine, because it will be called as a fortran subroutine. 
 *  General rule: avoid character variables and life is fairly simple in F2C.
 *
 *  In VMS the C names stay the same
 *  in UNICOS the C name is to be in UPPER CASE.
 *  In most (BSD) unix versions C name gets an underscore (_) appended
 */

#if !defined(VMS) && !defined(aix)
# if defined(unicos)
#  define nbody0  NBODY0
#  define inbods  INBODS
#  define inpars  INPARS
#  define outbods OUTBODS
#  define outene  OUTENE
# else
#  define nbody0  nbody0_
#  define inbods  inbods_
#  define inpars  inpars_
#  define outbods outbods_
#  define outene  outene_
# endif
#endif


#endif

