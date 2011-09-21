/*
 * VECTMATH.H: include file for vector/matrix operations.
 *  1-Dec-86    Created         		Joshua Barnes  Princeton, NJ.
 * 16-sep-90	merged Nemo/Starlab header files 	PJT/PxH
 *              (note starlab does not have them anymore)
 * 22-oct-90    added some impure vector operations	PJT
 * 27-apr-92    added SMULVV				PJT
 * 16-feb-97    sqrt() from math.h			pjt
 * 23-jun-01    ZENOisms
 * 23-oct-07    added ZENOism for Koda's etude code     pjt
 * 21-sep-11    added SMINV and SMAXV                   pjt
 *
 */

#ifndef _vectmath_h
#define _vectmath_h

#include <math.h>

#ifndef THREEDIM
#  ifndef TWODIM
#    ifndef NDIM
#      define THREEDIM
#    endif
#  endif
#endif

#ifdef TWODIM
#  define NDIM 2
#endif

#ifdef THREEDIM
#  define NDIM 3
#endif

#ifndef NOTYPEDEF
    typedef real vector[NDIM], matrix[NDIM][NDIM];
#endif

/*
 * Vector operations.
 */

#define CLRV(v)			/* CLeaR Vector */			\
{									\
    register int _i;							\
    for (_i = 0; _i < NDIM; _i++)					\
        (v)[_i] = 0.0;							\
}

#define UNITV(v,j)		/* UNIT Vector */			\
{									\
    register int _i;							\
    for (_i = 0; _i < NDIM; _i++)					\
        (v)[_i] = (_i == (j) ? 1.0 : 0.0);				\
}

#define SETV(v,u)		/* SET Vector */			\
{ 									\
    register int _i; 							\
    for (_i = 0; _i < NDIM; _i++) 					\
        (v)[_i] = (u)[_i]; 						\
}

#ifdef THREEDIM

#define ADDV(v,u,w)		/* ADD Vector */			\
{									\
    register real *_vp = (v), *_up = (u), *_wp = (w);			\
    *_vp++ = (*_up++) + (*_wp++);					\
    *_vp++ = (*_up++) + (*_wp++);					\
    *_vp   = (*_up  ) + (*_wp  );					\
}

#define SUBV(v,u,w)		/* SUBtract Vector */			\
{									\
    register real *_vp = (v), *_up = (u), *_wp = (w);			\
    *_vp++ = (*_up++) - (*_wp++);					\
    *_vp++ = (*_up++) - (*_wp++);					\
    *_vp   = (*_up  ) - (*_wp  );					\
}

#define MULVS(v,u,s)		/* MULtiply Vector by Scalar */		\
{									\
    register real *_vp = (v), *_up = (u);				\
    *_vp++ = (*_up++) * (s);						\
    *_vp++ = (*_up++) * (s);						\
    *_vp   = (*_up  ) * (s);						\
}

#else

#define ADDV(v,u,w)		/* ADD Vector */			\
{									\
    register int _i;							\
    for (_i = 0; _i < NDIM; _i++)					\
        (v)[_i] = (u)[_i] + (w)[_i];					\
}

#define SUBV(v,u,w)		/* SUBtract Vector */			\
{									\
    register int _i;							\
    for (_i = 0; _i < NDIM; _i++)					\
        (v)[_i] = (u)[_i] - (w)[_i];					\
}

#define MULVS(v,u,s)		/* MULtiply Vector by Scalar */		\
{									\
    register int _i;							\
    for (_i = 0; _i < NDIM; _i++)					\
        (v)[_i] = (u)[_i] * (s);					\
}

#endif

#define DIVVS(v,u,s)		/* DIVide Vector by Scalar */		\
{									\
    register int _i;							\
    for (_i = 0; _i < NDIM; _i++)					\
        (v)[_i] = (u)[_i] / (s);					\
}

#ifdef THREEDIM

#define DOTVP(s,v,u)		/* DOT Vector Product */		\
{									\
    register real *_vp = (v), *_up = (u);				\
    (s)  = (*_vp++) * (*_up++);						\
    (s) += (*_vp++) * (*_up++);						\
    (s) += (*_vp  ) * (*_up  );						\
}

#else

#define DOTVP(s,v,u)		/* DOT Vector Product */		\
{									\
    register int _i;							\
    (s) = 0.0;								\
    for (_i = 0; _i < NDIM; _i++)					\
        (s) += (v)[_i] * (u)[_i];					\
}

#endif

#define ABSV(s,v)		/* ABSolute value of a Vector */	\
{									\
    register int _i;							\
    (s) = 0.0;								\
    for (_i = 0; _i < NDIM; _i++)					\
        (s) += (v)[_i] * (v)[_i];					\
    (s) = sqrt(s);                                                      \
}

#define DISTV(s,u,v)		/* DISTance between Vectors */        	\
{									\
    register int _i;							\
    (s) = 0.0;								\
    for (_i = 0; _i < NDIM; _i++)					\
        (s) += ((u)[_i]-(v)[_i]) * ((u)[_i]-(v)[_i]);		        \
    (s) = sqrt(s);                                                      \
}

#ifdef TWODIM

#define CROSSVP(s,v,u)		/* CROSS Vector Product */		\
{									\
    (s) = (v)[0]*(u)[1] - (v)[1]*(u)[0];				\
}

#endif

#ifdef THREEDIM

#define CROSSVP(v,u,w)		/* CROSS Vector Product */		\
{									\
    (v)[0] = (u)[1]*(w)[2] - (u)[2]*(w)[1];				\
    (v)[1] = (u)[2]*(w)[0] - (u)[0]*(w)[2];				\
    (v)[2] = (u)[0]*(w)[1] - (u)[1]*(w)[0];				\
}

#endif

#define SADDV(v,u)             /* INCrementally ADD Vector */         \
{									\
    register int _i;                                                    \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (v)[_i] += (u)[_i];                                             \
}

#define SSUBV(v,u)             /* INCrementally SUBtract Vector */    \
{									\
    register int _i;                                                    \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (v)[_i] -= (u)[_i];                                             \
}

#define SMULVS(v,s)	/* INCrementally MULtiply Vector by Scalar */	\
{									\
    register int _i;                                                    \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (v)[_i] *= (s);                                                 \
}

#define SDIVVS(v,s)	/* INCrementally DIVide Vector by Scalar */	\
{									\
    register int _i;                                                    \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (v)[_i] /= (s);                                                 \
}

#define SMINV(v,u)     /* Min value of a vector */                      \
{                                                                       \
    register int _i;                                                    \
    for (_i = 0; _i < NDIM; _i++)                                       \
      (v)[_i] = MIN( (v)[_i], (u)[_i]);                                 \
}									

#define SMAXV(v,u)     /* Max value of a vector */                      \
{                                                                       \
    register int _i;                                                    \
    for (_i = 0; _i < NDIM; _i++)                                       \
      (v)[_i] = MAX((v)[_i], (u)[_i]);                                  \
}									

/* For compatiblity with older Nemo (1.3 and before) the INC macros */
/* are equated to the Self-macros (which came from starlab)	    */
#define INCADDV  SADDV
#define INCSUBV  SSUBV
#define INCMULVS SMULVS
#define INCDIVVS SDIVVS

/*
 * Matrix operations.
 */

#define CLRM(p)			/* CLeaR Matrix */			\
{									\
    register int _i, _j;						\
    for (_i = 0; _i < NDIM; _i++)					\
        for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = 0.0;						\
}

#define SETMI(p)		/* SET Matrix to Identity */		\
{									\
    register int _i, _j;						\
    for (_i = 0; _i < NDIM; _i++)					\
        for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (_i == _j ? 1.0 : 0.0);			\
}

#define SETM(p,q)		/* SET Matrix */			\
{									\
    register int _i, _j;						\
    for (_i = 0; _i < NDIM; _i++)					\
        for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (q)[_i][_j];					\
}

#define TRANM(p,q)		/* TRANspose Matrix */			\
{									\
    register int _i, _j;						\
    for (_i = 0; _i < NDIM; _i++)					\
        for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (q)[_j][_i];					\
}

#define ADDM(p,q,r)		/* ADD Matrix */			\
{									\
    register int _i, _j;						\
    for (_i = 0; _i < NDIM; _i++)					\
        for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (q)[_i][_j] + (r)[_i][_j];			\
}

#define SUBM(p,q,r)		/* SUBtract Matrix */			\
{									\
    register int _i, _j;						\
    for (_i = 0; _i < NDIM; _i++)					\
        for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (q)[_i][_j] - (r)[_i][_j];			\
}

#define MULM(p,q,r)		/* Multiply Matrix */			\
{									\
    register int _i, _j, _k;						\
    for (_i = 0; _i < NDIM; _i++)					\
	for (_j = 0; _j < NDIM; _j++) {					\
	    (p)[_i][_j] = 0.0;						\
            for (_k = 0; _k < NDIM; _k++)				\
		(p)[_i][_j] += (q)[_i][_k] * (r)[_k][_j];		\
        }								\
}

#define MULMS(p,q,s)		/* MULtiply Matrix by Scalar */		\
{									\
    register int _i, _j;						\
    for (_i = 0; _i < NDIM; _i++)				        \
        for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (q)[_i][_j] * (s);				\
}

#define DIVMS(p,q,s)		/* DIVide Matrix by Scalar */		\
{									\
    register int _i, _j;						\
    for (_i = 0; _i < NDIM; _i++)					\
        for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (q)[_i][_j] / (s);				\
}

#define MULMV(v,p,u)		/* MULtiply Matrix by Vector */		\
{									\
    register int _i, _j;						\
    for (_i = 0; _i < NDIM; _i++) {					\
	(v)[_i] = 0.0;							\
	for (_j = 0; _j < NDIM; _j++)					\
	    (v)[_i] += (p)[_i][_j] * (u)[_j];				\
    }									\
}

#define OUTVP(p,v,u)		/* OUTer Vector Product */		\
{									\
    register int _i, _j;						\
    for (_i = 0; _i < NDIM; _i++)					\
        for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (v)[_i] * (u)[_j];				\
}

#define TRACEM(s,p)		/* TRACE of Matrix */			\
{									\
    register int _i;							\
    (s) = 0.0;								\
    for (_i = 0.0; _i < NDIM; _i++)					\
	(s) += (p)[_i][_i];						\
}

/*
 * Scalar-valued vector functions.
 */

#if defined(DOUBLEPREC)

#define dotvp(v,u)	(_dotvp_d((v),(u),NDIM))
#define absv(v)		(_absv_d((v),NDIM))
#define distv(v,u)	(_distv_d((v),(u),NDIM))
#define tracem(p)	(_tracem_d((p),NDIM))

extern   double _dotvp_d(real *, real *, int);
extern   double _absv_d(real *, int);
extern   double _distv_d(real *, real *, int);
extern   double _tracem_d(real *, int);

#endif


#if defined(SINGLEPREC)

#define dotvp(v,u)	(_dotvp_f((v),(u),NDIM))
#define absv(v)		(_absv_f((v),NDIM))
#define distv(v,u)	(_distv_f((v),(u),NDIM))
#define tracem(p)	(_tracem_f((p),NDIM))

extern   double _dotvp_f(real *, real *, int);
extern   double _absv_f(real *, int);
extern   double _distv_f(real *, real *, int);
extern   double _tracem_f(real *, int);

#endif

/*
 * Misc. impure operations.
 */

#define SETVS(v,s)		/* SET Vector to Scalar */		\
{									\
    register int _i;							\
    for (_i = 0; _i < NDIM; _i++)					\
        (v)[_i] = (s);							\
}

#define ADDVS(v,u,s)		/* ADD Vector and Scalar */		\
{									\
    register int _i;							\
    for (_i = 0; _i < NDIM; _i++)					\
        (v)[_i] = (u)[_i] + (s);					\
}

#define SETMS(p,s)		/* SET Matrix to Scalar */		\
{									\
    register int _i, _j;						\
    for (_i = 0; _i < NDIM; _i++)					\
        for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (s);						\
}

#define MULVV(w,v,u)		/* MUL Vector and Vector */		\
{									\
    register int _i;							\
    for (_i = 0; _i < NDIM; _i++)					\
        (w)[_i] = (v)[_i] * (u)[_i];                                    \
}

#define SMULVV(v,u)   /* INCrementally MULtiply Vector and Vector */    \
{									\
    register int _i;                                                    \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (v)[_i] *= (u)[_i];                                             \
}


#define SQRTV(u)		/* SQRT of a Vector */			\
{									\
    register int _i;							\
    for (_i = 0; _i < NDIM; _i++)					\
        (u)[_i] = sqrt( (u)[_i] );                                      \
}

/* ZENO things to be added properly later */

#if defined(THREEDIM)

#define ADDMULVS(v,u,s)         /* MUL Vect by Scalar, ADD to vect */   \
{                                                                       \
    (v)[0] += (u)[0] * (s);                                             \
    (v)[1] += (u)[1] * (s);                                             \
    (v)[2] += (u)[2] * (s);                                             \
}

#define DOTPSUBV(s,v,u,w)       /* SUB Vectors, form DOT Prod */        \
{                                                                       \
    (v)[0] = (u)[0] - (w)[0];    (s)  = (v)[0] * (v)[0];                \
    (v)[1] = (u)[1] - (w)[1];    (s) += (v)[1] * (v)[1];                \
    (v)[2] = (u)[2] - (w)[2];    (s) += (v)[2] * (v)[2];                \
}

#define DOTPMULMV(s,v,p,u)      /* MUL Mat by Vect, form DOT Prod */    \
{                                                                       \
    DOTVP(v[0], p[0], u);    (s)  = (v)[0] * (u)[0];                    \
    DOTVP(v[1], p[1], u);    (s) += (v)[1] * (u)[1];                    \
    DOTVP(v[2], p[2], u);    (s) += (v)[2] * (u)[2];                    \
}

#define ADDMULVS2(v,u,s,w,r)    /* 2 times MUL V by S, ADD to vect */   \
{                                                                       \
    (v)[0] += (u)[0] * (s) + (w)[0] * (r);                              \
    (v)[1] += (u)[1] * (s) + (w)[1] * (r);                              \
    (v)[2] += (u)[2] * (s) + (w)[2] * (r);                              \
}

#define ADDVMULVS(v,u,w,s)      /* MUL Vect by Scalar, ADD to Vects */ \
{                                                                      \
    (v)[0] = (u)[0] + (w)[0] * (s);                                    \
    (v)[1] = (u)[1] + (w)[1] * (s);                                    \
    (v)[2] = (u)[2] + (w)[2] * (s);                                    \
}


#endif

#endif
