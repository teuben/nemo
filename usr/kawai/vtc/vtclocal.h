#ifndef __VTCLOCAL_H__
#define __VTCLOCAL_H__

/*
 * global definitions
 */
#ifdef DEBUG
#undef DEBUG
#define DEBUG (1)
#endif

#if DEBUG
#define Cprintf printf
#define Cfprintf fprintf
#else /* !DEBUG */
#define Cprintf
#define Cfprintf
#endif /* DEBUG */

#define PR(x,c)  fprintf(stderr, #x" = %"#c" ", x);
#define PRL(x,c)  fprintf(stderr, #x" = %"#c"\n", x);

#define TRUE (1)
#define FALSE (0)

#if (GRAPE5+MDGRAPE2+NOGRAPE) != 1
#error specify one of GRAPE5, MDGRAPE2, and NOGRAPE
#endif

#define NBOARDMAX (8) /* max # of GRAPE processor board */

/* max # of j-particles can be stored on a single board */

#if GRAPE5
#include <gp5util.h>
#define NJMAX (JMEMSIZE)
#endif /* GRAPE5 */

#if MDGRAPE2
#include <mdgp2defs.h>
#define NJMAX (MDGP2_BOARD_MAXJ)
#endif /* MDGRAPE2 */

#if NOGRAPE
#define NJMAX (20000)
#endif /* JMEMSIZE */

/* body attrs */
#if defined(__APPLE__) || defined(__i386__) || defined(__SVR4) /* 32-bit long platform (ex. Darwin, Linux, Solaris) */
typedef long long Mortonkey;
#else /* 64-bit long platform (ex. Tru64, UXP/V on VPP5k) */
typedef long Mortonkey;
#endif
#define MORTONKEY_LENGTH (16)
#define MORTONKEY_OFFSET (((Mortonkey)1)<<(MORTONKEY_LENGTH-1))
typedef int Body;

/* cell attrs */
#define NOCELL (-1)
typedef int Cell;

/* opaque (Body or Cell) */
typedef int Node;

#ifndef _FORCE_C_
extern void (*vtc_force_calculator[])(int ni, double (*xi)[3], int nj,
				      double (*xj)[3], double *mj,
				      double eps, double (*a)[3], double *p);
#endif /* _FORCE_C_ */

extern int vtc_get_grape(void);

/* timing analysis */
void init_cputime(void);
double get_cputime(void);

/* P2M2 */
void load_design(int order, double ss, int full_dof, int negativemass,
		 int *npp, double (**pppos)[3]);
void plgndr0(int n, double x, double *pln);
double plgndr(int l, int m, double x);
void ylm(int l, int m, double theta, double phi, double *re, double *im);
int jacobi(double (*a)[3], double d[3], double (*v)[3], int *nrot);
void eigenvalsort(double d[3], double (*v)[3]);

/* sorting */
void sort_body(Body *b, Mortonkey *key, int n);

#endif /* __VTCLOCAL_H__ */
