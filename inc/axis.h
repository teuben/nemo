/*
 * AXIS.H: definitions for axis routines.
 */
#ifndef _axis_h
#define _axis_h
/*
 * XAXISVAR, YAXISVAR: collected parameters for x, y axies.
 */

struct axisvars {
    int ndig;			/* max no. of digits after dec. pnt. */
    real tikup, tikdn;		/* axis tickmark params: plus and minus */
    real numdn, sznum;		/* tick-value positioning and size params */
    real labdn, szlab;		/* label positioning and size params */
};

extern struct axisvars xaxisvar;
extern struct axisvars yaxisvar;

/*
 * FORMALAXIS: if TRUE, draw axies in formal style.
 */

extern bool formalaxis;

typedef real (*axis_proc)(real);

extern void xaxis (real, real, real, real *, int, axis_proc, string);
extern void yaxis (real, real, real, real *, int, axis_proc, string);
extern void axvar (real, real, real);
extern void xaxvar(int, real, real, real, real);
extern void yaxvar(int, real, real, real, real);



#endif


/* axis.c */
