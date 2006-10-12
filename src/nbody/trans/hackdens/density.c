/*
 * DENSITY.C: routines to compute local density.
 * Public routines: hackden().
 *
 *      18-jul-92  PJT  replaced many if(debug)printf(...) by dprintf(1,...)
 *	 1-apr-01  PJT  compiler warnings
 *      15-sep-06  WD   compiler error (in gcc-3.4.5)/warning (otherwise)
 *      20-oct-06  PJT  removed all old style declarations, all local routines
 */

#include "defs.h"

local void walksub(real, nodeptr, vector, real, real *, int *);
local bool subdivp(real, vector, real);
local real distcount(real *, int, int);

real directden(p, nb, dis, ra, base, nbody)
    bodyptr p;
    int nb;
    real dis;
    real *ra;
    bodyptr base;
    int nbody;
{
    int total, i;
    bodyptr q;
    real rn, nbr, den;
    vector disp;
    for(i=0, q=base; i<nbody; i++, q++){
	SUBV(disp, Pos(p), Pos(q));
	DOTVP(ra[i], disp, disp);
    }
    rn=distcount(ra, nbody, nb);
    nbr = nb-2.0;
    den=nbr/(rn*sqrt(rn)*FRTHRD_PI);
#ifdef DEBUG
    dprintf(0,"Directden= %f\n", den);
#endif
    return (den);
}

real hackden(p, nb, dis, newdis, ra)
    bodyptr p;
    int nb;
    real dis;
    real *newdis;
    real *ra;
{
    int total;
    real rn, nbr, den;
    real dis0, dismax;
#ifdef DEBUG
    dprintf(0,"Hackden: nb, dis: %d %f\n", nb, dis);
#endif    
    total=0; dis0=0; dismax=1.1e30;
    while(total < nb || total > nb*10){
	hackcount(p, dis, ra, &total);
	if(total < nb) {
	   dis0=dis;
	   if(dismax < 1e30){
	      dis=(dismax+dis0)*0.5;
	   }else{
	      dis=dis*1.5;
	   }
	}
	if(total > nb*10){
	   dismax=dis;
	   dis=(dis+dis0)*0.5;
	}
#ifdef DEBUG
	dprintf(0,"nb, newdis: %d %f\n", total, dis);
#endif

     }
    nbr=nb;
    *newdis=1.5*dis*pow(nbr/total,0.333333);
#ifdef DEBUG
    dprintf(0,"Hackcount returns: %d\n", total);
#endif	
    rn=distcount(ra, total, nb);
#ifdef DEBUG
    dprintf(0,"distcount returns: %f\n", rn);
#endif
    if (Qdensity) {
      nbr = nb-2.0;
      den=nbr/(rn*sqrt(rn)*FRTHRD_PI);
    } else
      den = rn;

#ifdef DEBUG
    dprintf(0,"Hackden= %f\n", den);
#endif
    return den;
}

local real distcount(ra ,total, nb)
    real *ra;
    int total;
    int nb;
{
    register int i,j;
    register real tmp;
    int jmin;
#ifdef DEBUG
    dprintf(0,"distcount: distances--");
    for(i=0; i<total; i++)dprintf(0," %f", ra[i]);
    puts("");
#endif    
    for(i=0; i<nb; i++){
	tmp=1e20;
	for(j=i; j<total; j++){
	    if(ra[j] < tmp){
		tmp=ra[j];
		jmin=j;
	    }
	}
	ra[jmin]=ra[i];
	ra[i]=tmp;
    }
#ifdef DEBUG
    dprintf("after sort: distances--");
    for(i=0; i<nb; i++)dprintf(0," %f", ra[i]);
    puts("");
#endif    
    return(ra[nb-1]);
}
/*
 * HACKCOUNT: count the number of particles in a given radius.
 */

local bodyptr pskip;			/* body to skip in force evaluation */
local vector pos0;			/* point to evaluate field at */

hackcount(p, dis, ra, total)
bodyptr p;
real dis;
real *ra;
int * total;
{
    pskip = p;					/* exclude p from f.c.      */
    SETV(pos0, Pos(p));				/* set field point          */
    hackwalk(dis, ra, total);			/* recursively compute      */
}

/*
 * HACKWALK: walk the tree opening cells too close to a given point.
 */


hackwalk(dis, ra, total)
real dis;
real *ra;
int * total;
{
    vector croot;
    int i;
    for(i=0; i<NDIM; i++)croot[i]=rmin[i]+rsize*0.5;
    *total=0;
    walksub(dis, troot, croot, rsize, ra, total);
}

/*
 * WALKSUB: recursive routine to do hackwalk operation.
 */

local void walksub(dis, p, cpos, d, ra, total)
real dis;			        /* critical displacement */
register nodeptr p;                     /* pointer into body-tree */
vector cpos;			        /* geometoric center of the node */
real d;                                 /* size of box  */
real * ra;			/* array to store distances to */
				/* particles within sphere */
int *total;			/* number of particles in the sphere */
{
    register nodeptr *pp;
    register int i,j;
    register int k;
    real offset, r2;
    vector cpossub, disp;
    offset = d*0.25;
    dprintf(1,"walksub: p = %o  d = %f\n", p, d);
    if (Type(p) == BODY){
	r2=0.0;
	SUBV(disp, Pos(p), pos0);               /* compute displacement     */
	DOTVP(r2, disp, disp);                  /* and find dist squared    */
	if(r2 < dis*dis){
	    *total +=1;
	    ra[*total-1]=r2;
	}
    }else if (subdivp(dis, cpos, d)) {          /* should p be opened?      */
        pp = & Subp(p)[0];                      /*   point to sub-cells     */
        for (k = 0; k < NSUB; k++) {            /*   loop over sub-cells    */
	    for(i=NDIM-1, j=1; i>=0; i--, j*=2){
		if(j&k){
		    cpossub[i]=cpos[i]+offset;
		}else{
		    cpossub[i]=cpos[i]-offset;
		}
	    }
            if (*pp != NULL)                    /*     does this one exist? */
                walksub(dis,*pp, cpossub, d*0.5, ra, total);
            pp++;                               /*     point to next one    */
        }
    }
}

/*
 * SUBDIVP: decide if a node should be opened.
 * true if need to subdivide
 */

local bool subdivp(dis, cpos, d)
real dis;			        /* critical separation  */
vector cpos;			        /* geometrical center of the node */
real d;                                 /* size of cell squared */
{
    int i;
    vector dr;
    real drsq, lcrit;
    SUBV(dr, cpos, pos0);                     /* compute displacement     */
    for (i=0; i<3; i++){
	if(ABS(dr[i]) > dis+d*0.5){
	    return (0);
	}
    }
    DOTVP(drsq, dr, dr);                        /* and find dist squared    */
    lcrit= dis + 0.875*d;	                /* critical separation */
    lcrit = lcrit*lcrit;
    return (drsq < lcrit);                /* use geometrical rule     */
}
