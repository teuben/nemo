/*
 *  DISSIPATE:   dissipate energy: turn some random motion into bulk motion
 *
 *		Peter Teuben - february 1989
 *          2-nov-91    doc + more output       PJT
 *	   17-jun-92    attempt to conserve energy better	PJT
 *			(what we had so far, was momentum conservation)	PJT
 *         29-sep-92    added 'grid' parameter to set largest possible grid PJT
 *         25-oct-92    added fheat
 *	   10-apr-96	fixed  calloc() delcaration
 *         10-apr-01    gcc warnings
 */

#include "defs.h"
#ifndef HUGE
#define HUGE 1E20
#endif

static int *c = NULL;
static int    size = 0;
static int    entry = 0;

#define ECONS  1		/* flag energy conservation */
#define USE_MAXGRID 1           /* fix max allowed grid */


void dissipate (btab, nb, ndim, dr, eta, grid, fheat)
Body *btab;
int   nb;
int   ndim;
real *dr;
real  eta;
real  grid;
real  fheat;
{
    real rmin[NDIM], rmax[NDIM], *pos;
    int  nbin[NDIM], ndis=0;
    Body *b, *b1;
    int   i, k, n, ix, iy, iz, nx, ny, nz, nxyz, ind;
    bool Qheat, Qangle, Qkappa, scanopt();
    vector  velsum, veldif;
    real t_before, t_after, kappa;
    real angle, p, ss, cc, velsig, vx, vy;
    
    if (eta==0.0) return;          /* no work to do ... */
    Qheat = (fheat > 0) ;
    Qangle = scanopt(options,"angle");
    Qkappa = scanopt(options,"kappa");
    entry++;                                    /* debug counter of entries */
    	
    for (i=0; i<NDIM; i++) {                    /* init min and max of cube */
        rmin[i] = HUGE;
	rmax[i] = -HUGE;
    }
    
    for (b=btab; b<btab+nb; b++) {              /* get min and max of cube */
        pos = Pos(b);
	for (i=0; i<NDIM; i++) {
	    rmin[i] = MIN(rmin[i], pos[i]);
            rmax[i] = MAX(rmax[i], pos[i]);
	}
    }

#ifdef USE_MAXGRID
    if (grid > 0)
        for (i=0; i<NDIM; i++) {
            if (rmin[i] < -grid) {
                dprintf(1,"Reset rmin[%d] from %g to %g\n",i,rmin[i],-grid);
                rmin[i] = -grid;
            }
            if (rmax[i] >  grid) {
                dprintf(1,"Reset rmax[%d] from %g to %g\n",i,rmax[i],grid);
                rmax[i] =  grid;
            }
        }
#endif
    
    for (i=0; i<NDIM; i++)                      /* get cell size */
        if (dr[i] > 0)
	    nbin[i] = (rmax[i]-rmin[i])/dr[i] + 1;
	else
	    nbin[i] = 1;
	    
    nx = nbin[0];
    ny = nbin[1];
    nz = nbin[2];
    nxyz = nx*ny*nz;
    if (nxyz > size) {                      /* need more space !! */
        if (c)
	    free(c);                        /* allocate old one */
        c = (int *) calloc(nxyz, sizeof(int));
	if (c==NULL) {
            warning("No memory for %d * %d * %d cube; no dissipation",nx,ny,nz);
	    return;       /* error: not enough memory */
        }
	size = nxyz;                        /* and remember new space */
	dprintf(1,"Allocated %d on entry # %d\n",size,entry);
    } else
	for (i=0; i<nxyz; i++)
	    c[i] = 0;                     /* reset index cube */

    for (b=btab, i=0; b<btab+nb; b++, i++) {    /* build interaction list */
        ix = (Phase(b)[0][0] - rmin[0])/dr[0];  /* get cell location */
        iy = (Phase(b)[0][1] - rmin[1])/dr[1];
        iz = (Phase(b)[0][2] - rmin[2])/dr[2];
        if (ix<0 || ix>=nx) continue;       /* allow particles outside grid */
        if (iy<0 || iy>=ny) continue;
        if (iz<0 || iz>=nz) continue;
        ind = ix + nx*iy + nx*ny*iz;        /* index into c[] array */
	if ((k = c[ind]) == 0) {          /* no star in cell yet */
	    c[ind] = i+1;                 /* nonzero 1..nb index */
	} else {                            /* already star in cell */
            b1 = btab + k - 1;              /* this is that star */
	    while (Key(b1))
	        b1 = btab + (Key(b1) - 1);  /* next star */
	    Key(b1) = i + 1;                /* point to current star */
	}
	Key(b) = 0;                         /* terminate list */
    }


    for (ix=0; ix<nx; ix++)                 /* walk through cube */
    for (iy=0; iy<ny; iy++)
    for (iz=0; iz<nz; iz++) {
        ind = ix+nx*iy+nx*ny*iz;            /* offset in cube */
	k = c[ind];                         /* ask first star in cell */
        if (!k)
            continue;                       /* empty cell */

	CLRV(velsum);	                    /* reset */
        t_before = t_after = 0.0;
	n = 0;
	while (k) {
	    b = btab + k - 1;               /* point to star */
	    ADDV(velsum,velsum,Vel(b));     /* accumulate mean cell velo */
	    t_before += dotvp(Vel(b),Vel(b));/* kinetic before */
	    n++;                            /* number of stars in cell */
	    k = Key(b);                     /* next star in cell ? */
	}
	if (n<2)
	    continue;                       /* no need to average */
        MULVS(velsum,velsum,1.0/n);         /* get average velocity in cell */
	ndis++;				    /* count dissipative cells */	

	k = c[ind];                         /* start again in cell */
        velsig = 0;
	while (k) {
	    b = btab + k - 1;
	    SUBV(veldif,velsum,Vel(b));     /* get difference from mean */
            if(Qheat) velsig += dotvp(veldif,veldif);
	    MULVS(veldif,veldif,eta);       /* take a fraction */
	    ADDV(Vel(b),Vel(b),veldif);     /* add it to velocity */
            t_after += dotvp(Vel(b),Vel(b));/* kinetic after */
	    k = Key(b);                     /* point to next star */
        }

        if (Qheat) {
            angle = fheat * sqrt(velsig/n);
            k = c[ind];
            while(k) {
                b=btab+k-1;
		if (Qangle) Aux(b) = angle;
                p=grandom(0.0,angle);
                ss=sin(p);   cc=cos(p);
                vx = Vel(b)[0] * cc  -  Vel(b)[1] * ss;   /* rotate the vector */
                vy = Vel(b)[0] * ss  +  Vel(b)[1] * cc;
                Vel(b)[0] = vx;
                Vel(b)[1] = vy;
		k=Key(b);
            }
        }

        kappa = sqrt(t_before/t_after);
/**/	kappa=1;	/**PJT**/
        k = c[ind];                         /* final renormalization */
	while (k) {			    /* loop again over all */
	    b = btab + k - 1;		    /* stars in the cell */
            if (Qkappa) Aux(b) = kappa;
            SMULVS(Vel(b),kappa);	    /* correct amplitude */
	    k = Key(b);                     /* point to next star */
	    Key(b) = 0;                     /* reset Key - not needed */
        }
    }  /* end loop cube (ix,iy,iz) */
    dprintf(1,"Dissipating %d in %d*%d*%d=%d cells\n",ndis,nx,ny,nz,nxyz);
    if (ndis==0) 
        warning("No dissipation done, cell=%g or nbody=%d too small?",
                    dr, nb);
}
