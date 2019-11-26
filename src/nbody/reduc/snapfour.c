/*
 * SNAPFOUR.C:   fourier coefficients of an N-body distribution
 *
 *      30-nov-90       V1.0    Created - after Kevin Long's talk       PJT
 *	20-dec-90	V1.0b	added printed header			PJT
 *      17-feb-92       V1.1    added weight=                           PJT
 *	22-feb-92	V1.1b   usage
 *       9-nov-93       V1.2    times=
 *       7-may-02       minor code cleanup
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>  
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

string defv[] = {
    "in=???\n              Input snapshot",
    "radii=0:2:0.1\n       Set of radii denoting edges of cylinders",
    "cos=0:4:1\n	   List of noj-zero cos(m.phi) terms",
    "sin=1:4:1\n	   List of non-zero sin(m.phi) terms",
    "xvar=x\n              X variable",
    "yvar=y\n              Y variable",
    "fvar=vy\n             Fourier Observable to be decomposed",
    "weight=1\n            Weight applied to observable",
    "amode=t\n             Display sin/cos amps or amp/phase if possible?",
    "times=all\n           Snapshots to select",
    "VERSION=1.2c\n        25-nov-2019 PJT",
    NULL,
};

string usage = "Fourier coefficients of an N-body distribution";

#define TIMEFUZZ        0.0001  /* tolerance in time comparisons */

#define MAXRAD 513
#define MAXORDER 8

nemo_main()
{
    stream instr;
    string times, vstr;
    Body   *btab = NULL, *bp;
    int    i, n, nbody, bits, nrad, maxorder, tmpi[MAXORDER+1];
    real rad2[MAXRAD], tsnap;
    bool Qcos[MAXORDER+1], Qsin[MAXORDER+1], amode;
    rproc btrtrans(), xproc, yproc, fproc, wproc;

    times = getparam("times");
    nrad = nemoinpr(getparam("radii"),rad2,MAXRAD);     /* get radii */
    for (i=0; i<nrad; i++)
        rad2[i] = sqr(rad2[i]);             /* but actually save the square */

    for (i=0; i<=MAXORDER; i++) {       /* initially set all coefs to false */
        Qcos[i] = FALSE;
        Qsin[i] = FALSE;        /* Qsin[0] also set but never used though */
    }
    n = nemoinpi(getparam("cos"),tmpi,MAXORDER+1);  /* get true cos coefs */
    for (i=0; i<n; i++)
        if (tmpi[i]<0 || tmpi[i]>MAXORDER)
            warning("Illegal value %d for cos= skipped",tmpi[i]);
        else
            Qcos[tmpi[i]] = TRUE;
    n = nemoinpi(getparam("sin"),tmpi,MAXORDER+1);  /* get true sin coefs */
    for (i=0; i<n; i++)
        if (tmpi[i]<0 || tmpi[i]>MAXORDER)
            warning("Illegal value %d for sin= skipped",tmpi[i]);
        else
            Qsin[tmpi[i]] = TRUE;
    for (i=0, maxorder=-1; i<=MAXORDER; i++)
        if (Qsin[i] || Qcos[i]) maxorder=i;
    if (maxorder<0) error("No true sin or cos coefficients supplied");
    xproc = btrtrans(getparam("xvar"));
    yproc = btrtrans(getparam("yvar"));
    fproc = btrtrans(getparam("fvar"));
    wproc = btrtrans(getparam("weight"));
    amode = getbparam("amode");
        
    instr = stropen(getparam("in"), "r");           /* open input file */
    get_history(instr);                         /* get history */
    for (;;) {                          /* loop through snapshots */
        if (!get_tag_ok(instr, SnapShotTag))
                break;                           /* until done */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (!streq(times,"all") && !within(tsnap,times,0.0001))
            continue;                   /* skip work on this snapshot */
        if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0) {
            dprintf (2,"Time= %f auto skipping ",tsnap);
            continue;       /* just skip - it maybe diagnostics */
        }
        dprintf (2,"Time= %f ",tsnap);
        snap_four(btab,nbody,tsnap,xproc,yproc,fproc,wproc,
                maxorder,Qcos,Qsin,rad2,nrad,amode);
    }
}


snap_four(btab,nbody,tsnap,xproc,yproc,fproc,wproc,
          maxorder,Qcos,Qsin,rad,nrad,amode)
Body *btab;                 /* pointer to snspshot with nbody Bodie's */
real rad[];                 /* radii for shells */
real tsnap;                 /* time of snapshot */
rproc xproc, yproc, fproc;  /* procedures to compute radius, angle and f */
rproc wproc;                /* weight factor per data point */
bool Qcos[], Qsin[];        /* designate if coef to be used */
int nbody, maxorder, nrad;     
bool amode;                 /* TRUE=amps only FALSE=amp+phase if all available */
{
    real   th, r2, v, rsum, vsum, cosk, amp, pha, radius, w;
    int    i,k,m,cnt,dim,ip;
    Body *bp;
    real mat[2*(MAXORDER+1)*(MAXORDER+1)],vec[2*(MAXORDER+1)];
    real sol[2*(MAXORDER+1)], a[2*(MAXORDER+1)+1];
    permanent bool first=TRUE;

    for (m=0, dim=0; m<=maxorder; m++)  /* count dimension of matrix needed */
        if (Qcos[m]) dim++;
    for (m=1; m<=maxorder; m++)
        if (Qsin[m]) dim++;
    if (dim==0) {
       warning("snap_four: dim=0");
       return 0;
    }
    dprintf(1,"snap_four: maxorder=%d dim=%d\n",maxorder,dim);
    if (!amode) {       /* check if OK to do phases and amplitudes */
        for (m=1; m<=maxorder; m++) {
            if (Qcos[m] && !Qsin[m]) amode=TRUE;
            if (!Qcos[m] && Qsin[m]) amode=TRUE;
        }
        if (amode) 
            warning("amode=f requested, but missing cos/sin terms");
    }
    if (first) {
        print_header(maxorder,Qcos,Qsin,amode);
        first = FALSE;
    }

    for (i=1; i<nrad; i++) {            /* foreach ring */
        lsq_zero(dim,mat,vec);          /* reset accum. matrix and vector */
        cnt = 0;                        /* count points in this ring */
        for (ip=0, bp=btab; ip<nbody; ip++, bp++) {  /* loop for all bodies */
            w = wproc(bp,tsnap,ip);
            r2 = sqr(xproc(bp,tsnap,ip)) + sqr(yproc(bp,tsnap,ip));
            if (r2<rad[i-1] || r2>rad[i])          /* if not in ring:  */
                continue;                           /* skip this particle */
            cnt++;
            th = atan2(yproc(bp,tsnap,ip) , xproc(bp,tsnap,ip));
            k=0;                            /* always count how many coefs */
            for(m=0; m<=maxorder; m++) {     /* cos(m.theta) */
                if (!Qcos[m]) continue;
                a[k++] = cos(m*th);
            }
            for(m=1; m<=maxorder; m++) {     /* sin(m.theta) */
                if (!Qsin[m]) continue;
                a[k++] = sin(m*th);
            }
            a[k] = fproc(bp,tsnap,ip);
            dprintf(1,"adding %d: r^2=%g th=%g, fvar=%g wt=%g\n",
			      cnt,r2,th*180/PI,a[k],w);
            if (k!=dim) error("snapfour: Counting error dim=%d k=%d",dim,k);
            lsq_accum(dim,mat,vec,a,w);  /* accumulate for LSQ normal matrix */
        } /* bp */
    
        radius = 0.5*(sqrt(rad[i-1])+sqrt(rad[i]));
	if (cnt<dim) {
            dprintf(0,"radius %g has %d points: skipped\n",radius,cnt);
            continue;
        }
        lsq_solve(dim,mat,vec,sol);
        printf("%g %d",radius,cnt);
	if (amode) {				/* only print amplitudes */
            for (k=0; k<dim; k++)
                printf(" %g",sol[k]);
        } else {			   /* figure out amp/phase stuff */
            k=0;        /* pointer which 'sol' has been printed */
            if (Qcos[0])		/* if offset wanted, print it now */
                printf(" %g",sol[k++]);
            for( ; k < (dim+1)/2; k++) {	/* go over all amp/phase */
                amp = sqrt(sqr(sol[k]) + sqr(sol[k+dim/2]));
                pha = atan2(sol[k+dim/2],sol[k]) * 180/PI;
                printf(" %g %g",amp,pha);
            }
	}
        printf("\n");
    } /* i */
    return 0;
}

print_header(maxorder,Qcos,Qsin,amode)
int maxorder;
bool Qcos[], Qsin[], amode;
{
    int m;
    
    dprintf(0,"<R> N ");
    if (amode) {
        for(m=0; m<=maxorder; m++)
            if (Qcos[m]) dprintf(0,"A%d ",m);
        for(m=1; m<=maxorder; m++)
            if (Qsin[m]) dprintf(0,"B%d ",m);
    } else {
        if (Qcos[0]) dprintf(0,"C0 ");
        for(m=1; m<=maxorder; m++)
            if (Qsin[m]) dprintf(0,"C%d P%d ",m,m);
    }
    dprintf(0,"\n");
}

