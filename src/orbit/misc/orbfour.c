/*
 *  orbfour: get some fourier coeficients from an orbit
 *
 *      1-sep-90       V1.0 created     - only even cos-terms   PJT
 *      19-dec-90      V1.1 similar method to snapfour - quick hack PJT
 *	 7-mar-92 	gcc happy				PJT
 *	24-may-92      V1.1c added <potential.h>		PJT
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <potential.h>
#include <orbit.h>

string defv[] = {
    "in=???\n               Input filename (an orbit)",
    "cos=0:4:1\n           List of noj-zero cos(m.phi) terms",
    "sin=1:4:1\n           List of non-zero sin(m.phi) terms",
    "xvar=x\n              X variable (**cannot change yet**)",
    "yvar=y\n              Y variable (**cannot change yet**)",
    "fvar=vy\n             Fourier Observable to be decomposed (vx,vy,vr,vt)",
    "amode=t\n             Display amps or amp/phase if possible?",
    "VERSION=1.1c\n        24-may-92 PJT",
    NULL,
};

string usage = "get some fourier coeficients from an orbit";

string  infile;                 /* file names */
stream  instr;                  /* file streams */

orbitptr o_in=NULL;                     /* pointer to input orbit */

#define MAXORDER  32
int maxorder;

void orb_four(), print_header();


void nemo_main ()
{
    bool Qcos[MAXORDER+1], Qsin[MAXORDER+1], amode;
    int    i, n, maxorder, tmpi[MAXORDER+1];
    string xvar, yvar, fvar;
    rproc fproc;
    real otr_vr(), otr_vt(), otr_vy(), otr_vx();

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

    fvar = getparam("fvar");
    if (streq(fvar,"vr"))
        fproc=otr_vr;
    else if (streq(fvar,"vt"))
        fproc=otr_vt;
    else if (streq(fvar,"vx"))
        fproc=otr_vx;
    else if (streq(fvar,"vy"))
        fproc=otr_vy;
    else
        error("fvar=%s not available for orbit transformation",fvar);
    amode = getbparam("amode");

    xvar = getparam("xvar");
    if (!streq(xvar,"x")) error("Cannot handle this xvar yet");
    yvar = getparam("yvar");
    if (!streq(yvar,"y")) error("Cannot handle this yvar yet");


    instr = stropen(getparam("in"),"r");
    while (read_orbit(instr,&o_in))
    	orb_four(o_in, maxorder, fproc, Qcos, Qsin, amode);
    strclose(instr);
}

void orb_four(optr, maxorder, fproc, Qcos, Qsin, amode)
orbitptr optr;
int maxorder;
rproc fproc;
bool Qcos[], Qsin[], amode;
{
    real   th, r2, amp, pha, radius;
    int    k,m,cnt,dim,ip;
    real mat[2*(MAXORDER+1)*(MAXORDER+1)],vec[2*(MAXORDER+1)];
    real sol[2*(MAXORDER+1)], a[2*(MAXORDER+1)+1];
    void lsq_zero(), lsq_accum(), lsq_solve();
    permanent bool first=TRUE;

    dprintf(1,"[Analyzing %d steps\n",Nsteps(optr));


    for (m=0, dim=0; m<=maxorder; m++)  /* count dimension of matrix needed */
        if (Qcos[m]) dim++;
    for (m=1; m<=maxorder; m++)
        if (Qsin[m]) dim++;
    if (dim==0) {
       warning("orb_four: dim=0");
       return;
    }
    dprintf(1,"orb_four: maxorder=%d dim=%d\n",maxorder,dim);
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

        lsq_zero(dim,mat,vec);          /* reset accum. matrix and vector */
        cnt = 0;                        /* count points in this ring */
        for (ip=0; ip<Nsteps(optr); ip++) {  /* loop for all steps */
            cnt++;
            th = atan2(Yorb(optr,ip), Xorb(optr,ip));
            k=0;                            /* always re-count how many coefs */
            for(m=0; m<=maxorder; m++) {     /* cos(m.theta) */
                if (!Qcos[m]) continue;
                a[k++] = cos(m*th);
            }
            for(m=1; m<=maxorder; m++) {     /* sin(m.theta) */
                if (!Qsin[m]) continue;
                a[k++] = sin(m*th);
            }
            a[k] = fproc(Xorb(optr,ip),Yorb(optr,ip),Zorb(optr,ip),
                         Uorb(optr,ip),Vorb(optr,ip),Worb(optr,ip));
            dprintf(1,"adding %d: r^2=%g th=%g, vy=%g\n",
                              cnt,r2,th*180/PI,a[k]);
            if (k!=dim) error("orb_four: Counting error dim=%d k=%d",dim,k);
            lsq_accum(dim,mat,vec,a,1.0);  /* accumulate for LSQ normal matrix */
        } /* ip */
    
        if (cnt<dim) {
            dprintf(0,"orbit has %d points: skipped\n",cnt);
            return;
        }
        lsq_solve(dim,mat,vec,sol);
        radius = MAX(Xorb(optr,0),Yorb(optr,0));  /* assign arbitrary radius */
        printf("%g %d",radius,cnt);
        if (amode) {                            /* only print amplitudes */
            for (k=0; k<dim; k++)
                printf(" %g",sol[k]);
        } else {                           /* figure out amp/phase stuff */
            k=0;        /* pointer which 'sol' has been printed */
            if (Qcos[0])                /* if offset wanted, print it now */
                printf(" %g",sol[k++]);
            for( ; k < (dim+1)/2; k++) {        /* go over all amp/phase */
                amp = sqrt(sqr(sol[k]) + sqr(sol[k+dim/2]));
                pha = atan2(sol[k+dim/2],sol[k]) * 180/PI;
                printf(" %g %g",amp,pha);
            }
        }
        printf("\n");

}

void print_header(maxorder,Qcos,Qsin,amode)
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

/*  functions for orbit transformations */
real otr_vx(x,y,z,vx,vy,vz)
real x,y,z,vx,vy,vz;
{
    return(vx);

}

real otr_vy(x,y,z,vx,vy,vz)
real x,y,z,vx,vy,vz;
{
    return(vy);
}

real otr_vr(x,y,z,vx,vy,vz)
real x,y,z,vx,vy,vz;
{
    return( (x*vx+y*vy)/sqrt(x*x+y*y) );
}

real otr_vt(x,y,z,vx,vy,vz)
real x,y,z,vx,vy,vz;
{
    return( sqrt(vx*vx+vy*vy -  sqr(x*vx+y*vy)/(x*x+y*y)) );
}

