/*
 * ROTCURVES:	draw rotation curves along a major axis (x,y,z) of a 
 *		composite potential(NEMO5) model. 
 *		An input table is assumed to be a data file with
 *		in column 1 radius and column 2 velocity, in the same
 *		units as in program. For option=vel or ome the data is
 *		also overlayed on the composite model, and some
 *		statistic is used to estimate 'fit' qualitity
 *
 *		V1.0 25-oct-90 	derived from potlist
 *		V1.1 28-oct-90  What the heck, have optional Lindblad too
 *              V1.2  2-nov-90  Table of Vrad(r_from_sun)
 *		V1.2a 6-nov-90  input table with data
 *		    b 8-nov-90  Show chi^2 between input table and models
 *              V1.3 15-nov-90  Using new get_atable() table interface
 *                  a17-nov-90  Allow input of errors in rotcurv in= too
 *		    b 7-mar-92  Happy gcc2.0  				PJT
 *	       V1.4  11-aug-92  merged two divergent versions		PJT
 *             V1.5  12-nov-93  allow other resonances (n=2 being default) PJT
 *	       V1.5a 26-mar-95  proto
 *		   b  8-feb-96  bigger tables
 *		   c 20-feb-97  fixed for SINGLEPREC
 *	       V1.6   8-apr-97  find resonances, more SINGLEPREC needed PJT
 *                 a 19-feb-01  bad bug in mode=LV, used the wrong radius
 *      	   b 13-sep-01  better prototype for proc
 *             V1.7   7-jul-02  added format=
 *                 a  9-jan-03  more C++ friendly
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <yapp.h>
#include <axis.h>
#include <spline.h>
#include <potential.h>

string defv[] = {
    "name1=???\n         Name of potential_1",
    "pars1=\n            Parameters for potential_1 (1st one is pattern speed)",
    "file1=\n            Optional data file associated with potential_1",
    "name2=\n            --same for potential_2",
    "pars2=\n            --same for potential_2",
    "file2=\n            --same for potential_2",
    "name3=\n            --same for potential_3",
    "pars3=\n            --same for potential_3",
    "file3=\n            --same for potential_3",
    "name4=\n            --same for potential_4",
    "pars4=\n            --same for potential_4",
    "file4=\n            --same for potential_4",
    "radii=0:2:0.1\n     Radii to sample rotation curve",
    "axis=x\n            Axis along which to sample rotation curve",
    "mode=velocity\n     Output mode: {velocity, omega, lv}",
    "n=2\n               Resonance (Omega +/- Kappa/n) to plot if mode=omega",
    "r0l=1,90,0:1:0.1\n  Solar radius and longitude in lv-mode; sample radii",
    "plot=t\n            Make Plot (t|f)?",
    "tab=f\n             Make table (t|f)?",
    "format=%f\n         Format for table output",
    "xrange=\n           X-range for plot (radius)",
    "yrange=\n           Y-range for plot (velocity/omega)",
    "headline=\n         Optional plot label for identification",
    "in=\n               Optional input rotation curve table",
    "cols=1,2\n          Columns for r, v, dr, dv (use 0 when not present)",
    "VERSION=1.7a\n      9-jan-03 PJT",
    NULL,
};

string usage = "display rotation curve of composite potentials";

#ifndef MAXPT
#define MAXPT 10000
#endif

potproc_double mypot1, mypot2, mypot3, mypot4;

real xplot[2], yplot[2];
char plotmsg[256];


local real xtrans(real), ytrans(real);
local void lindblad(int, real *, real *, real *, real *, real *, real *, int);
local void goodness(int ,real *, real *, int, real *, real *, real *);
local int  read_table(stream, int, real *, real *, real *, real *, int *);
local void lv(int nrad, real *, real *, real, real, int, real *);
local int  read_table(stream, int, real *, real *, real *, real *,int *);
local void peak(int, real *, real *, int, int, real *, real *);

extern int get_atable(stream ,int, int *, real **, int);
extern void lsq_zero(int n, real *mat, real *vec);
extern void lsq_accum(int n, real *mat, real *vec, real *a, real w);
extern void lsq_solve(int n, real *mat, real *vec, real *sol);


void nemo_main()
{
    int    i, dir, nrad, npots=0, ltype, ndim = NDIM, nx, ny, ns, ndat, nret;
    int    cols[4], n, idx, idx_max;
    real   pmax, symsize, rr, omk_max = 0.0, omk_rmax;
    real   rad[MAXPT], *vel, *vel1, *vel2, *vel3, *vel4, *curve;
    real   *ome, *kap, *opk, *omk, r0l[MAXPT+2], omega, *f;
    real   inrad[MAXPT], invel[MAXPT], inrade[MAXPT], invele[MAXPT];
    double pos[3], acc[3], pot, time = 0.0;
/*    char   *fmt, s[20], pfmt[256];    */
    char   headline[256], fmt1[80];
    string axis, mode, infile, plotlabel;
    stream instr;
    bool   Qtab, Qplot, Qome, Qvel, Qlv, Qin, QoILR;

    mode = getparam("mode");
    n = getiparam("n");
    plotlabel = getparam("headline");
    sprintf(fmt1,"%s ",getparam("format"));
    Qome = (*mode == 'o');      /*  options are: velocity|omega|lv */
    Qlv = (*mode == 'l');
    Qvel = (*mode == 'v');
    Qtab = getbparam("tab");
    Qplot = getbparam("plot");
    infile = getparam("in");
    Qin =  (*infile != 0);
    if (Qin) {
        nret = nemoinpi(getparam("cols"),cols,4);
        if (nret<0 || nret > 4) error("cols= requires 4 numbers");
        for (i=nret; i<4; i++)
            cols[i] = 0;
        instr = stropen(infile,"r");
        ndat = read_table(instr,MAXPT,inrad,invel,inrade,invele,cols);
        strclose(instr);
    }
    
    mypot1 = get_potential(getparam("name1"),getparam("pars1"),getparam("file1"));
    omega = get_pattern();
    dprintf(0,"Pattern speed: %f\n",omega);
    mypot2 = get_potential(getparam("name2"),getparam("pars2"),getparam("file2"));
    mypot3 = get_potential(getparam("name3"),getparam("pars3"),getparam("file3"));
    mypot4 = get_potential(getparam("name4"),getparam("pars4"),getparam("file4"));
    headline[0] = '\0';         /* accumulate headline */
    if (mypot1) {
        strcat(headline,getparam("name1"));
        strcat(headline,"(");
        strcat(headline,getparam("pars1"));
        strcat(headline,")");
        npots++;
    } 
    if (mypot2) {
        strcat(headline,getparam("name2"));
        strcat(headline,"(");
        strcat(headline,getparam("pars2"));
        strcat(headline,") ");
        npots++;
    }
    if (mypot3) {
        strcat(headline,getparam("name3"));
        strcat(headline,"(");
        strcat(headline,getparam("pars3"));
        strcat(headline,") ");
        npots++;
    }
    if (mypot4) {
        strcat(headline,getparam("name4"));
        strcat(headline,"(");
        strcat(headline,getparam("pars4"));
        strcat(headline,")");
        npots++;
    }

    nrad = nemoinpr(getparam("radii"),rad,MAXPT);   /* get radii */
    if (nrad <= 0)
        warning("Using %d radii is not very productive",nrad);
    vel  = (real *) allocate(sizeof(real) * nrad);  /* allocate stuff */
    vel1 = (real *) allocate(sizeof(real) * nrad);
    vel2 = (real *) allocate(sizeof(real) * nrad);
    vel3 = (real *) allocate(sizeof(real) * nrad);
    vel4 = (real *) allocate(sizeof(real) * nrad);
    if (Qome) {
        ome = (real *) allocate(4 * sizeof(real) * nrad);  /* plus spline */
        kap = (real *) allocate(sizeof(real) * nrad);
        opk = (real *) allocate(sizeof(real) * nrad);
        omk = (real *) allocate(sizeof(real) * nrad);
    } 

    axis = getparam("axis");
    dir = 0;
    if (*axis == 'x') dir=0;
    if (*axis == 'y') dir=1;
    if (*axis == 'z') dir=2;
    if (dir>NDIM) error("Axis %s not supported in NDIM=%d",axis,NDIM);

    pmax = 0.0;

    for (i=0; i<nrad; i++) {            /* loop to compute */
        CLRV(pos);                      /* clear positions */
        pos[dir] = rad[i];              /* set the right axis */
        vel[i] = 0.0;
        if (mypot1) {
            CLRV(acc);
            (*mypot1) (&ndim,pos,acc,&pot,&time);
            vel1[i] = -rad[i] * acc[dir];
            vel[i] += vel1[i];
            vel1[i] = sqrt(vel1[i]);        
        }
        if (mypot2) {
            CLRV(acc);
            (*mypot2) (&ndim,pos,acc,&pot,&time);
            vel2[i] = -rad[i] * acc[dir];
            vel[i] += vel2[i];
	    vel2[i] = sqrt(vel2[i]);        
        }
        if (mypot3) {
            CLRV(acc);
            (*mypot3) (&ndim,pos,acc,&pot,&time);
            vel3[i] = -rad[i] * acc[dir];
            vel[i] += vel3[i];
	    vel3[i] = sqrt(vel3[i]);        
        }
        if (mypot4) {
            CLRV(acc);
            (*mypot4) (&ndim,pos,acc,&pot,&time);
            vel4[i] = -rad[i] * acc[dir];
            vel[i] += vel4[i];
            vel4[i] = sqrt(vel4[i]);        
        }
        vel[i]  = sqrt(vel[i]);        
    }
    if (Qome) {
	lindblad(nrad,rad,vel,ome,kap,opk,omk,n);
        if (omega> 0.0) {                               /* compute resonances */
            f = opk;
            idx = nrad-1;
            if (omega < f[idx]) {
                warning("Radii not far enough out for OLR: %g",f[idx]);
                f = ome;
                if (omega < f[idx]) {
                    warning("Radii not far enough out for CR: %g",f[idx]);
                    f = omk;
                }
            }
            QoILR = FALSE;
            for(; idx>0; idx--) {
                if (omk[idx] > omk_max) {
                    idx_max = idx;
                    omk_max = omk[idx];
                }
                if (f==omk) {
                    if (QoILR) {
                        if (omega < f[idx]) continue;
                    } else {
                        if (omega > f[idx]) continue;
                    }
                } else {
                    if (omega > f[idx]) continue;
                }
                
                /* found a resonance: */

                rr = rad[idx] + (rad[idx+1]-rad[idx])*
                                (omega-f[idx])/(f[idx+1]-f[idx]);
                if (f == omk) {
#if 0                    
                    if (QoILR) {
                        dprintf(0,"iILR: %g\n",rr);
                        break;
                    } else {
                        dprintf(0,"oILR: %g\n",rr);
                        QoILR = TRUE;
                    }
#endif                    
                } else if (f == ome) {
                    dprintf(0,"CR: %g\n",rr);
                    f = omk;
                } else if (f == opk) {
                    dprintf(0,"OLR: %g\n",rr);
                    f = ome;
                } else
                    error("impossble resonance");
            }
            peak(nrad,rad,omk,idx_max,1, &omk_rmax, &omk_max);
            dprintf(0,"OMK_max: %g\n",omk_max);
            dprintf(0,"OMK_rmax: %g\n",omk_rmax);

            if (omega < omk_max) {			/* search for ILR */
            	for (idx=idx_max; idx<nrad; idx++) {
                    if (omega > omk[idx]) {
                        rr = rad[idx-1] + (rad[idx]-rad[idx-1])*
                                (omega-f[idx-1])/(f[idx]-f[idx-1]);
                        dprintf(0,"oILR: %g\n",rr);
                        break;
                    }
            	}
                for (idx=idx_max; idx>0; idx--) {
                    if (omega > omk[idx]) {
                        rr = rad[idx] + (rad[idx+1]-rad[idx])*
                               (omega-f[idx])/(f[idx+1]-f[idx]);
                        dprintf(0,"iILR: %g\n",rr);
                        break;
                    }
            	}
            }
        }
    }
    for (i=0; i<nrad; i++) {                            /* loop to print */
        if (Qtab) {
	  printf(fmt1,rad[i]);
	  printf(fmt1,vel[i]);
	}
	if (Qtab && npots>1 && !Qome) {
	    if (mypot1) printf(fmt1,vel1[i]);
	    if (mypot2) printf(fmt1,vel2[i]);
	    if (mypot3) printf(fmt1,vel3[i]);
	    if (mypot4) printf(fmt1,vel4[i]);
        }
        if (Qtab && Qome) {
	  printf(fmt1,ome[i]);
	  printf(fmt1,kap[i]);
	  printf(fmt1,opk[i]);
	  printf(fmt1,omk[i]);
	}
	if (Qtab) printf("\n");
        if (Qome)
            pmax = MAX(pmax,opk[i]);
        else
            pmax = MAX(pmax,vel[i]);
    }
    if (Qin && Qvel) 
        goodness(nrad,rad,vel,ndat,inrad,invel,(cols[3]>0?invele:NULL));
    if (Qplot) {
        plinit("***",0.0,20.0,0.0,20.0);                /* open device */
        nx = nemoinpr(getparam("xrange"),xplot,2);      /* get xrange in plot */
        switch(nx) {
         case 0:
            xplot[0] = rad[0];
         case 1:
            xplot[1] = rad[nrad-1];
            break;
         case 2:
            break;
         default:
            warning("xrange= only accepts two values");
            break;
        }
        ny = nemoinpr(getparam("yrange"),yplot,2);      /* get yrange in plot */
        switch(ny) {
         case 0:
            yplot[0] = 0.0;
            yplot[1] = 1.1 * pmax;      /* extra 10% for egde */
            break;
         case 1:
            yplot[1] = 1.1 * pmax;      /* extra 10% for egde */
            break;
         case 2:
            break;
         default:
            warning("yrange= only accepts two values");
            break;
        }
        xaxis ( 2.0, 2.0, 16.0, xplot, -7, xtrans, "R");    /* plot axes */
        xaxis ( 2.0,18.0, 16.0, xplot, -7, xtrans, NULL);
        if (Qome)
            yaxis ( 2.0, 2.0, 16.0, yplot, -7, ytrans, "[V/R]");
        else
            yaxis ( 2.0, 2.0, 16.0, yplot, -7, ytrans, "V");
        yaxis (18.0, 2.0, 16.0, yplot, -7, ytrans, NULL);
        if (*plotlabel)
            pltext(plotlabel,2.0,18.5,0.5,0.0);
        else
            pltext(headline,2.0,18.5,0.35,0.0);
        if (*plotmsg)
            pltext(plotmsg,8.0,2.5,0.25,0.0);

        curve = (Qome ? ome : vel);            /* assign first curve */
        plltype(3,1);                                 /* thick solid line */
        plmove(xtrans(rad[0]),ytrans(curve[0]));
        for (i=1; i<nrad; i++)
            plline(xtrans(rad[i]),ytrans(curve[i]));
        if (Qome) {                   /* if Lindblad - plot omk, opk */
            plltype(1,1);                      /* all regular solid lines */
            plmove(xtrans(rad[0]), ytrans(omk[0]));
            for (i=1; i<nrad; i++)
                plline(xtrans(rad[i]),ytrans(omk[i]));
            plmove(xtrans(rad[0]), ytrans(opk[0]));
            for (i=1; i<nrad; i++)
                plline(xtrans(rad[i]),ytrans(opk[i]));
        } else if (npots>1) {            /* if velocity and > 1 component */
            ltype = 1;
            if (mypot1) {
                plltype(1,++ltype);
                plmove(xtrans(rad[0]),ytrans(vel1[0]));
                for (i=1; i<nrad; i++)
                    plline(xtrans(rad[i]),ytrans(vel1[i]));
            }
            if (mypot2) {
                plltype(1,++ltype);
                plmove(xtrans(rad[0]),ytrans(vel2[0]));
                for (i=1; i<nrad; i++)
                    plline(xtrans(rad[i]),ytrans(vel2[i]));
            }
            if (mypot3) {
                plltype(1,++ltype);
                plmove(xtrans(rad[0]),ytrans(vel2[0]));
                for (i=1; i<nrad; i++)
                    plline(xtrans(rad[i]),ytrans(vel3[i]));
            }
            if (mypot4) {
                plltype(1,++ltype);
                plmove(xtrans(rad[0]),ytrans(vel2[0]));
                for (i=1; i<nrad; i++)
                    plline(xtrans(rad[i]),ytrans(vel4[i]));
            }
        }
	plltype(1,1); 
        symsize = 0.1;
        if (Qin && Qvel) {           /* if input file with velocities */
            for (i=0; i<ndat; i++)
                plbox(xtrans(inrad[i]),ytrans(invel[i]),symsize);
            if (cols[3]>0) {        /* if error bars in radius */
                for (i=0; i<ndat; i++) {
                    plmove(xtrans(inrad[i]-inrade[i]),ytrans(invel[i]));
                    plline(xtrans(inrad[i]+inrade[i]),ytrans(invel[i]));
                }
            }
            if (cols[4]>0) {        /* if error bars in velocity */
                for (i=0; i<ndat; i++) {
                    plmove(xtrans(inrad[i]),ytrans(invel[i]-invele[i]));
                    plline(xtrans(inrad[i]),ytrans(invel[i]+invele[i]));
                }
            }
        } else if (Qin && Qome) {       /* if input file with omega */
            for (i=0; i<ndat; i++)
                plbox(xtrans(inrad[i]),ytrans(invel[i]/inrad[i]),symsize);
        }
        plstop();
    }  /* if plot vel/ome */
    if (Qlv) {
        ns = nemoinpr(getparam("r0l"),r0l,MAXPT+2) - 2;
        if (ns < 0)
            error("r0l= needs at least two values: r0 and l");
        else if (ns==0)
            warning("r0l= no lv-radii array supplied");
        lv(nrad,rad,vel,r0l[0],r0l[1],ns,&r0l[2]);
    }
}

real xtrans(real x)
{
    return  2.0 + 16.0*(x-xplot[0])/(xplot[1]-xplot[0]);
}
	 
real ytrans(real y)
{
    return  2.0 + 16.0*(y-yplot[0])/(yplot[1]-yplot[0]);
}
		   

void lindblad(int nrad, real *rad, real *vel, 
                        real *ome, real *kap, real *opk, real *omk, int n)
{
    int i, nbad=0;
    real fac;

    fac = 1.0/(real)n;

    for (i=0; i<nrad; i++) {        /* get ome */
        if (rad[i]==0.0) continue;
        ome[i] = vel[i]/rad[i];
        if (i>0 && rad[i] <= rad[i-1])
            error("lindblad: Need sorted data for spline fit");
    }
    for (i=0; i<nrad; i++) {        /* fix up when rad=0 && set ome^2 */
        if (rad[i]==0.0)
            if (i==0)
                ome[0] = ome[1];
            else
                error("Cannot handle zero radius here");
        else
            ome[i] = sqr(ome[i]);
    }
    spline (&ome[nrad], rad, ome, nrad);    /* spline through ome^2 */
    for (i=0; i<nrad; i++) {
        kap[i] = rad[i] * spldif(rad[i],rad,ome,&ome[nrad],nrad) + 4*ome[i];
        if (kap[i] < 0) {
            nbad++;
            dprintf(1,"r=%f v=%f kappa^2 = %f\n",rad[i],vel[i],kap[i]);
            kap[i] = -sqrt(-kap[i]);
        }
            kap[i] = sqrt(kap[i]);
    }
    for (i=0; i<nrad; i++) {
        ome[i] = sqrt(ome[i]);
        opk[i] = ome[i] + fac * kap[i];
        omk[i] = ome[i] - fac * kap[i];
    }
    if (nbad>0) warning("There were %d radii with badly formed kappa",nbad);
}

void lv(int nrad, real *rad, real *vel, real r0, real l, int ns, real *rs)
{
    int i;
    real *vel3=NULL, sinl, cosl, v0, r, v, r1;
    
    vel3 = (real *) allocate(3*nrad*sizeof(real));

    spline (vel3, rad, vel, nrad);    /* spline through vel(rad) */
    v0 = seval(r0,rad,vel,vel3,nrad);   /* get solar velocity */
    dprintf(1,"Solar radius=%g velocity=%g\n",r0,v0);
    sinl = sin(l/57.2958);
    cosl = cos(l/57.2958);
    for (i=0; i<ns; i++) {
        r1=rs[i];
        if (r1==0.0) 
            printf("%f %f\n",r1,r1);
        else {
            r = sqrt(sqr(r1)+sqr(r0)-2*r0*r1*cosl);     /* radius from GC */
            if (r>rad[nrad-1]) 
                warning("Insufficient points in rotation curve, r=%g > %g",
				r,rad[nrad-1]);
            v = seval(r,rad,vel,vel3,nrad);
            printf("%f %f\n",r1,(v*r0/r-v0)*sinl);
        }
    }
    free(vel3);
} /* lv */


local int read_table(stream instr, int nmax, real *x, real *y, 
                     real *dx, real *dy,
                     int cols[4])
{
    int  colnr[4], nread = 0, ncols=0;
    real *coldat[4];

    if (cols[0]>0) {
        dprintf(1,"Reading radii from column %d\n",cols[0]);
        colnr[ncols] = cols[0];
        coldat[ncols] = x;
        ncols++;
    } else
        error("read_table: Radii not requested");

    if (cols[1]>0) {
        dprintf(1,"Reading velocity from column %d\n",cols[1]);
        colnr[ncols] = cols[1];
        coldat[ncols] = y;
        ncols++;
    } else
        error("read_table: velocities not requested");

    if (cols[2]>0) {
        dprintf(1,"Reading radii errors from column %d\n",cols[2]);
        colnr[ncols] = cols[2];
        coldat[ncols] = dx;
        ncols++;
    } 

    if (cols[3]>0) {
        dprintf(1,"Reading velocity errors from column %d\n",cols[3]);
        colnr[ncols] = cols[3];
        coldat[ncols] = dy;
        ncols++;
    } 

    dprintf(2,"Reading table: \n");
    nread = get_atable(instr,ncols,colnr,coldat,nmax);
    dprintf(2,"Read %d entries from table\n",nread);
    return nread;
}


local void goodness(int n,real *x, real *y, int m, 
                    real *xdat, real *ydat, real *dydat)
{
    real  *x3, chi2, e, o, s;
    int    i, k=0;

    x3 = (real *) allocate(3*n*sizeof(real));
    spline(x3,x,y,n);           /* get spline coefs for interpolation */
    chi2 = 0.0;
    if (dydat==NULL) {
        for (i=0; i<m; i++) {
            o = ydat[i];
            e = seval(xdat[i],x,y,x3,n);
            if (e==0.0) {
                warning("Skipping datapoint %d (%g %g)",i+1,xdat[i],ydat[i]);
                continue;
            }
            chi2 += sqr(e-o)/e;
            k++;
        }
        sprintf(plotmsg,"For %d data points chi^2 = %g\n",k,chi2);
    } else {
        for (i=0; i<m; i++) {
            o = ydat[i];
            e = seval(xdat[i],x,y,x3,n);
            s = dydat[i];
            if (s==0.0) {
                warning("Skipping datapoint %d (%g %g %g)",
                                i+1,xdat[i],ydat[i],dydat[i]);
                continue;
            }
            chi2 += sqr((e-o)/s);
            k++;
        }
        chi2 /= k;      /* degrees of freedom is n-#-1 though */
        if (chi2<1)
            sprintf(plotmsg,"N=%d chi^2/N = %g (Too good a fit?)",k,chi2);
        else if (chi2<10)
            sprintf(plotmsg,"N=%d chi^2/N = %g!",k,chi2);
        else if (chi2<100)
            sprintf(plotmsg,"N=%d chi^2/N = %g (Bad fit I'd say)",k,chi2);
	else
            sprintf(plotmsg,"N=%d chi^2/N = %g (You can't be serious)",k,chi2);
    }
    printf("%s\n",plotmsg);
}




local void peak(int n, real *x, real *y, int idx, int mode,
               real *xmax, real *ymax)
{
    int i, k, order = 2;
    real mat[9], vec[3], sol[3], a[4];

    if (idx < 1 || idx > n-2) error("peak: too close to the edge");
    if (ABS(mode) != 1)error("peak: mode should be 1 or -1");

    lsq_zero(order+1, mat, vec);
    for (i=idx-1; i<=idx+1; i++) {
        a[0] = 1.0;
        for (k=0; k<order; k++)
            a[k+1] = a[k] * (x[i]-x[idx]);
        a[order+1] = mode * y[i];
        lsq_accum(order+1,mat,vec,a,1.0);
    }
    lsq_solve(order+1,mat,vec,sol);
    dprintf(1,"Peak (poly2) fit near idx=%d (%g,%g)\n",idx+1,x[idx],y[idx]);
    *xmax = x[idx] - sol[1]/(2*sol[2]);
    *ymax = sol[0]+sqr(sol[1])/(2*sol[2]);
    
}
