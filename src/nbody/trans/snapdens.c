/*
 *  SNAPDENS: get density estimator for a snapshot, using Kth-nearest
 *            neighbor technique, and also phase space density.
 *
 * See also:   eq (II.2) in Casertano & Hut 1985 ApJ 298, 80.
 *             and references therein.
 *
 *	1-Nov-88	V1.0 created      		PJT
 *     10-Nov-88        V1.1 output correct now         PJT
 *     15-nov-90        V1.2 helpvec                    PJT
 *      9-apr-91        V1.3 adding debug sigma2        PJT
 *     14-apr-91        V1.4 6D search + tfactor        PJT
 *     18-apr-91            a   more general non-equal mass model   PJT
 *	1-apr-01 	    b   compiler warning
 */

#include <stdinc.h>
#include <getparam.h>
#include <math.h>
#include <vectmath.h>		/* otherwise NDIM undefined */
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {
    "in=???\n			  input file name (snapshot) ",
    "out=\n			  optional output file (snapshot) ",
    "kmax=6\n                     number of nearest neighbours ",
    "dens=t\n                     density or phase space density ",
    "tab=f\n                      get an ascii table of all ",
    "format=%e\n                  format for numbers in table ",
    "tfactor=-1.0\n               conversion factor v->r [virial=sqrt(2)]",
    "VERSION=1.4b\n		  1-apr-01 PJT",
    NULL,
};

#define FAC1   4.188790203	/* 3.pi/4 */
#define FAC2   15.74960994	/* (2.pi)^(3/2) */
#define FAC3   1.644934067	/* pi^2/6 */

#define MAXK   128

int   iindex[MAXK+1];              /* pointer to nearest neighbours */
Body  *bindex[MAXK+1];            /* pointer to body */
real  r[MAXK+1];                  /* radius squared to nearest neighbours */

Body  *btab = NULL;               /* pointer to snapshot Body datastructure */
int   nbody, kmax, klen;

bool  Qdens, Qtab;
char  *fmt;
real  tfactor;



nemo_main()
{
    stream instr, outstr;
    real   tsnap, dm;
    string headline=NULL, outfile;
    int i, bits;
    double sqr();
    Body   *bi, *bj;
					
    instr =  stropen(getparam("in"),  "r");
    outfile = getparam("out");
    if (outfile == NULL || *outfile==0)
        outstr = NULL;
    else
        outstr = stropen(outfile, "w");   
    kmax = getiparam("kmax");
    if (kmax > MAXK)
        error("parameter kmax too large");
    Qdens = getbparam("dens"); 
    Qtab = getbparam("tab");  
    fmt = getparam("format");
    tfactor=getdparam("tfactor");
    if (Qdens && tfactor>0)
        warning("tfactor & Qdens incomplete");

/* only do one (the first) snapshot */

    	get_history(instr);
        while (get_tag_ok(instr, HeadlineTag))
        	headline = get_string(instr, HeadlineTag);
        if (!get_tag_ok(instr, SnapShotTag))
		printf ("Input file not a snapshot\n");
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
	if ( (bits & PhaseSpaceBit)==0)
		error("need phasespace in snapshot");
	if ( (bits & MassBit)==0) {
	    warning("no masses in snapshot, assume M=1, m_i=1/%d",nbody);
            dm = 1.0 / (double) nbody;
	    for (i=0, bi=btab; i<nbody; i++, bi++)
                Mass(bi) = dm;
            bits |= MassBit;        /* turn Mass bit on */
        }
        density();   /* compute density for all stars, put result in Aux(bp) */
        if (outstr) {
            put_history(outstr);
            bits |= AuxBit;         /* turn Aux bit on */
            put_snap(outstr,&btab,&nbody,&tsnap,&bits);
        }
                
    strclose(instr);
    if (outstr)
        strclose(outstr);
}


density()
{
    double tmp2, drmin, mtot, m2tot, rdtot, com[NDIM], rmtot[NDIM], mmax;
    Body  *bi, *bj;
    int    i, j, k, kk;
    double raddif();
    
    drmin = HUGE;       /* init minimum interparticle distance */
    mmax = -HUGE;       /* init maximum density */
    rdtot = 0.0;
    for (j=0; j<NDIM; j++)
        rmtot[j] = 0.0;
    mtot = 0.0;
    for (i=0, bi=btab; i<nbody; i++, bi++) {
        klen = 0;                   /* reset nearest neighbours list length */
        for (k=0; k<=kmax; k++) {          /* .. and index pointers etc */
            iindex[k] = -1;
            bindex[k] = NULL;
            r[k] = HUGE;
        }
        for (j=0, bj=btab; j<nbody; j++, bj++) { /* look at all other stars */
            if (i==j)                            /* except itself */
                continue;
            tmp2 = raddif(bi,bj);                /* radius diff squared */
            if (tmp2 < drmin)
                drmin = tmp2;
            if (tmp2 > r[klen])          /* if already larger than largest */
                continue;                        /* goto next star 'j' */
            for (k=0; k<=klen; k++) {     /* check where to insert in list */
                if (tmp2 < r[k]) {
                    if (klen<kmax)               /* increase list length */
                        klen++;
                    for (kk=klen-1; kk>k; kk--) {/* shift higher values */
                        r[kk] = r[kk-1];
                        iindex[kk] = iindex[kk-1];
                        bindex[kk] = bindex[kk-1];
                    } /* for-kk */
                    r[k] = tmp2;                 /* and insert in list */
                    iindex[k] = j;
                    bindex[k] = bj;
                    break;                       /* done with k-loop */
                } /* if */
            } /* for-k */
        } /* for-j */
        stat_nn(bi);                       /* statistics of NN list stars */
        for (j=0; j<NDIM; j++) {
            rmtot[j] += Aux(bi) * Pos(bi)[j]; /* (phase space) density weight */
        }
        mtot +=  Aux(bi);
	m2tot += Aux(bi) * Aux(bi);
        if (Aux(bi) > mmax)
            mmax = Aux(bi);
    } /* for-i */
/* Table header in debug mode */
    dprintf(1,"Weighted_c_o_m[%d]  ",NDIM);
    dprintf(1,"Nearest_neighbor_distance   ");
    dprintf(1,"Max_%s  ",(Qdens ? "dens" : "phase_space_dens"));
    dprintf(1,"R_d  ");
    dprintf(1,"rho_d ");
    dprintf(1,"K_max  \n");
    for (j=0; j<NDIM; j++) {
        com[j] = rmtot[j]/mtot;
        dprintf (0," %f ",com[j]);
    }
    dprintf(0," %f ",sqrt(drmin));
    dprintf (0," %f ",mmax);
    for (i=0, bi=btab; i<nbody; i++, bi++) {
        DISTV(tmp2,Pos(bi),com);
        rdtot += tmp2 * Aux(bi);
    }
    rdtot /= mtot;
    dprintf (0," %f ",rdtot);
    dprintf (0," %f ",m2tot/mtot);
    dprintf (0," %d ",kmax);
    dprintf (0,"\n");
        
}

real raddif(bi,bj)              /* calculate distance squared between 2 stars */
Body *bi, *bj;                  /* saving a sqrt() for now if DISTV was used */
{
    int i, j;
    real rtmp;
    double sqr();
    
    rtmp = 0.0;
    for (i=0; i<NDIM; i++)
        rtmp += sqr(Pos(bi)[i] - Pos(bj)[i]);
    if (tfactor > 0.0) {
        for (i=0; i<NDIM; i++)
            rtmp += sqr((Vel(bi)[i] - Vel(bj)[i])*tfactor);
    }
    
    return(rtmp);
}

/*  stat_nn:   some statistics on the K nearest neighbors of a star
 *
 */
stat_nn(bi)
Body *bi;
{
    real sigma, sigma2, rad, dens, fc, fc2, radius;
    real v1[NDIM], v2[NDIM], s[NDIM];
    int i, k;
    double sqr(), qbe();
    Body *bp;
    
    dens = sigma = sigma2 = 0.0;
    for (i=0; i<NDIM; i++) {
        v1[i] = v2[i] = s[i] = 0.0;
    }
    for (k=0; k<klen; k++) {        /* loop over nearest neighbors */
        bp = bindex[k];
        if (k<klen-1) dens += Mass(bp); /* eq (II.2) in CH 1985 ApJ 298,80) */
        for (i=0; i<NDIM; i++) {
            v1[i] += Vel(bp)[i];
            v2[i] += Vel(bp)[i] * Vel(bp)[i];
        }
    }
    rad = sqrt(r[klen-1]);      /* radius of K-th nearest neighbor */
    dens = dens / (rad*rad*rad*FAC1);        /* space density estimate */
    if (klen!=kmax) {       /* should never occur */
        error ("DENSITY  %f %f %f %f klen=%d\n",rad,dens,0.0,0.0,klen);
        return;
    }
    for (i=0; i<NDIM; i++) {
        s[i] = v2[i]/(double)klen - sqr(v1[i]/(double)klen);
        sigma += s[i];
	sigma2 += v2[i]/(double)klen;
    }
    sigma = sqrt(sigma);
    fc = dens / (sigma*sigma*sigma*FAC2);  /* phase space density estimate */
    fc2 = dens / qbe(rad/tfactor) / FAC3;    /* new 6D phase space estimate */


    Aux(bi) = (Qdens ? dens : fc);        /* replace (phase space) density */
    if (tfactor>0 && !Qdens)  Aux(bi) = fc2;	/* new new */
    if (Qtab) {
        ABSV(radius,Pos(bi));
        printf(fmt,radius);       printf(" ");
        printf(fmt,dens);         printf(" ");
        printf(fmt,fc);           printf(" ");
        printf(fmt,rad);          printf(" ");
        printf(fmt,sigma);        printf (" ");
        printf(fmt,sqrt(sigma2)); printf (" ");	/* debug */
        if (tfactor > 0) printf(fmt,fc2);          
	printf ("\n");

    }
}
