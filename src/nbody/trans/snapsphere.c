/*
 *  SNAPSPHERE:  copy a selected number of particles from an N-body system
 *
 *       8-Apr-87       V1.0 created                            PJT
 *      22-Apr-87       V1.1 headline bug in read_snap          PJT
 *      21-May-87       V1.2 select= keyword                    PJT
 *      18-Jun-87       V1.3 didn't read in 'cs' before         PJT
 *       1-Jul-87       V1.4 didn't check for time              PJT
 *      xx-apr-88       V1.5 readhist() stuff                   PJT
 *       6-jun-88       V1.6 new filestruct                     PJT
 *       1-aug-88       V1.6a bug: mass not copied              PJT
 *      20-dec-88       V1.7 allow also shells to be copied     PJT
 *      14-feb-89       V1.7a buf rrad/vrad ; get_snap          PJT
 *      19-oct-90       V1.8 fixed bug when diagnostics         PJT
 *	 5-may-94	V1.8a fix to keep up with the times	PJT
 *	16-feb-97           b SINGLEPREC fix                    pjt
 *      23-oct-97       V1.9  added central particle select     pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n                     snapshot input file name ",
    "out=???\n                    snapshot output file name ",
    "rcen=0,0,0\n                 pos center ",
    "vcen=0,0,0\n                 vel center ",
    "rrad=0:infinity\n            inner and outer radius in pos",
    "vrad=0:infinity\n            inner and outer radius in vel",
    "times=all\n                  times to select",
    "icen=\n                      overrid center by selecting a particle (1=first)",
    "VERSION=1.9\n                23-oct-97 PJT",
    NULL,
};

string usage="copy a selected number of particles from an N-body system";

#define TIMEFUZZ        0.0001  /* tolerance in time comparisons */

#if !defined(HUGE)
#define HUGE 1e20
#endif

nemo_main()
{
    stream instr, outstr;
    real   tsnap, radsqr, rrad[2], vrad[2];
    vector rcen, vcen;
    string radstr, velstr, times;
    char   *cp;
    Body   *btab = NULL, *b1, *b2;
    int    i, nbody, nbody_out, bits, nrcen, nvcen, icen;
    bool   keep, usebody, scanopt();

    usebody = hasvalue("icen");

    radstr = getparam("rrad");
    cp = strchr(radstr,':');
    if (cp==NULL) {                     /* rrad=rmax (i.e. rmin=0.0) */
        rrad[0] = 0.0;
        rrad[1] = atof(radstr);
    } else {                            /* rrad=rmin:rmax */
       rrad[0] = atof(radstr);
       cp++;            /* point to string after the ':' */
       if (strcmp(cp,"infinity")==0)
          rrad[1] = sqrt(HUGE);
       else
          rrad[1] = atof(cp);
    }
    dprintf (1,"configuration space from r=%f to r=%f\n",rrad[0],rrad[1]);
    rrad[0] = sqr(rrad[0]);
    rrad[1] = sqr(rrad[1]);

    velstr = getparam("vrad");
    cp = strchr(velstr,':');
    if (cp==NULL) {                     /* vrad = vmax i.e. vmin=0.0 */
        vrad[0] = 0.0;
        vrad[1] = atof(velstr);
    } else {                            /* vrad = vmin:vmax */
        vrad[0] = atof(velstr);         
        cp++;                           /* point to string after ':' */
        if (strcmp(cp,"infinity")==0)
           vrad[1] = sqrt(HUGE);
        else                                            
           vrad[1] = atof(cp);
    }
    dprintf (1,"velocity space from r=%f to r=%f\n",vrad[0],vrad[1]);
    vrad[0] = sqr(vrad[0]);
    vrad[1] = sqr(vrad[1]);

    if (usebody) {
        icen = getiparam("icen");
    } else {
        nrcen = nemoinpr(getparam("rcen"), rcen, 3);
        if (nrcen!=NDIM)
            error("rcen= needs %d values", NDIM);
    
        nvcen = nemoinpr(getparam("vcen"),vcen,3);
        if (nvcen!=NDIM)
            error("vcen= needs %d values", NDIM);
    }
        
    times = getparam("times");

                                        /* open files */
    instr = stropen(getparam("in"), "r");
    outstr = stropen(getparam("out"), "w");
    get_history(instr);
    put_history(outstr);                
    
    for (;;) {                  /* do the work in an infinite loop */
        get_history(instr);
        if (!get_tag_ok(instr, SnapShotTag))
                break;                  /* done with work in loop */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if ((bits & PhaseSpaceBit) == 0) /* skip if no coords */
                continue;                /* e.g. when only diagnostics */
        if ((bits & TimeBit) == 0)
                tsnap = 0.0;
        else if (!streq(times,"all") && !within(tsnap, times, TIMEFUZZ))
                continue;
        if (usebody) {
            if (icen < 1 || icen >= nbody)
                error("bad icen=%d ; nbody=%d",icen,nbody);
                b1 = &btab[icen-1];
                SETV(rcen,Pos(b1));
                SETV(vcen,Vel(b1));
        }                
        for (b1=btab, b2=btab; b1 < btab+nbody; b1++) {
            keep = TRUE;
            radsqr=0;
            for (i=0; i<NDIM; i++)
                radsqr += sqr(Pos(b1)[i]-rcen[i]);
            keep = (radsqr>=rrad[0] && radsqr<=rrad[1]);
            if (keep) {
                radsqr=0;
                for (i=0; i<NDIM; i++)
                    radsqr += sqr(Pos(b1)[i]-rcen[i]);
                keep = (radsqr>=vrad[0] && radsqr<=vrad[1]);
            }
            if (!keep)
                continue;           /* don't copy */
            if (b1==b2) {           /* no need to cp, still in sync */
                b2++;
                continue;
            }
            bcopy (b1, b2, sizeof(Body));
            b2++;
        }
        nbody_out = b2 - btab;
        dprintf (1,"Snapshot: time=%f copied %d particles\n",tsnap,nbody_out);
        if (nbody_out > 0)
            put_snap(outstr, &btab, &nbody_out, &tsnap, &bits);
        else
            warning("Time = %f no bodies to copy",tsnap);
    }
    strclose (instr);
    strclose (outstr);
}
