/*
 *  UNBIND:     iteratively finds unbound stars to an N-body system
 *
 *      12-nov-87       V1.0 created                            PJT
 *      15-nov-87       V1.2 output and ecutoff added           PJT
 *      19-nov-87       V1.3 optional output unbound stars      PJT
 *       8-mar-88       V1.4 LSM implementation + hist-data     PJT
 *      xx-mar-88       V1.6 enhanced options (iter,map)        PJT
 *       6-jun-88       V1.7 new filestruct (should be V2.0)    PJT
 *      24-oct-88       V1.8 added output of key field as star# PJT
 *      10-may-89       V1.9 kept potentials in there if they were      PJT
 *      12-feb-90       V2.0 multi-time snapshots  PJT
 *      15-nov-90       V2.1 NEMO 2.x stuff     PJT
 *      14-aug-91       V2.2a     renamed time->stime  PJT
 *      21-sep-91       V2.3 fixed bug not writing key;         PJT
 *	 7-jan-92	V2.3a don't output & warn when no stars left 	PJT
 *      30-mar-92            fixed typo: time= should have been times=  PJT
 *                           but was never used
 *	22-dec-92	V2.4 again write out 0 length snapshots	PJT
 *      28-dec-92       V2.4a - fixed cases where Mass output negative  PJT/SF
 *	15-aug-96       V2.5 code cleaned (old version crashed on linux)  PJT
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
    "in=???\n           Input file name",
    "out=???\n          Output file name",
    "exact=f\n          Exact N-squared potential ?",
    "eps=0.025\n        Softening length in case exact potentials",
    "ecutoff=0.0\n      Cutoff for (un)binding",
    "bind=t\n           Output bound(t) or unbound(f) stars",
    "map=f\n            Print map of bound/unbound",
    "times=all\n        Times of shapshots to copy",
    "VERSION=2.5\n      15-aug-96 PJT",
    NULL,
};

string usage="iteratively finds unbound stars to an N-body system";

local string headline;                /* random text message */
local string iname,oname;             /* file  names */
local stream instr,outstr;            /* file streams */

local int    nbody;                   /* number of particles in input snapshot */
local double stime;                   /* current time in snapshot */
local Body   *btab;                   /* pointer to input snapshot */
local int    bits;                    /* bitmap of output quantities */
 
local double  ecutoff;                /* cutoff energy */
local int     nesc;                   /* counter how many flagged as escaped */

local double sqreps;                  /* square of softening length */
local bool   Qexact;                  /* exact potential ? */
local bool   Qbind;                   /* true=keep bound   false=keep escapers */
local bool   Qmap;                    /* true=make map of bound/unnound */


nemo_main()
{
    int  i;
    double ekin;
    Body  *bp;

    instr = stropen(getparam("in"), "r");       /* get parameters */
    outstr = stropen(getparam("out"),"w");
    sqreps = sqr(getdparam("eps"));
    ecutoff = getdparam("ecutoff");
    Qexact = getbparam("exact");
    Qbind = getbparam("bind");
    Qmap = getbparam("map");
    dprintf (1,"Stars with binding energy above %f will be ",ecutoff);
    if (Qbind)
        dprintf(1,"removed\n");
    else
        dprintf(1,"kept\n");

    while (read_snap()) {             /* read snapshot */
        nesc = 0;
        for (bp=btab; bp<btab+nbody; bp++) {
                ekin = 0;
                for (i=0; i<NDIM; i++)
                        ekin += sqr(Vel(bp)[i]);
                ekin *= 0.5;
                Phi(bp) += ekin;
                if (Phi(bp) >= ecutoff) {
                        Mass(bp) *= -1;         /* flag as escaper */
                        nesc++;
                } 
        }
        if (Qbind)
            dprintf (0,"%d out of %d stars bound and written to file\n",
                            nbody-nesc,nbody);
        else
            dprintf (0,"%d out of %d stars unbound and written to file\n",
                            nesc,nbody);
        if (Qmap)
           map();
        write_snap();                   /* write out (un)bound stars */
    } 
    strclose(instr);
}


read_snap()
{                               
    int    i;
    Body   *b;
        
    for(;;) {         /* 'infinite' loop until we find a ParticleTag */
        get_history(instr);
        if (!get_tag_ok(instr,SnapShotTag))     /* we DO need a SnapShot !! */
                return(0);
        get_snap(instr, &btab, &nbody, &stime, &bits);
        if ((bits & MassBit) == 0 || (bits & PhaseSpaceBit) == 0)
                error("missing essential data");
        if (Qexact) {
            dprintf (0,"Doing an exact potential calculation\n");
            exact();            /* fill in newtonian potentials */
        } else if ((bits & PotentialBit)==0)            
            error("missing potentials in snapshot, use hackforce or exact=t");
        else
           dprintf (1,"Using potentials in snapshot for energy calculation\n");
        if ((bits & KeyBit) == 0) {
            warning ("Keys (re)set according to their order in file");
            for (i=0, b=btab; i<nbody; i++, b++)
                Key(b) = i;
            bits |= KeyBit;     /* turn KeyBit on, since it was off */
        }
        return(1);
    } /* for(;;) */
} /* readsnap() */

write_snap()
{
        Body *b1,*b2;
        int  nbody_out;
        permanent bool first=TRUE;

        for (b1=btab, b2=btab; b1<btab+nbody; b1++) {
            if (Qbind && Mass(b1)<0)
                continue;               /* no copy */
            else if (!Qbind && Mass(b1)>0)
                continue;                /* no copy */

            if (Mass(b1) < 0)           /* negative mass was used as a signal */
                Mass(b1) = -Mass(b1);
            if (b1==b2) {
                b2++;
                continue;       /* no need to copy yet, still in sync */
            }
            bcopy(b1, b2, sizeof(Body));
            b2++;
        }
        nbody_out = b2-btab;
        dprintf (1,"Writing %d particles to output\n",nbody_out);
        if (first) {
            put_history(outstr);
            first = FALSE;
        }
#if 0
	if (nbody_out > 0)
            put_snap (outstr,&btab,&nbody_out,&stime,&bits);
        else
            warning("No stars to output");
#else
        put_snap (outstr,&btab,&nbody_out,&stime,&bits);
#endif
}

map()
{
    Body *bp;
    int i;
  
    for (bp=btab, i=0; i<nbody; bp++, i++) {
        if (i%50 == 0) printf("\n");
        if (Phi(bp) >= ecutoff)
            printf("*");    /*       a '*' for an unbound star */
        else
            printf(".");    /*      a  '.' for a bound star */
    }
}

/*
 * newton_potential exact 
 */
 
exact()
{
    Body *bi, *bj;
    int i,j,k;
    real rij;
    
    for (bi=btab; bi<btab+nbody; bi++)
        Phi(bi) = 0.0;
    for (i=1, bi=btab+1; i<nbody; i++, bi++) {
        for (j=0, bj=btab; j<i; bj++, j++) {
            rij = 0.0;
            for (k=0; k<NDIM; k++)
                rij += sqr(Pos(bi)[k] - Pos(bj)[k]);
            rij = 1.0/sqrt(sqreps + rij);
            Phi(bi) -= Mass(bj)*rij;            /* G==1 */
            Phi(bj) -= Mass(bi)*rij;            /* G==1 */
        }
    }
}
