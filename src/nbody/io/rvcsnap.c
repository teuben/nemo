/*
 *  convert rvc to snapshot
 *  8-may-94  V1
 * 19-mar-99  V1.1  fixed set_headline args
 *  6-apr-99  V1.1a merged two versions
 */
 
#define SinglePrec

#include <nemo.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {
    "in=???\n       Input rvc file",
    "out=???\n      Output snapshot file",
    "headline=\n    Random verbiage",
    "VERSION=1.1a\n 6-apr-99 pjt",
    NULL,
};

string usage = "convert rvc files to snapshot";

 
nemo_main()
{
    stream instr, outstr;
    char buf[400];
    char *cp = buf;
    int *ip = (int *) cp;
    float *fp = (float *) cp;
    long nread,indx,foffset;
    int i, nbody;
    real *phase, *pp, tsnap;
    float *rv, *rvp;

    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");
    if (hasvalue("headline"))
        set_headline(getparam("headline"));
    put_history(outstr);

    for (;;) {    

        /* header 1 */ 
        nread = unfread(instr, buf, 400);
        if (nread<1) break;
        dprintf(0,"nobj %d ips %d aexp %f p %f dpnp1 %f dpn %f\n",
            ip[0],  ip[1], fp[2],  fp[3], fp[4], fp[5]);
        nbody = ip[0];
        tsnap = fp[3];
        phase = (real *) allocate(NDIM*2*nbody*sizeof(real));
        rv = (float *) allocate(NDIM*nbody*sizeof(float));

        /* header 2 */
        nread = unfread(instr, buf, 400);
        if (nread<1) error("cannot read 2nd header");

        /* pos */
        nread = unfread(instr, rv, nbody*NDIM*sizeof(float));
        if (nread<1) error("cannot read 3rd header");
        rvp = rv;
        pp = phase;
        for (i=0; i<nbody; i++) {
            *pp++ = *rvp++;     /* X */
            *pp++ = *rvp++;     /* Y */
            *pp++ = *rvp++;     /* Z */
            *pp += NDIM;        /* skip vel */
        }

        /* vel */
        nread = unfread(instr, rv, nbody*NDIM*sizeof(float));
        if (nread<1) error("cannot read 4th header");
        rvp = rv;
        pp = phase;
        for (i=0; i<nbody; i++) {
            *pp += NDIM;        /* skip pos */
            *pp++ = *rvp++;     /* VX */
            *pp++ = *rvp++;     /* VY */
            *pp++ = *rvp++;     /* VZ */
        }
        put_set(outstr,SnapShotTag);
        put_set(outstr,ParametersTag);
        put_data(outstr, NobjTag, IntType, &nbody, 0);
        put_data(outstr, TimeTag, RealType, &tsnap, 0);
        put_tes(outstr,ParametersTag);
        put_set(outstr,ParticlesTag);
        put_data(outstr,PhaseSpaceTag,RealType,phase,nbody,2,NDIM,0);
        put_tes(outstr,ParticlesTag);
        put_tes(outstr,SnapShotTag);
        
        
        free(phase);
        free(rv);
    }
}




