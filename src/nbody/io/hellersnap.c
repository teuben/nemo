/*
 *  convert Heller SPH data to standard (single precision) snapshot
 *
 *  25-jun-98   PJT     Created         
 */
 
#define SINGLEPREC

#include <nemo.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {
    "in=???\n       Input rvc file",
    "out=???\n      Output snapshot file",
    "headline=\n    Random verbiage",
    "VERSION=1.0\n  26-jun-98 pjt",
    NULL,
};

string usage = "convert Heller's SPH dumps to snapshot";


local void scatter(int n, real *from, real *to, int skip);
 
nemo_main()
{
    stream instr, outstr;
    char buf[400];
    char *cp = buf;
    int *ip = (int *) cp;
    float *fp = (float *) cp;
    double *dp = (double *) buf;
    long nread, indx, foffset;
    int i, nbody, ndim, dmpindx, n1, n2, *key;
    real *phase, *pp, tsnap, *mass, *aux, *phi, *acc;
    float *rv, *rvp, *rbuf;

    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");
    if (hasvalue("headline"))
        set_headline(getparam("headline"));
    put_history(outstr);

    for (;;) {    

        /* header 1 */ 
        nread = unfread(instr, buf, 400);
        if (nread<1) break;

        dmpindx = ip[0];
        ndim = ip[7];
        if (ndim != NDIM)
            error("No support for ndim=%d (NDIM=%d)\n",ndim,NDIM);

        /* header 2 */
        nread = unfread(instr, buf, 400);
        if (nread<1) error("cannot read 2nd header");
        
        dprintf(0,"Time %g N1 %d N2 %d  DMPINDX=%d  NDIM=%d\n",
               dp[0],  ip[2], ip[3], dmpindx, ndim);

        n1 = ip[2];
        n2 = ip[3];            
        nbody = n1 + n2;
        tsnap = dp[0];
        
        rbuf = (real *) allocate(nbody*sizeof(real));

        key = (int *) allocate(nbody*sizeof(int));        

        nread = unfread(instr, key, nbody * sizeof(int));
        if (nread<1) error("cannot read CLASS");

        nread = unfread(instr, rbuf, nbody * sizeof(real));     /* ignore */
        if (nread<1) error("cannot read STATE");

        mass = (real *) allocate(nbody*sizeof(real));
        nread = unfread(instr, mass, nbody * sizeof(real));
        if (nread<1) error("cannot read MASS");

        phase = (real *) allocate(NDIM*2*nbody*sizeof(real));

        for (i=0; i<ndim; i++) {
            nread = unfread(instr, rbuf, nbody * sizeof(real));
            if (nread<1) error("cannot read POS(%d)",i+1);
            scatter(nbody,rbuf,&phase[i],2*ndim);
        }
        for (i=0; i<ndim; i++) {
            nread = unfread(instr, rbuf, nbody * sizeof(real));
            if (nread<1) error("cannot read VEL(%d)",i+1);
            scatter(nbody,rbuf,&phase[i+ndim],2*ndim);            
        }

        if (dmpindx > 1) {

            nread = unfread(instr, rbuf, nbody * sizeof(real));
            if (nread<1) error("cannot read EPS");

            phi = (real *) allocate(nbody*sizeof(real));        
            nread = unfread(instr, phi, nbody * sizeof(real));
            if (nread<1) error("cannot read PHI");

            acc = (real *) allocate(NDIM*nbody*sizeof(real));
            for (i=0; i<ndim; i++) {
                nread = unfread(instr, rbuf, nbody * sizeof(real));
                if (nread<1) error("cannot read ACC(%d)",i+1);
                scatter(nbody,rbuf,&acc[i],ndim);                
            }

        }

        if (dmpindx > 2 && n2 > 0) {

            warning("Discarding U, D, DIVV, UDOT1, UDOT2, H");
            nread = unfread(instr, rbuf, n2 * sizeof(real));
            if (nread<1) error("cannot read U");

            nread = unfread(instr, rbuf, n2 * sizeof(real));
            if (nread<1) error("cannot read D");

            nread = unfread(instr, rbuf, n2 * sizeof(real));
            if (nread<1) error("cannot read DIVV");

            nread = unfread(instr, rbuf, n2 * sizeof(real));
            if (nread<1) error("cannot read UDOT1");
            
            nread = unfread(instr, rbuf, n2 * sizeof(real));
            if (nread<1) error("cannot read UDOT2");

            nread = unfread(instr, rbuf, n2 * sizeof(real));
            if (nread<1) error("cannot read H");
        }


        put_set(outstr,SnapShotTag);
        put_set(outstr,ParametersTag);
        put_data(outstr, NobjTag, IntType, &nbody, 0);
        put_data(outstr, TimeTag, RealType, &tsnap, 0);
        put_tes(outstr,ParametersTag);
        put_set(outstr,ParticlesTag);
        put_data(outstr,MassTag,RealType,mass,nbody,0);
        put_data(outstr,PhaseSpaceTag,RealType,phase,nbody,2,NDIM,0);
        put_data(outstr,KeyTag,IntType,key,nbody,0);
        if (dmpindx > 1) {
            put_data(outstr,PotentialTag,RealType,phi,nbody,0);
            put_data(outstr,AccelerationTag,RealType,acc,nbody,NDIM,0);
        }
        put_tes(outstr,ParticlesTag);
        put_tes(outstr,SnapShotTag);
        
        
    }
}



local void scatter(int n, real *from, real *to, int skip)
{
    int i, j;

    for (i=0, j=0; i<n; i++, j+=skip)
        to[j] = from[i];
}
