/*
 *  convert Martel 'pos-vel' dumps to snapshot - in single precision
 */
 
#define SinglePrec

#include <nemo.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {
    "in=???\n       Input rvc file",
    "out=???\n      Output snapshot file",
    "nbody=262144\n  Max number of particles",
    "headline=\n    Random verbiage",
    "VERSION=1.0\n  19-mar-99 pjt",
    NULL,
};

string usage = "convert Martel \"pos-vel\" dump files to snapshot";

 
nemo_main()
{
    stream instr, outstr;
    real *fbuf, tsnap;
    long nread;
    int nbody, nmax;

    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");
    nmax = getiparam("nbody");
    if (hasvalue("headline"))
        set_headline(getparam("headline"));
    put_history(outstr);

    fbuf = (real *) allocate(nmax*6*sizeof(real));

    for (;;) {    
        nread = unfread(instr, fbuf, nmax*6*sizeof(real));
        if (nread<1) break;
        nbody = nread/6/sizeof(real);
        if (nread % (6*sizeof(real)))
            error("nbody=%d does not divide into %s",nbody,6*sizeof(real));

        tsnap = 0.0;
        put_set(outstr,SnapShotTag);
        put_set(outstr,ParametersTag);
        put_data(outstr, NobjTag, IntType, &nbody, 0);
        put_data(outstr, TimeTag, RealType, &tsnap, 0);
        put_tes(outstr,ParametersTag);
        put_set(outstr,ParticlesTag);
        put_data(outstr,PhaseSpaceTag,RealType,fbuf,nbody,2,NDIM,0);
        put_tes(outstr,ParticlesTag);
        put_tes(outstr,SnapShotTag);
    }
    free(fbuf);
    strclose(instr);
    strclose(outstr);
}




