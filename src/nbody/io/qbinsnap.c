/*
 *  convert generic binary dumps to snapshot
 *     26-jun-2012   V0.1 - Q&D for QMOND (jun 7, 2012)
 *
 */
 
#include <nemo.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {
    "in=???\n                        Input rvc file",
    "out=???\n                       Output snapshot file",
    "nbody=0\n                       Number of bodies",
    "body=pos,vel\n                  Body variables",
    "np=256\n                        Np",
    "lbox=512\n                      Box   (MPc/H)",
    "h=72\n                          Hubble",
    "headline=\n                     Random verbiage",
    "VERSION=0.3\n                   26-jun-2012 PJT",
    NULL,
};

string usage = "convert binary files to snapshot";

 
nemo_main()
{
    stream instr, outstr;
    char buf[400];
    char *cp = buf;
    short  *sp = (short *)  cp;
    int    *ip = (int *)    cp;
    float  *fp = (float *)  cp;
    double *dp = (double *) cp;
    long nread,indx,foffset;
    short l;
    int i, nbody, nbuf;
    real *p1, *m1, *pp, tsnap, *t1, *t2;
    real *phase, *mass;
    real pscale, vscale, size, pshift;
    float *rv, *rvp;
    int np;
    real lbox, hub, ng;
    

    warning("Program only parses the QMOND 2012 data");

    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");
    nbody = getiparam("nbody");
    np = getiparam("np");
    lbox = getrparam("lbox");
    hub = getrparam("h");
    if (nbody==0) nbody = np*np*np;
    ng = np + 1.0;
    pshift = -1.0;
    pscale = lbox/ng;
    vscale = pscale*hub;
    dprintf(0,"Pscale=%g  Vscale=%g\n",pscale,vscale);
    if (hasvalue("headline"))
        set_headline(getparam("headline"));
    put_history(outstr);
    unfsize(4);
    nbuf = 48;    /* 6 doubles per particle */

    if (nbody==0) error("nbody=0  not supported yet");

    nread = unfread(instr, buf, 16);
    dp = (double *) buf;
    dprintf(0,"Header: %g   %g\n",dp[0],dp[1]);
    tsnap = 1.0/dp[0] - 1.0;

    phase = p1 = (real *) allocate(NDIM*2*nbody*sizeof(real));
    pp = (real *) allocate(128*sizeof(real));


    for (i=0; i<nbody; i++) {
      nread = unfread(instr, buf, nbuf);
      if (nread<1) error("cannot read particle %d",i);
      pp = (double *) buf;
      dprintf(1,"%d: pos %g %g %g\n",i+1,pp[0],pp[1],pp[2]);
      dprintf(1,"%d: vel %g %g %g\n",i+1,pp[3],pp[4],pp[5]);
      *p1++ = (*pp++ + pshift) * pscale;  /* X */
      *p1++ = (*pp++ + pshift) * pscale;
      *p1++ = (*pp++ + pshift) * pscale;
      *p1++ = *pp++ * vscale;             /* VX */
      *p1++ = *pp++ * vscale;
      *p1++ = *pp++ * vscale;
    }
    printf("done\n");
    put_set(outstr,SnapShotTag);
    put_set(outstr,ParametersTag);
    put_data(outstr, NobjTag, IntType, &nbody, 0);
    put_data(outstr, TimeTag, RealType, &tsnap, 0);
    put_tes(outstr,ParametersTag);
    put_set(outstr,ParticlesTag);
    put_data(outstr,PhaseSpaceTag,RealType,phase,nbody,2,NDIM,0);
    put_tes(outstr,ParticlesTag);
    put_tes(outstr,SnapShotTag);
}




