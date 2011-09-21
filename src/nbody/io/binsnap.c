/*
 *  convert generic binary dumps to snapshot
 *     20-sep-2011   V0.1 - Q&D for KxM2007 models
 *
 * @todo
 *      types= needs to be implemented, 1,2,4,8 means skip that many bytes
 *      body=  needs to be parsed
 */
 
#include <nemo.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {
    "in=???\n                        Input rvc file",
    "out=???\n                       Output snapshot file",
    "nbody=0\n                       Number of bodies",
    "body=pos,vel,den,Tr,Tg,l\n      Body variables",
    "type=ddddddddds\n               Types (double,float,int,short,1,2,4,8)",
    "mscale=1e-21\n                  Scale densities by",
    "pscale=1e17\n                   Scale position by",
    "vscale=1e5\n                    Scale velocities by",
    "dr=1.448e16\n                   Cell size for level=0",
    "mass=t\n                        Use density (f) or mass (t)",
    "headline=\n                     Random verbiage",
    "VERSION=0.2\n                   21-sep-2011 PJT",
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
    real mscale, pscale, vscale, dr, size;
    float *rv, *rvp;
    bool Qmass;

    warning("Program only parses the KxM2007 data");

    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");
    nbody = getiparam("nbody");
    mscale = getrparam("mscale");
    pscale = getrparam("pscale");
    vscale = getrparam("vscale");
    dr = getrparam("dr");
    Qmass = getbparam("mass");
    if (hasvalue("headline"))
        set_headline(getparam("headline"));
    put_history(outstr);
    unfsize(0);
    nbuf = 74;    /* 9 doubles and a short */

    if (nbody==0) error("nbody=0  not supported yet");

    phase = p1 = (real *) allocate(NDIM*2*nbody*sizeof(real));
    mass =  m1 = (real *)  allocate(nbody*sizeof(real));
    t1 = (real *)  allocate(nbody*sizeof(real));
    t2 = (real *)  allocate(nbody*sizeof(real));
    pp = (real *) allocate(128*sizeof(real));



    for (i=0; i<nbody; i++) {
      nread = unfread(instr, buf, nbuf);
      if (nread<1) error("cannot read particle %d",i);
      pp = (double *) buf;
      dprintf(1,"%d: pos %g %g %g\n",i+1,pp[0],pp[1],pp[2]);
      dprintf(1,"%d: vel %g %g %g\n",i+1,pp[3],pp[4],pp[5]);
      dprintf(1,"%d: mtt %g %g %g\n",i+1,pp[6],pp[7],pp[8]);
      *p1++ = *pp++ / pscale;  /* X */
      *p1++ = *pp++ / pscale;
      *p1++ = *pp++ / pscale;
      *p1++ = *pp++ / vscale;  /* VX */
      *p1++ = *pp++ / vscale;
      *p1++ = *pp++ / vscale;
      *m1++ = *pp++ / mscale;  /* mass */
      *t1++ = *pp++;  /* t1 */
      *t2++ = *pp++;  /* t2 */
      sp    =  pp;
      l     = (short ) *sp;
      dprintf(1,"%d: l=%d (0x%x)\n",i+1,l,l);
      if (Qmass) {
	size  = pow(2.0,(double)l);
	mass[i] /= size;
      }
      dprintf(0,"mass %d %g\n",i,mass[i]);
    }
    printf("done\n");
    put_set(outstr,SnapShotTag);
    put_set(outstr,ParametersTag);
    put_data(outstr, NobjTag, IntType, &nbody, 0);
    put_data(outstr, TimeTag, RealType, &tsnap, 0);
    put_tes(outstr,ParametersTag);
    put_set(outstr,ParticlesTag);
    put_data(outstr,MassTag,RealType,mass,nbody,0);
    put_data(outstr,PhaseSpaceTag,RealType,phase,nbody,2,NDIM,0);
    put_tes(outstr,ParticlesTag);
    put_tes(outstr,SnapShotTag);
}




