/*
 *  Convert NBODYx output datasets ('OUT3' files) to snapshots
 *
 *   7-apr-93   1.0  Created
 *   9-apr-93   1.1  NB3 to include NK length of A() array      PJT
 *  29-mar-94   1.2  optional alen= override; also ansi         pjt
 *  23-may-95   1.2a changed routine names in nbody_io.f        pjt
 *  10-jun-95   1.2c fixed default of a blank alen=		pjt
 *   8-aug-95   1.3  added a swap= keyword if to swap bytes     pjt
 *  16-jun-97   1.3a warning swap/FIO
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <history.h>
#include <snapshot/snapshot.h>
#include "nbody_io.h"

string defv[] = {
    "in=???\n       Input file (in NBODY's OUT3 format)",
    "out=???\n      Output file (snapshot(5NEMO) format)",
    "frame=0\n      Frames to read (0=all)",
    "alen=\n        If given, override length of A array (0=do not read A)",
    "swap=f\n       Swap byte (big vs. little endian only)",
    "integer=4\n    Size of integers in dataset (2 or 4)",
    "VERSION=1.3a\n 16-jun-97 PJT",
    NULL,
};

string usage = "Convert NBODY output to snapshot";


#define MAXFRAME  512
#define MAXHEADER  32

extern int file_size(string);


void nemo_main(void)
{
    string fname = getparam("in");
    int ilen, alen, saved=0, alen_fix=-1;
    int frame[MAXFRAME], nframe, iframe;
    int nbody, model, run, *name, i, j, k, ibody;
    int coordsys = CSCode(Cartesian, 3, 2);
    float *mass, *pos, *vel, a[MAXHEADER];
    real *rmass, *rphase, tsnap;
    bool Qswap;
    stream outstr;

    nframe = nemoinpi(getparam("frame"),frame,MAXFRAME);
    if(nframe<0) error("Parsing frame=%s",getparam("frame"));
    if(nframe==1 && frame[0]==0) nframe=0;
    if (hasvalue("alen")) alen_fix = getiparam("alen");
    Qswap = getbparam("swap");
    ilen = getiparam("integer");

    if (file_size(fname) < 0) error("File %s does not exist",fname);
    outstr = stropen(getparam("out"),"w");
    put_history(outstr);
#ifdef FIO
    nb3open_(fname,strlen(fname));
    if (Qswap) warning("-DFIO compiled cannot swap");
#else
    nb3open_c(fname, ilen, Qswap);
#endif

    for(iframe=0;;) {                         /* Loop over all frames */
      alen = alen_fix;
#ifdef FIO
      nb3header_(&nbody, &model, &run, &alen);
#else
      nb3header_c(&nbody, &model, &run, &alen);
#endif
      if (nbody == 0) break;     /* end of file */
      if (nbody<0 || nbody > 1000000) 
        warning("Strange value for nbody=%d; perhaps need to change to swap=%c",
                nbody, Qswap ? "false" : "true");
      if (alen_fix > 0) alen = alen_fix;		/* override alen */
      dprintf(1,"Header: nbody=%d model=%d run=%d nk=%d Snapshotsize=%d bytes\n",
                nbody,model,run,alen,
                32 + sizeof(float)*(alen+7*nbody) + sizeof(int)*nbody);
      if (nbody < 0) break;      /* something bad surely */
      if (alen<1 || alen>MAXHEADER)
         error("Bad headerlength nk=%d\n",alen);

      mass  = (float *) allocate(nbody*sizeof(float));
      pos   = (float *) allocate(nbody*sizeof(float)*3);
      vel   = (float *) allocate(nbody*sizeof(float)*3);
      name  = (int *)   allocate(nbody*sizeof(int));
#ifdef FIO
      nb3data_(&nbody,&alen,a,mass,pos,vel,name);
#else
      nb3data_c(&nbody,&alen,a,mass,pos,vel,name);
#endif

      dprintf(1,"Data  : a(%d)=",alen);
      for (i=0; i<alen; i++)
        dprintf(1," %g",a[i]);
      dprintf(1,"\n");

      if (nframe==0 || model==frame[iframe]) {
        tsnap = a[0];
        dprintf(0,"Saving model %d (nbody=%d time=%g)\n",model,nbody,tsnap);
        rmass  = (real *) allocate(nbody*sizeof(real));
        rphase = (real *) allocate(nbody*sizeof(real)*6);
        for (ibody=0, i=0, j=0, k=0; ibody<nbody ; ibody++) {
            rmass[ibody] = mass[ibody];
            rphase[k++] = pos[i++];
            rphase[k++] = pos[i++];
            rphase[k++] = pos[i++];
            rphase[k++] = vel[j++];
            rphase[k++] = vel[j++];
            rphase[k++] = vel[j++];
        }
        put_set(outstr, SnapShotTag);
        put_set(outstr, ParametersTag);
        put_data(outstr, NobjTag, IntType, &nbody, 0);
        put_data(outstr, TimeTag, RealType, &tsnap, 0);
        put_tes(outstr, ParametersTag);
        put_set(outstr, ParticlesTag);
        put_data(outstr, CoordSystemTag, IntType, &coordsys, 0);
        put_data(outstr, MassTag, RealType, rmass, nbody, 0);
        put_data(outstr,PhaseSpaceTag,RealType,rphase,nbody,2,3,0);
        put_data(outstr,KeyTag,IntType,name,nbody,0);
        put_tes(outstr, ParticlesTag);
        put_tes(outstr, SnapShotTag);
        free(rmass);
        free(rphase);
        iframe++;
        saved++;
      }

      free(mass);
      free(pos);
      free(vel);
      free(name);

      if(nframe>0 && iframe>=nframe) break;      /* early stop */
    }

#ifdef FIO
    nb3close_();
#else
    nb3close_c();
#endif
    strclose(outstr);

    dprintf(0,"Saved %d snapshots\n",saved);

}
    
