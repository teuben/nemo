/*
 *  Convert NBODYx output datasets (unit4 files) to snapshots
 *
 *  7-apr-93 	Created for NBODY1		pjt
 * 29-mar-94	ansi				pjt
 * 16-jun-97    fixed _ portability in name	pjt
 *  2-mar-06    header=                         PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <history.h>
#include <snapshot/snapshot.h>
#include "nbody_io.h"

string defv[] = {
    "in=\n          Input file (in UNIT 4 format)",
    "out=???\n      Output file (snapshot(5NEMO) format)",
    "nbody=\n       Input Number of particles",
    "header=\n      if used, force unfio header size (4 or 8)",
    "VERSION=1.2\n  4-mar-06 PJT",
    NULL,
};

string usage = "Convert NBODY unit-4 file to snapshot";

string cvsid="$Id$";

#define SIZEPP  36      /* 36 bytes per particle on fort.4 = 28 data + 8 header */

extern int nemo_file_size(string);

void nemo_main(void)
{
    string fname;
    int nbody, i, j, k, ibody;
    int coordsys = CSCode(Cartesian, 3, 2);
    float *mass, *vel, *pos;
    real *rmass, *rphase;
    stream outstr;

    if (hasvalue("header"))
      unfsize(getiparam("header"));

    if (hasvalue("in")) {
        fname = getparam("in");
        nb4open_(fname,strlen(fname));
    } else
        fname = "fort.4";

    if (hasvalue("nbody")) {
	nbody = getiparam("nbody");
    } else {
        nbody = nemo_file_size(fname);
        if (nbody%SIZEPP != 0) 
            error("Strange filesize %s; try setting nbody=",fname);
        nbody /= SIZEPP;
    }

    outstr = stropen(getparam("out"),"w");
    put_history(outstr);

    mass  = (float *) allocate(nbody*sizeof(float));
    pos   = (float *) allocate(nbody*sizeof(float)*3);
    vel   = (float *) allocate(nbody*sizeof(float)*3);

    nb4get_(&nbody,mass,pos,vel);
    if (nbody==0) error("Some error on unit 4");

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
    put_tes(outstr, ParametersTag);
    put_set(outstr, ParticlesTag);
    put_data(outstr, CoordSystemTag, IntType, &coordsys, 0);
    put_data(outstr, MassTag, RealType, rmass, nbody, 0);
    put_data(outstr,PhaseSpaceTag,RealType,rphase,nbody,2,3,0);
    put_tes(outstr, ParticlesTag);
    put_tes(outstr, SnapShotTag);

    free(rmass);
    free(rphase);

    free(mass);
    free(pos);
    free(vel);

    strclose(outstr);
    nb4close_();

}
    
