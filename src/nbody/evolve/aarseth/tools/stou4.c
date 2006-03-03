/*
 * Convert snapshot to unit4 for NBODYx
 *	7-apr-93  created
 *     29-mar-94  ansi
 *     16-jun-97  _ naming portability
 *     11-jul-97  output of nbody 
 *      2-mar-06  header=                        
 */

#include <stdinc.h>
#include <assert.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <snapshot/snapshot.h>
#include "nbody_io.h"

string defv[] = {          
    "in=???\n                     input file name",
    "out=\n                       output file name [fort.4]",
    "nbody=0\n                    test if this nbody is ok (0=no test)",
    "header=4\n                   unfio header size (4 or 8)",
    "VERSION=1.2\n		  2-mar-06 PJT",
    NULL,
};

string usage = "convert snapshot to unit-4 for NBODYx";

extern int nemo_file_size(string);

void nemo_main(void)
{
    stream instr;
    string fname;
    int ibody, nbody, maxbody, i, j, k;
    real *mbuf, *cbuf;
    float *mass, *vel, *pos;
    int hdr_size = getiparam("header");

    maxbody = getiparam("nbody");
    instr = stropen(getparam("in"), "r");
    unfsize(hdr_size);
    get_history(instr);
    if (hasvalue("out")) {
        fname = getparam("out");
        if (streq(fname,".") || streq(fname,"-"))
            error("Cannot use special NEMO I/O names for out=");
        nb4open_(fname,strlen(fname));
    } else
        fname = "fort.4";
    if (nemo_file_size(fname) > 0) warning("Overwriting %s",fname);

	get_set(instr, SnapShotTag);
	get_set(instr, ParametersTag);
	get_data_coerced(instr, NobjTag, IntType, &nbody, 0);
	if (maxbody>0 && nbody != maxbody)
            warning("Incorrect nbody=%d, should be %d",nbody,maxbody);
	get_tes(instr, ParametersTag);

	get_set(instr, ParticlesTag);
	mbuf = (real *) allocate(nbody * sizeof(real));
	get_data_coerced(instr, MassTag, RealType, mbuf, nbody, 0);

	cbuf = (real *) allocate(nbody * 2 * NDIM * sizeof(real));
        get_data_coerced(instr, PhaseSpaceTag, RealType, cbuf,
			     nbody, 2, NDIM, 0);

        mass = (float *)allocate(nbody * sizeof(float));
        pos = (float *)allocate(nbody * sizeof(float)*3);
        vel = (float *)allocate(nbody * sizeof(float)*3);

        for (ibody=0, i=0, j=0, k=0; ibody<nbody; ibody++) {
            mass[ibody] = mbuf[ibody];
            pos[i++] = cbuf[k++];
            pos[i++] = cbuf[k++];
            pos[i++] = cbuf[k++];
            vel[j++] = cbuf[k++];            
            vel[j++] = cbuf[k++];            
            vel[j++] = cbuf[k++];
        }
        nb4put_(&nbody,mass,pos,vel);

	get_tes(instr, ParticlesTag);
	get_tes(instr, SnapShotTag);

    strclose(instr);
    nb4close_();
}

