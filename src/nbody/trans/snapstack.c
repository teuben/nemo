/*
 * SNAPSTACK: stack two N-body systems on top of each other.
 *
 *	22-jan-89  1.1 JEB
 *	 1-jul-90  1.1a added helpvec PJT
 *	 4-feb-93  1.1b nemo_main (also to return 0 exit)  -- PJT
 *			malloc -> allocate
 *	mar94 - ansi
 *	 6-aug-96  1.1d printf -> dprintf
 *	30-dec-97  1.1e ansi 
 */
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <snapshot/snapshot.h>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in1=???\n			  input file name ",
    "in2=???\n			  input file name ",
    "out=???\n			  output file name ",
    "deltar=0.0,0.0,0.0\n	  position of in1 wrt in2 ",
    "deltav=0.0,0.0,0.0\n	  velocity of in1 wrt in2 ",
    "zerocm=true\n		  zero center of mass ",
    "headline=\n		  random verbiage ",
    "VERSION=1.1e\n		  30-dec-97 PJT",
    NULL,
};

string usage="stack two N-body systems on top of each other";

extern string *burststring(string,string);

nemo_main()
{
    readdata();
    snapstack();
    writedata();
}

int nbody, nbody1, nbody2;

real *mass, *mass1, *mass2;

real *phase, *phase1, *phase2;

readdata()
{
    stream instr1, instr2;

    instr1 = stropen(getparam("in1"), "r");
    get_history(instr1);
    instr2 = stropen(getparam("in2"), "r");
    get_history(instr2);
    get_set(instr1, SnapShotTag);
    get_set(instr2, SnapShotTag);
    get_set(instr1, ParametersTag);
    get_set(instr2, ParametersTag);
    get_data(instr1, NobjTag, IntType, &nbody1, 0);
    get_data(instr2, NobjTag, IntType, &nbody2, 0);
    get_tes(instr1, ParametersTag);
    get_tes(instr2, ParametersTag);
    dprintf(1,"nbody1 = %d    nbody2 = %d\n", nbody1, nbody2);
    nbody = nbody1 + nbody2;
    mass = (real *) allocate(sizeof(real) * nbody);
    phase = (real *) allocate(sizeof(real) * 2*NDIM * nbody);
    mass1 = mass;
    mass2 = mass + nbody1;
    phase1 = phase;
    phase2 = phase + 2 * NDIM * nbody1;
    get_set(instr1, ParticlesTag);
    get_set(instr2, ParticlesTag);
    get_data_coerced(instr1, MassTag, RealType, mass1, nbody1, 0);
    get_data_coerced(instr2, MassTag, RealType, mass2, nbody2, 0);
    get_data_coerced(instr1, PhaseSpaceTag, RealType, phase1,
		     nbody1, 2, NDIM, 0);
    get_data_coerced(instr2, PhaseSpaceTag, RealType, phase2,
		     nbody2, 2, NDIM, 0);
    get_tes(instr1, ParticlesTag);
    get_tes(instr2, ParticlesTag);
    get_tes(instr1, SnapShotTag);
    get_tes(instr2, SnapShotTag);
}

snapstack()
{
    vector deltar, deltav;
    real *pp;

    setvect(deltar, getparam("deltar"));
    setvect(deltav, getparam("deltav"));
    for (pp = phase1; pp < phase2; ) {
	ADDV(pp, pp, deltar);
	pp += NDIM;
	ADDV(pp, pp, deltav);
	pp += NDIM;
    }
    if (getbparam("zerocm"))
	zerocms(phase, 2*NDIM, mass, nbody, nbody);
}

setvect(vec, str)
vector vec;
string str;
{
    string *vcp;
    int i;

    vcp = burststring(str, ", ");
    for (i = 0; i < NDIM; i++)
	vec[i] = (*vcp != NULL ? atof(*vcp++) : 0.0);
}

writedata()
{
    stream outstr;
    int cscode = CSCode(Cartesian, NDIM, 2);

    if (! streq(getparam("headline"), ""))
	set_headline(getparam("headline"));
    outstr = stropen(getparam("out"), "w");
    put_history(outstr);
    put_set(outstr, SnapShotTag);
    put_set(outstr, ParametersTag);
    put_data(outstr, NobjTag, IntType, &nbody, 0);
    put_tes(outstr, ParametersTag);
    put_set(outstr, ParticlesTag);
    put_data(outstr, CoordSystemTag, IntType, &cscode, 0);
    put_data(outstr, MassTag, RealType, mass, nbody, 0);
    put_data(outstr, PhaseSpaceTag, RealType, phase, nbody, 2, NDIM, 0);
    put_tes(outstr, ParticlesTag);
    put_tes(outstr, SnapShotTag);
}
