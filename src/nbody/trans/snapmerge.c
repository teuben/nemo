/*
 * SNAPMERGE: merge N-body snapshots from multiple times (from 1 file) together
 *            into one snapshot - only masses and phase space coordinates
 *            survive.
 *
 *	28-jun-92  1.0 Created
 *	22-jul-92  1.0a    added history output		PJT
 *	mar-94 ansi
 *	19-jan-98  1.1  adding nbody= and nsnap=	pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <snapshot/snapshot.h>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n		input file names",
    "out=???\n		output file name",
    "nbody=\n           Size of one (or more) snapshot(s)",
    "nsnap=\n           Number of pieces to cut a snapshot into",
    "VERSION=1.1\n	19-jan-98 PJT",
    NULL,
};

string usage="merge N-body snapshots together";

stream instr, outstr;
int nbodytot=0, nbody;
real *masstot=NULL, *mass;
real *phasetot=NULL, *phase;
int  *keytot=NULL, *key;
bool Qmass;
real tsnap;


void nemo_main() 
{
    instr = stropen(getparam("in"), "r");
    outstr = stropen(getparam("out"), "w");

    readdata();
    writedata();
}

readdata()
{
    int i, offset=0, n=0;

    for(;;) {                   /* loop forever until no more data */
        get_history(instr);
        if (!get_tag_ok(instr,SnapShotTag)) break;

        get_set(instr, SnapShotTag);
        get_set(instr, ParametersTag);
        get_data(instr, NobjTag, IntType, &nbody, 0);
        if(get_tag_ok(instr,TimeTag))
             get_data(instr, TimeTag, RealType, &tsnap, 0);
        else
             tsnap=0.0;
        get_tes(instr, ParametersTag);

        if (!get_tag_ok(instr,ParticlesTag)) {
            get_tes(instr, SnapShotTag);
	    continue;
	}

        get_set(instr,ParticlesTag);

        if (nbodytot==0) {
            nbodytot = nbody;
            masstot = (real *) allocate(sizeof(real) * nbodytot); 
            for (i=0; i<nbody; i++)
                masstot[i] = 0.0;
            phasetot = (real *) allocate(sizeof(real) * 2*NDIM * nbodytot);
	    if (get_tag_ok(instr,KeyTag))
	      keytot = (int *) allocate(sizeof(int) * nbodytot); 
        } else {
            nbodytot += nbody;
            masstot = (real *) reallocate(masstot, sizeof(real) * nbodytot); 
            phasetot = (real *) reallocate(phasetot,
                                    sizeof(real) * 2*NDIM * nbodytot);
	    if (keytot)
	      keytot = (int *) reallocate(keytot, sizeof(int) * nbodytot); 
        }



        mass = masstot + offset;
        phase = phasetot + 2 * NDIM * offset;
	if (keytot)
	  key = keytot + offset;

        if (get_tag_ok(instr, MassTag)) {
            get_data_coerced(instr, MassTag, RealType, mass, nbody, 0);
            Qmass = TRUE;
        } else {
            for (i=0; i<nbody; i++) mass[i] = masstot[i];   /* no mass change */
        }
        get_data_coerced(instr, PhaseSpaceTag, RealType, phase,
		     nbody, 2, NDIM, 0);
	if (keytot)
	  get_data(instr, KeyTag, IntType, key, nbody, 0);
        get_tes(instr, ParticlesTag);
        get_tes(instr, SnapShotTag);
	offset = nbodytot;
    }
}


writedata()
{
    permanent int cscode = CSCode(Cartesian, NDIM, 2);

    put_history(outstr);

    put_set(outstr, SnapShotTag);
      put_set(outstr, ParametersTag);
        put_data(outstr, NobjTag, IntType, &nbodytot,0);
        put_data(outstr, TimeTag, RealType, &tsnap,0);
      put_tes(outstr, ParametersTag);
      put_set(outstr, ParticlesTag);
        put_data(outstr, CoordSystemTag, IntType, &cscode,0);
	if (Qmass)
          put_data(outstr, MassTag, RealType, masstot, nbodytot,0);
        put_data(outstr, PhaseSpaceTag, RealType, phasetot, nbodytot,2,NDIM,0);
	if (keytot)
	  put_data(outstr, KeyTag, IntType, keytot, nbodytot,0);
      put_tes(outstr, ParticlesTag);
    put_tes(outstr, SnapShotTag);
}
