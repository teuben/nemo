/*
 * aid functions for stod, but in C to avoid too many collisions
 * between NEMO and STARLAB and vice versa.
 *
 * The C++ declarations live in stod_subs.h
 */


#include <stdinc.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <snapshot/body.h>
#include <snapshot/snapshot.h>

void check_real(int size)
{
    if (sizeof(real) != size) 
        error("NEMO is using %d-byte real, STARLAB is using %d real",
            sizeof(real), size);
}

int get_snap_c(string fname, real *time, real **mass, real **pos, real **vel)
{
    int i, nbody;
    real *mptr, *pptr, *vptr, *phase;
    stream instr;

    *time = -1.0;    /* for now */

    instr = stropen(fname,"r");                     /* open file */
    get_history(instr);                             /* add data history */
    if (!get_tag_ok(instr, SnapShotTag)) return 0;
    get_set(instr, SnapShotTag);                    /* and get first dataset */
    get_set(instr, ParametersTag);
    get_data(instr, NobjTag, IntType, &nbody, 0);

    get_tes(instr, ParametersTag);

    if (!get_tag_ok(instr, ParticlesTag)) return 0;
    get_set(instr, ParticlesTag);
    if (get_tag_ok(instr, MassTag)) {
        *mass = mptr  = (real *) allocate(nbody * sizeof(real));
        get_data_coerced(instr, MassTag, RealType, mptr,nbody, 0);
    } else
        *mass = NULL;

    if (get_tag_ok(instr, PhaseSpaceTag)) {
        phase = (real *) allocate(2 * NDIM * nbody * sizeof(real));
        *pos  = pptr  = (real *) allocate(NDIM * nbody * sizeof(real));
        *vel  = vptr  = (real *) allocate(NDIM * nbody * sizeof(real));
        get_data_coerced(instr, PhaseSpaceTag, RealType, phase,nbody,2,NDIM,0);
        for (i=0; i<nbody; i++) {
            *pptr++ = *phase++;
            *pptr++ = *phase++;
            *pptr++ = *phase++;
            *vptr++ = *phase++;       
            *vptr++ = *phase++;       
            *vptr++ = *phase++;
        }       
    }
    return nbody;
}


void put_snap_c(string fname, string headline, int nbody, real time, 
	real *mass, real *pos, real *vel, real *acc, real *aux, real *phi, int *key)
{
    int i;
    real *mptr, *pptr, *vptr, *ptr, *phase;
    static stream instr;
    static int n=0;
    int cs = CSCode(Cartesian, NDIM, 2);

    if (n == 0) {
	dprintf(1,"Opening %s to write\n",fname);
	instr = stropen(fname,"w");                     /* open file */
	if (strlen(headline))
	    set_headline(headline);
	put_history(instr);                             /* add data history */
    } else
	dprintf(1,"Appending more data\n");
    n++;

    put_set(instr, SnapShotTag);                    /* and put first dataset */

    put_set(instr, ParametersTag);
    put_data(instr, NobjTag, IntType, &nbody, 0);
    put_data(instr, TimeTag, RealType, &time, 0);
    put_tes(instr, ParametersTag);


    put_set(instr, ParticlesTag);
    put_data(instr, CoordSystemTag, IntType, &cs, 0);
    put_data(instr, MassTag, RealType, mass, nbody, 0);

    ptr = phase = (real *) allocate(2 * NDIM * nbody * sizeof(real));
    pptr = pos;
    vptr = vel;
    for (i=0; i<nbody; i++) {
            *ptr++ = *pptr++;
            *ptr++ = *pptr++;
            *ptr++ = *pptr++;
            *ptr++ = *vptr++;
            *ptr++ = *vptr++;
            *ptr++ = *vptr++;
    }       
    put_data(instr, PhaseSpaceTag, RealType, phase, nbody,2,NDIM,0);
    put_data(instr, AuxTag, RealType, aux, nbody, 0);
    put_data(instr, AccelerationTag, RealType, acc, nbody,NDIM, 0);
    put_data(instr, PotentialTag, RealType, phi, nbody, 0);
    put_data(instr, KeyTag, IntType, key, nbody, 0);
    put_tes(instr, ParticlesTag);
    put_tes(instr, SnapShotTag);

}


