/*
 * simple_read:   this code can be linked with a "uNEMO" version
 * 		  and run ok. So, no need for NEMO's user interface.
 *		  It need fake_getparam.c to give dummy code for
 *		  getparam() and finiparam()
 *
 *  Notice that this routine is the EXACT counterpart (apart from
 *  the optional get_history()) of simple_write(), with the following
 *  replacement:
 *	outstr  -> instr
 *	put_XXX -> put_XXX
 *
 *	21-may-2001	Written in Mexico, for Walter Dehnen		PJT
 */

#include <stdinc.h>
#include <filestruct.h>
#include <snapshot.h>

#define N     3
#define NDIM  3

main()
{
    int nbody = N;
    float tsnap = 1.0, cput = 0.0;
    int cs = CSCode(Cartesian, NDIM, 2);
    stream instr = stropen("junk.dat","r");
    float mass[N],
          phase[N*2*NDIM], 
          phi[N],
          etot[3],
          keten[NDIM][NDIM], 
          peten[NDIM][NDIM], 
          amten[NDIM][NDIM], 
          cmphase[2][NDIM];
    
    get_history(instr);			/* optional */
    get_set(instr,SnapShotTag);
      get_set(instr,ParametersTag);
        get_data(instr, NobjTag, IntType, &nbody, 0);
        get_data(instr, TimeTag, FloatType, &tsnap, 0);
      get_tes(instr, ParametersTag);
      printf("Found snapshot with Nbody=%d, Time=%g\n",nbody,tsnap);

      get_set(instr, ParticlesTag);
        get_data(instr, CoordSystemTag, IntType, &cs, 0);
        get_data(instr, MassTag, FloatType, mass, nbody, 0);
        get_data(instr, PhaseSpaceTag, FloatType, phase, nbody, 2, NDIM, 0);
        get_data(instr, PotentialTag, FloatType, phi, nbody, 0);
      get_tes(instr, ParticlesTag);

      get_set(instr, DiagnosticsTag);
        get_data(instr, EnergyTag, FloatType, etot, 3, 0);
        get_data(instr, KETensorTag, FloatType, keten, NDIM, NDIM, 0);
        get_data(instr, PETensorTag, FloatType, peten, NDIM, NDIM, 0);
        get_data(instr, AMTensorTag, FloatType, amten, NDIM, NDIM, 0);
        get_data(instr, CMPhaseSpaceTag, FloatType, cmphase, 2, NDIM, 0);
        get_data(instr, "cputime", FloatType, &cput, 0);
      get_tes(instr, DiagnosticsTag);


    get_tes(instr,SnapShotTag);

    strclose(instr);
}
