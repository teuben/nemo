/*
 * simple_write:  this code can be linked with a "uNEMO" version
 * 		  and run ok. So, no need for NEMO's user interface.
 *		  It need fake_getparam.c to give dummy code for
 *		  getparam() and finiparam()
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
    stream outstr = stropen("junk.dat","w");
    float mass[N],
          phase[N*2*NDIM], 
          phi[N],
          etot[3],
          keten[NDIM][NDIM], 
          peten[NDIM][NDIM], 
          amten[NDIM][NDIM], 
          cmphase[2][NDIM];
    
    put_set(outstr,SnapShotTag);
      put_set(outstr,ParametersTag);
        put_data(outstr, NobjTag, IntType, &nbody, 0);
        put_data(outstr, TimeTag, FloatType, &tsnap, 0);
      put_tes(outstr, ParametersTag);

      put_set(outstr, ParticlesTag);
        put_data(outstr, CoordSystemTag, IntType, &cs, 0);
        put_data(outstr, MassTag, FloatType, mass, nbody, 0);
        put_data(outstr, PhaseSpaceTag, FloatType, phase, nbody, 2, NDIM, 0);
        put_data(outstr, PotentialTag, FloatType, phi, nbody, 0);
      put_tes(outstr, ParticlesTag);

      put_set(outstr, DiagnosticsTag);
        put_data(outstr, EnergyTag, FloatType, etot, 3, 0);
        put_data(outstr, KETensorTag, FloatType, keten, NDIM, NDIM, 0);
        put_data(outstr, PETensorTag, FloatType, peten, NDIM, NDIM, 0);
        put_data(outstr, AMTensorTag, FloatType, amten, NDIM, NDIM, 0);
        put_data(outstr, CMPhaseSpaceTag, FloatType, cmphase, 2, NDIM, 0);
        put_data(outstr, "cputime", FloatType, &cput, 0);
      put_tes(outstr, DiagnosticsTag);


    put_tes(outstr,SnapShotTag);

    strclose(outstr);
}
