/*
 * SNAPSHOT.H: definitions for N-body snapshot files.
 *	sep-1990	added Story, Pos and Vel options
 */

#ifndef _snapshot_h
#define _snapshot_h

/*
 * Item tags for SnapShot components.
 */

#define SnapShotTag		"SnapShot"

#define   ParametersTag		"Parameters"
#define     NobjTag             "Nobj"
#define     TimeTag             "Time"

#define   ParticlesTag		"Particles"
#define     CoordSystemTag      "CoordSystem"
#define     MassTag             "Mass"
#define     PhaseSpaceTag       "PhaseSpace"
#define     PosTag              "Pos"
#define     VelTag              "Vel"
#define	    AccelerationTag	"Acceleration"
#define	    PotentialTag	"Potential"
#define     AuxTag              "Aux"
#define     KeyTag              "Key"
#define     StoryTag		"Story"

#define   DiagnosticsTag	"Diagnostics"
#define     EnergyTag           "Energy"
#define     KETensorTag         "KETensor"
#define     PETensorTag		"PETensor"
#define     AMTensorTag		"AMTensor"
#define     CMPhaseSpaceTag	"CMPhaseSpace"

/*
 * Symbolic names for input/output bit flags.
 *	Assumed 32 bits are available
 */

#define TimeBit          (1 <<  0)

#define MassBit          (1 <<  1)
#define PhaseSpaceBit    (1 <<  2)
#define PotentialBit     (1 <<  3)
#define AccelerationBit  (1 <<  4)
#define AuxBit           (1 <<  5)
#define KeyBit           (1 <<  6)

#define EnergyBit        (1 <<  7)
#define KETensorBit      (1 <<  8)
#define PETensorBit      (1 <<  9)
#define AMTensorBit      (1 << 10)
#define CMPhaseSpaceBit  (1 << 11)

#define StoryBit         (1 << 12)
#define PosBit		 (1 << 13)
#define VelBit		 (1 << 14)

/*
 * Coordinate system codes; these assume 32-bit ints.
 */

#define CSCode(type,ndim,nder) ((type) | (ndim) << 8 | (nder))

#define CSType(code) ( (code) & (0377 << 16))
#define CSNdim(code) (((code) >> 8) & 0377)
#define CSNder(code) ( (code) & 0377)

#define Cartesian   (1 << 16)
#define Spherical   (2 << 16)
#define Scattering  (3 << 16)

#define TrueAnomaly (4 << 16)
#define EccAnomaly  (5 << 16)
#define MeanAnomaly (6 << 16)
#define PeriPassage (7 << 16)

#endif
