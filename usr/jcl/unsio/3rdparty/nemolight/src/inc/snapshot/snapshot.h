/*
 * SNAPSHOT.H: definitions for N-body snapshot files.
 *	sep-1990	added Story, Pos and Vel options
 *      may-2001        added Yanc stuff for Dehnen
 *      oct-2001        added SPH/Yanc (dens/eps)
 *      jan-2002        added ZENO compatible ones
 *      nov-2003        removed Yanc tags (YANC is also called gyrfalcON now)
 *      feb-2004        added some new SPH stuff (GasDensity, NPartners, NSPHPartners)
 *      may-2010        added in a few tag from atos.c for handling its' SPH 
 */

#ifndef _snapshot_h
#define _snapshot_h

/*
 * Item tags for SnapShot components. 
 * !! see also some ZENO compatibility components below !!
 * Handling of SPH particles is still done rather poorly 
 * the current NEMO convention is to store the SPH particles first
 */

#define SnapShotTag		"SnapShot"

#define   ParametersTag		"Parameters"
#define     NobjTag             "Nobj"
#define     NgasTag             "Ngas"        /* atos.c */
#define     TimeTag             "Time"

#define   ParticlesTag		"Particles"
#define     CoordSystemTag      "CoordSystem"
#define     MassTag             "Mass"
#define     PhaseSpaceTag       "PhaseSpace"
#define     PosTag              "Position"
#define     VelTag              "Velocity"
#define	    AccelerationTag	"Acceleration"
#define	    PotentialTag	"Potential"
#define     AuxTag              "Aux"
#define     KeyTag              "Key"
#define     DensityTag          "Density"
#define     GasDensTag          "GasDensity"
#define     TemperatureTag      "Temperature"   /* atos.c */
#define     SmoothLengthTag     "SmoothLength"  /* atos.c */
#define     EpsTag              "Eps"
#define     StoryTag		"Story"

#define     NumberTag           "NPartners"
#define     SPHNumberTag        "NSPHPartners"


#define   DiagnosticsTag	"Diagnostics"
#define     EnergyTag           "Energy"
#define     KETensorTag         "KETensor"
#define     PETensorTag		"PETensor"
#define     AMTensorTag		"AMTensor"
#define     CMPhaseSpaceTag	"CMPhaseSpace"
#define     CPUTimeTag          "cputime"

/* Some ZENO compatible tags */

#define	AccTag         	AccelerationTag
#define CMPhaseTag     	CMPhaseSpaceTag
#define RhoTag         	DensityTag
#define PhaseTag       	PhaseSpaceTag


/* 
 *  New ZENO tags that NEMO had not defined before 
 *  Note the NGas  vs. Ngas
 *           NBody vs. Nobj
 */

#define NBodyTag 	"NBody"
#define NGasTag 	"NGas"
#define UdotIntTag 	"UdotInternal"
#define UdotRadTag 	"UdotRadiation"
#define UinternTag 	"Uinternal"
#define EntFuncTag 	"EntropyFunc"
#define AMVectorTag 	"AMVector"
#define AuxVecTag 	"AuxVec"
#define BodyTag 	"Body"
#define SmoothTag       "SmoothLength"


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

#define DensBit          (1 << 15)
#define EpsBit           (1 << 16)

/* Note:  there should be a more clear separation of grav.softening (particularly if
 *	  variable) and what we probably mean as SPH smoothing length "Eps" in this
 *  	  context
 */

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
