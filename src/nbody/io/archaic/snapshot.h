/*
 * SNAPSHOT.H: structured binary file definitions for standard
 * N-body snapshot files.
 */

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
#define	    AccelerationTag	"Acceleration"
#define	    PotentialTag	"Potential"
#define     AuxTag              "Aux"
#define     KeyTag              "Key"

#define   DiagnosticsTag	"Diagnostics"
#define     EnergyTag           "Energy"
#define     KETensorTag         "KETensor"
#define     PETensorTag		"PETensor"
#define     AMTensorTag		"AMTensor"
#define     CMPhaseSpaceTag	"CMPhaseSpace"

/*
 * Accessor functions for SnapShot components.
 */

#define ParametersItem(s)   scanitem(s, ParametersTag, NULL)
#define NobjItem(s)         scanitem(s, ParametersTag, NobjTag, NULL)
#define TimeItem(s)         scanitem(s, ParametersTag, TimeTag, NULL)

#define ParticlesItem(s)    scanitem(s, ParticlesTag, NULL)
#define CoordSystemItem(s)  scanitem(s, ParticlesTag, CoordSystemTag, NULL)
#define MassItem(s)         scanitem(s, ParticlesTag, MassTag, NULL)
#define PhaseSpaceItem(s)   scanitem(s, ParticlesTag, PhaseSpaceTag, NULL)
#define AccelerationItem(s) scanitem(s, ParticlesTag, AccelerationTag, NULL)
#define PotentialItem(s)    scanitem(s, ParticlesTag, PotentialTag, NULL)
#define AuxItem(s)          scanitem(s, ParticlesTag, AuxTag, NULL)
#define KeyItem(s)          scanitem(s, ParticlesTag, KeyTag, NULL)

#define DiagnosticsItem(s)  scanitem(s, DiagnosticsTag, NULL)
#define EnergyItem(s)       scanitem(s, DiagnosticsTag, EnergyTag, NULL)
#define KETensorItem(s)     scanitem(s, DiagnosticsTag, KETensorTag, NULL)
#define PETensorItem(s)     scanitem(s, DiagnosticsTag, PETensorTag, NULL)
#define AMTensorItem(s)     scanitem(s, DiagnosticsTag, AMTensorTag, NULL)
#define CMPhaseSpaceItem(s) scanitem(s, DiagnosticsTag, CMPhaseSpaceTag, NULL)

#define NobjData(s)         ((int  *) ItemDat(NobjItem(s)))
#define TimeData(s)         ((real *) ItemDat(TimeItem(s)))

#define CoordSystemData(s)  ((int  *) ItemDat(CoordSystemItem(s)))
#define MassData(s)         ((real *) ItemDat(MassItem(s)))
#define PhaseSpaceData(s)   ((real *) ItemDat(PhaseSpaceItem(s)))
#define AccelerationData(s) ((real *) ItemDat(AccelerationItem(s)))
#define PotentialData(s)    ((real *) ItemDat(PotentialItem(s)))
#define AuxData(s)          ((real *) ItemDat(AuxItem(s)))
#define KeyData(s)          ((int  *) ItemDat(KeyItem(s)))

#define EnergyData(s)       ((real *) ItemDat(EnergyItem(s)))
#define KETensorData(s)     ((real *) ItemDat(KETensorItem(s)))
#define PETensorData(s)     ((real *) ItemDat(PETensorItem(s)))
#define AMTensorData(s)     ((real *) ItemDat(AMTensorItem(s)))
#define CMPhaseSpaceData(s) ((real *) ItemDat(CMPhaseSpaceItem(s)))

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
