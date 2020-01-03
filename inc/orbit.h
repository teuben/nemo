/*
 * ORBIT.H: structured binary file definitions for standard
 *      orbit paths
 *
 *      July 1987       Peter Teuben @ IAS, Princeton, NJ
 *
 *  13-Jul-87    V1.0  first attempt, still independant of lower level snapshot
 *  16-jul-87    V1.1  slight mod. in struct 
 *  28-jul-87    V2.0  orbit(5) changed
 *  18-dec-90    V2.1  save space allocate, for read_orbit's reallocing
 *   7-mar-92    V2.2  happy gcc2.0
 *  24-may-92    V3.0  added potential as part of an orbit
 *  10-mar-93          potential called a_potential, conflict with potential0.c
 *   6-dec-93    V3.1  added maxsteps, to be able to re-allocate
 *  22-feb-93    V3.1a ansi header
 *  11-apr-95    V3.1b no more ARGS, included more header files here
 *   1-mar-03    V3.3  added iom_err, errors in the integrals of motion
 *  25-jul-13    V4.0  added Key
 *  10-dec-2019  V4.1  added optional PHI and ACC (for BOOM)
 */

#include <filestruct.h>
#include <potential.h>
#include <history.h>

#define ORBIT_PHI

typedef struct {
        int     _ndim;
        int     _coordsys;
	int     _size;
        real    _mass;
        realptr _iom;
        realptr _iom_err;
        int     _key;
        int     _nsteps;
        int     _maxsteps;
        realptr _time;
        realptr _phase;
#ifdef ORBIT_PHI
        realptr _phi;
        realptr _acc;
#endif  
	a_potential _pot;
} orbit, *orbitptr;

#define Ndim(optr)      ((optr)->_ndim)
#define CoordSys(optr)  ((optr)->_coordsys)
#define Size(optr)      ((optr)->_size)
#define Masso(optr)     ((optr)->_mass)
#define Key(optr)       ((optr)->_key)
#define IOM(optr)       ((optr)->_iom)
#define IOMERR(optr)    ((optr)->_iom_err)
#define Nsteps(optr)    ((optr)->_nsteps)
#define MAXsteps(optr)  ((optr)->_maxsteps)
#define TimePath(optr)  ((optr)->_time)
#define PhasePath(optr) ((optr)->_phase)
#define PosPath(optr)   ((optr)->_phase)
#define VelPath(optr)   ((optr)->_phase+Ndim(optr))
#define PhiPath(optr)   ((optr)->_phi)
#define AccPath(optr)   ((optr)->_acc)
#define Potential(optr) ((optr)->_pot)

/*      a few dangerous (does not check for ndim) access functions */

#define I1(optr)        (*(IOM(optr)))
#define I2(optr)        (*(IOM(optr)+1))
#define I3(optr)        (*(IOM(optr)+2))
#define IE1(optr)       (*(IOMERR(optr)))
#define IE2(optr)       (*(IOMERR(optr)+1))
#define IE3(optr)       (*(IOMERR(optr)+2))
/*      OLD V1.x structure
#define Torb(optr,i)    *(TimePath(optr)+i)
#define Xorb(optr,i)    *(PhasePath(optr)+i)
#define Yorb(optr,i)    *(PhasePath(optr)+i+Nsteps(optr))
#define Zorb(optr,i)    *(PhasePath(optr)+i+2*Nsteps(optr))
#define Uorb(optr,i)    *(PhasePath(optr)+i+Ndim(optr)*Nsteps(optr))
#define Vorb(optr,i)    *(PhasePath(optr)+i+(Ndim(optr)+1)*Nsteps(optr))
#define Worb(optr,i)    *(PhasePath(optr)+i+(Ndim(optr)+2)*Nsteps(optr))
*/
/*      New V2.0 structure      */
#define Torb(optr,i)    (*(TimePath(optr)+(i)))
#define Posorb(o,i,j)   (*(PhasePath(o)+Ndim(o)*2*(i)+j))
#define Velorb(o,i,j)   (*(PhasePath(o)+Ndim(o)*2*(i)+Ndim(o)+j))
#define Xorb(optr,i)    (*(PhasePath(optr)+Ndim(optr)*2*(i)))
#define Yorb(optr,i)    (*(PhasePath(optr)+Ndim(optr)*2*(i)+1))
#define Zorb(optr,i)    (*(PhasePath(optr)+Ndim(optr)*2*(i)+2))
#define Uorb(optr,i)    (*(PhasePath(optr)+Ndim(optr)*2*(i)+Ndim(optr)))
#define Vorb(optr,i)    (*(PhasePath(optr)+Ndim(optr)*2*(i)+Ndim(optr)+1))
#define Worb(optr,i)    (*(PhasePath(optr)+Ndim(optr)*2*(i)+Ndim(optr)+2))
#define Porb(optr,i)    (*(PhiPath(optr)+(i)))
#define AXorb(optr,i)   (*(AccPath(optr)+Ndim(optr)*(i)))
#define AYorb(optr,i)   (*(AccPath(optr)+Ndim(optr)*(i)+1))
#define AZorb(optr,i)   (*(AccPath(optr)+Ndim(optr)*(i)+2))

#define PotName(optr)   (Potential(optr).name)
#define PotPars(optr)   (Potential(optr).pars)
#define PotFile(optr)   (Potential(optr).file)
/*
 * Item tags for Image components.
 */

#define OrbitTag                "Orbit"

#define ParametersTag         "Parameters"
#define     NdimTag             "Ndim"
#define     CoordSysTag         "CoordSys"
#define     MassTag             "Mass"
#define     KeyTag              "Key"
#define     IOMTag              "IOM"
#define     IOMERRTag           "IOM_error"

#define PotentialTag		"Potential"
#define     PotNameTag		"Name"
#define     PotParsTag		"Pars"
#define     PotFileTag		"File"

#define PathTag                 "Path"
#define     NstepsTag           "Nsteps"
#define     TimePathTag         "TimePath"
#define     PhasePathTag        "PhasePath"
#define     PosPathTag          "PosPath"
#define     VelPathTag          "VelPath"
#define     PhiPathTag          "PhiPath"
#define     AccPathTag          "AccPath"



/* external functions */
void write_orbit    ( stream, orbitptr );
int  read_orbit     ( stream, orbitptr * );
void free_orbit     ( orbitptr );
int  allocate_orbit ( orbitptr *, int, int );
void copy_orbit     ( orbitptr, orbitptr );
void list_orbit     ( orbitptr, double, double, int, string );
