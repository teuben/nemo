/* syst.h - Bodies, CBodies, EBodies, MBodies, Massmatrix, Nbody, Ncell, 
            Nmaxcell, PBodies, Regenergy, Reghamiltonian, Reglagrangian, Root, 
            Potmincorner, Potsizesq, Tnow, csyst, csystptr, esyst, esystptr,
            msyst, msystptr, psyst, psystptr, syst, systptr */

/*
 *  syst.h: for "syst", a substructure of a "state".
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
#include  "body.h"
#ifdef EXTRAPOLATION
#  include  "bodyext.h"
#endif
#ifdef TREE
#  include  "node.h"
#endif

/*-----------------------------------------------------------------------------
 *  syst, systptr  --  contains dynamical information of the system in the form
 *                     of an array of bodies, the length of the array, the time
 *                     and other variables which are essential for the force
 *                     calculations.
 *                     note: the extra components needed for tree or
 *                           regularized calculations are not reproduced in
 *                           xsystems (x out of {e,m,p,c}) on the next page;
 *                           this could be changed in the future if needed.
 *                     note: in the REGULARIZATION version, the unregularized
 *                           center-of-mass position and velocity are taken
 *                           out, and stored separately; the need for doing
 *                           this follows from the fact that the regularized 
 *                           description assumes that the system is descibed
 *                           in its center-of-mass frame.
 *-----------------------------------------------------------------------------
 */
typedef struct
    {
    real  tnow;                /* current time for all bodies                */
    int  nbody;                /* number of bodies                           */
    bodyptr  bodies;           /* points to first body in array              */
#ifdef EXTRAPOLATION
    orbitsegmentptr  orbitsegments; /* points to first orbitsegment in array */
#endif
#ifdef TREE
    int  ncell;                /* number of cells in use                     */
    nodeptr  root;             /* origin of the tree of cells and bodies.    */
		               /* root is declared as a node pointer to      */
                               /* allow the special case of a tree with      */
                               /* only one body attached                     */
    int  nmaxcell;             /* total number of cells allocated            */
    real  potmincorner[NDIM];  /* lower-left corner of outer coordinate box  */
    real  potsizesq; 	       /* square of side-length of integerized       */
                     	       /* outer coordinate box                       */
#endif
#ifdef REGULARIZATION
    real  com_position[NDIM];  /* position of center of mass                 */
    real  com_velocity[NDIM];  /* velocity of center of mass                 */
    realptr  massmatrix;       /* contains coefficients for the Hamiltonian  */
    real  regenergy;           /* intial total energy of regularized system  */
    real  reglagrangian;       /* Lagrangian of regularized system           */
    real  reghamiltonian;      /* Hamiltonian of regularized system          */
#endif
    } syst, *systptr;

/*-----------------------------------------------------------------------------
 *  macros to extract individual components from a system:
 *  note: Tnow always denotes the real, unregularized time, when part of
 *        System(some_state)  as well as of  Regsystem(some_state) .
 *  NOTE:  ptr should be a pointer of type "systptr" for all macros but the
 *         first two;  Tnow()  and  Nbody()  can be given a pointer of type
 *         "systptr", "esystptr", "msystptr", "psystptr", or "csystptr".
 *-----------------------------------------------------------------------------
 */
#define  Tnow(ptr)              ((ptr)->tnow)              /* type: real     */
#define  Nbody(ptr)             ((ptr)->nbody)             /* type: int      */
#define  Bodies(ptr)            ((ptr)->bodies)            /* type: bodyptr  */
#ifdef EXTRAPOLATION
#  define  Orbitsegments(ptr)   ((ptr)->orbitsegments)     /* type:          */
#endif                                                     /* orbitsegmentptr*/
#ifdef TREE
#  define  Ncell(ptr)           ((ptr)->ncell)             /* type: int      */
#  define  Root(ptr)            ((ptr)->root)              /* type: nodeptr  */
#  define  Nmaxcell(ptr)        ((ptr)->nmaxcell)          /* type: int      */
#  define  Potmincorner(ptr)    ((ptr)->potmincorner)      /* type: realptr  */
#  define  Potsizesq(ptr)       ((ptr)->potsizesq)         /* type: real     */
#endif
#ifdef REGULARIZATION
#  define  Com_Pos(ptr)         ((ptr)->com_position)      /* type: realptr  */
#  define  Com_Vel(ptr)         ((ptr)->com_velocity)      /* type: realptr  */
#  define  Massmatrix(ptr)      ((ptr)->massmatrix)        /* type: realptr  */
#  define  Regenergy(ptr)       ((ptr)->regenergy)         /* type: real     */
#  define  Reglagrangian(ptr)   ((ptr)->reglagrangian)     /* type: real     */
#  define  Reghamiltonian(ptr)  ((ptr)->reghamiltonian)    /* type: real     */
#endif

/*-----------------------------------------------------------------------------
 *  esyst, esystptr  --  contains a esystem in the form of an array of ebodies,
 *                       the length of the array and the time.
 *-----------------------------------------------------------------------------
 */
typedef struct
    {
    real  tnow;                /* current time for all ebodies               */
    int  nbody;                /* number of ebodies                          */
    ebodyptr  ebodies;         /* points to first ebody in array             */
#ifdef REGULARIZATION
    real  com_position[NDIM];  /* position of center of mass                 */
    real  com_velocity[NDIM];  /* velocity of center of mass                 */
#endif
    } esyst, *esystptr;

/*-----------------------------------------------------------------------------
 *  msyst, msystptr  --  contains a msystem in the form of an array of mbodies,
 *                       the length of the array and the time.
 *-----------------------------------------------------------------------------
 */
typedef struct
    {
    real  tnow;                /* current time for all mbodies               */
    int  nbody;                /* number of mbodies                          */
    mbodyptr  mbodies;         /* points to first mbody in array             */
#ifdef REGULARIZATION
    real  com_position[NDIM];  /* position of center of mass                 */
    real  com_velocity[NDIM];  /* velocity of center of mass                 */
#endif
    } msyst, *msystptr;

/*-----------------------------------------------------------------------------
 *  psyst, psystptr  --  contains a psystem in the form of an array of pbodies,
 *                       the length of the array and the time.
 *-----------------------------------------------------------------------------
 */
typedef struct
    {
    real  tnow;                /* current time for all pbodies               */
    int  nbody;                /* number of pbodies                          */
    pbodyptr  pbodies;         /* points to first pbody in array             */
#ifdef REGULARIZATION
    real  com_position[NDIM];  /* position of center of mass                 */
    real  com_velocity[NDIM];  /* velocity of center of mass                 */
#endif
    } psyst, *psystptr;

/*-----------------------------------------------------------------------------
 *  csyst, csystptr  --  contains a csystem in the form of an array of cbodies,
 *                       the length of the array and the time.
 *-----------------------------------------------------------------------------
 */
typedef struct
    {
    real  tnow;                /* current time for all cbodies               */
    int  nbody;                /* number of cbodies                          */
    cbodyptr  cbodies;         /* points to first cbody in array             */
    } csyst, *csystptr;

/*-----------------------------------------------------------------------------
 *  macros to extract individual components from a system: time and number of
 *  bodies can be extracted using  Tnow()  and  Nbody()  introduced above;
 *  the following macros are identical to  Bodies()  used above, but are here
 *  introduced to accentuate the use of different xbodies, x out of {e,m,p,c}
 *-----------------------------------------------------------------------------
 */
#define  EBodies(ptr)         ((ptr)->ebodies)           /* type: ebodyptr   */
#define  MBodies(ptr)         ((ptr)->mbodies)           /* type: mbodyptr   */
#define  PBodies(ptr)         ((ptr)->pbodies)           /* type: pbodyptr   */
#define  CBodies(ptr)         ((ptr)->cbodies)           /* type: cbodyptr   */

/* endof: syst.h */
