/* body.h - Acc, BODY, CBODY, Config, EBODY, MBODY, Mass, PBODY, PMass,
            PMom, PPos, PPot, PVel, PartType, PdPos_ds, PdVel_ds, Phase, Pos,
            Pot, Vel, body, bodyptr, cbody, cbodyptr, ebody, ebodyptr, 
            mbody, mbodyptr, pbody, pbodyptr */

/*
 *  body.h: for bodies and their xbody subsets, x out of {e,m,p,c}
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
/*-----------------------------------------------------------------------------
 *  body, bodyptr  --  contains dynamical information per particle in the form
 *                     of mass, and position, velocity and acceleration 
 *                     vectors, plus the gravitational potential.
 *                     The first entry gives the type of the structure,
 *                     identifying it as a body.
 *-----------------------------------------------------------------------------
 */
typedef struct
    {
    short  type;                           /* type of a body                 */
    real  phasespace[2][NDIM];             /* position and velocity vectors  */
    real  mass;                            /* mass of body                   */
    real  pot;                             /* potential at position of body  */
    real  acc[NDIM];                       /* acceleration vector, unless:   */
#ifdef REGULARIZATION                      /*   then the above is dQ/ds and: */
    real  dp_ds[NDIM];                     /* d/ds of Hamiltonian Momentum P */
#endif                                     /*   where s is pseudo-time       */
    } body, *bodyptr;

/*-----------------------------------------------------------------------------
 *  ebody, mbody, pbody, cbody  --  decremental subsets of a body
 *  ebodyptr, mbodyptr, pbodyptr, cbodyptr  --  their pointer equivalents
 *-----------------------------------------------------------------------------
 */
typedef struct
    {
    short  type;                              /* type of a ebody             */
    real  phasespace[2][NDIM];                /* position & velocity vectors */
    real  mass;                               /* mass                        */
    real  pot;                                /* potential                   */
    } ebody, *ebodyptr;             /* e for energy information              */

typedef struct
    {
    short  type;                              /* type of a mbody             */
    real  phasespace[2][NDIM];                /* position & velocity vectors */
    real  mass;                               /* mass                        */
    } mbody, *mbodyptr;             /* m for mass-point information          */

typedef struct
    {
    short  type;                              /* type of a pbody             */
    real  phasespace[2][NDIM];                /* position & velocity vectors */
    } pbody, *pbodyptr;             /* p for phase space information         */

typedef struct
    {
    short  type;                              /* type of a cbody             */
    real  configspace[NDIM];          /* only pos. vectors (or vel. or acc.) */
    } cbody, *cbodyptr;             /* c for configuration space information */

/*-----------------------------------------------------------------------------
 *  macros to extract individual components from a body, or one of its
 *  derivatives: ebody, mbody, pbody or cbody.
 *
 *  remarks:
 *       
 *   1: these macros are overloaded:  Pos()  and  Vel()  can be used on any
 *      type of body except cbodies; Mass()  on bodies, ebodies and mbodies; 
 *      Pot()  on bodies and ebodies;  Acc()  only on bodies.
 *      
 *   2: xbodies of smaller type than bodies (x being one of {e,m,p,c}) can
 *      be extracted from bodies in non-unique ways:
 *      
 *          ebodies    have a natural, unique meaning;
 *          mbodies    have a natural, unique meaning;
 *          pbodies    normally contain positions and velocities
 *                         (in which case Pos() and Vel() are natural),
 *                     but can also contain velocities and accelerations
 *                         (in which case Pos() and Vel() are still the
 *                          correct macros to use, although in this case Pos()
 *                          extracts the velocity and Vel() the acceleration)
 *          cbodies    normally contain positions,
 *                     but can also contain velocities,
 *                     and can also contain accelerations
 *                         (and therefore Config() is provided as a general
 *                          handle on cbodies)
 *      
 *   3: Phase(space vector) is formally equivalent to Pos(ition), but is 
 *      introduced to provide a recognizable handle for operations which 
 *      involve positions and velocities together as 2*NDIM - dimensional
 *      vectors.
 *
 *   4: PartType is chosen to identy the type of particle; Type would have
 *      been simpler, but this name may have been used elsewhere in a large
 *      environment, so PartType seems a more prudent choice.
 *  NOTE:  ptr could be a pointer of various types as described above:
 *         a "bodyptr" is valid for any of the macros below, while the 
 *         following pointers can be used for at least some of these macros:
 *         "ebodyptr, mbodyptr, pbodyptr, cbodyptr, nodeptr, cellptr".
 *-----------------------------------------------------------------------------
 */                        
#define  PartType(ptr)   ((ptr)->type)                      /* type: short   */
#define  Pos(ptr)        ((ptr)->phasespace[0])             /* type: realptr */
#define  Vel(ptr)        ((ptr)->phasespace[1])             /* type: realptr */
#define  Mass(ptr)       ((ptr)->mass)                      /* type: real    */
#define  Pot(ptr)        ((ptr)->pot)                       /* type: real    */
#define  Acc(ptr)        ((ptr)->acc)                       /* type: realptr */
#define  Phase(ptr)      ((real *)(ptr)->phasespace)        /* type: realptr */
#define  Config(ptr)     ((ptr)->configspace)               /* type: realptr */

/*-----------------------------------------------------------------------------
 *  encoding for different types of xbodies:
 *-----------------------------------------------------------------------------
 */
#define   BODY    010
#define   EBODY   011
#define   MBODY   012
#define   PBODY   013
#define   CBODY   014

#ifdef REGULARIZATION

/*-----------------------------------------------------------------------------
 *  macros to extract components from a regularized body:
 *
 *  in regularized representation, each regularized "particle" represents
 *  the information about the relative motion of a pair of unregularized 
 *  particles. To highlight the appearance of regularized particles in a piece
 *  of code, it is recommended to use different macros to extract the
 *  components of those particles: PPos() instead of Pos() , for example, 
 *  where  PPos  stands for  Pair Position . Note the use of  PMom()
 *  which stands for  Pair Momentum  instead of Vel(), since regularized
 *  coordinates are expressed in terms of positions and momenta instead of
 *  positions and velocities; they are Hamiltonian and not Lagrangian.
 *     A more drastic change is the replacement of the acceleration by separate
 *  derivatives of the PPos() and the PMom(), denoted by PdPos_ds() and
 *  PdMom_ds(), respectively. The reason is that in general the dynamic 
 *  information of a particle requires two vectors (position and momentum, say,
 *  or position and velocity) and the information of a particle and its time
 *  rate of change twice that: four NDIM-dimensional vectors. In the
 *  nonregularized representation Vel() just happens to play a double role,
 *  being the (d/dt) of Pos(), and letting Acc() be its (d/dt).
 *     Technically the definitions of PPos() and Pos() are identical, as are
 *  those of PMom() and Vel(), and of PdPos_ds and Acc(); the latter overlap
 *  is used when translating between regularized coordinates and unregularized
 *  coordinates. For simplicity, the same structure is used for the System()
 *  and the Regsystem() part of a state. Either system contains the same four
 *  fourdimensional vectors, with the difference that Pos(), Vel() and Acc()
 *  are used to access the information from a nonregularized system (in which
 *  case only the first three components of the vectors carry information),
 *  while PPos(), PMom(), PPosdot() and PMomdot() are used for a regularized
 *  system. Summing up:
 *  NOTE: CONVENTION: in regularized representation, the velocity and
 *                    acceleration slots are ALWAYS used to store ONLY momenta
 *                    and forces, as the macros below indicate;
 *                    in non-regularized representation, even in 4D with the
 *                    fourth component zero (or at least very small, if 
 *                    round-off errors sneak in), these slots are ALWAYS used
 *                    to store ONLY velocities and accelerations.
 *  NOTE:  ptr could be a pointer of type "bodyptr"; for some of the macros
 *         "ebodyptr", "mbodyptr", "pbodyptr" or  "cbodyptr" is a valid 
 *         alternative.
 *-----------------------------------------------------------------------------
 */                        
#define  PPos(ptr)       ((ptr)->phasespace[0])             /* type: realptr */
#define  PMom(ptr)       ((ptr)->phasespace[1])             /* type: realptr */
#define  PMass(ptr)      ((ptr)->mass)                      /* type: real    */
#define  PPot(ptr)       ((ptr)->pot)                       /* type: real    */
#define  PdPos_ds(ptr)   ((ptr)->acc)                       /* type: realptr */
#define  PdMom_ds(ptr)   ((ptr)->dp_ds)                     /* type: realptr */

#endif

/* endof: body.h */
