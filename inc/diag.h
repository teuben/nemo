/* diag.h - CURRENT, DIAGde_nsteps, DIAGekin, DIAGemax_err, DIAGepot, DIAGetot,
            DIAGmin_inert, DIAGmin_pair, DIAGn2bcalc, DIAGnbccalc, DIAGnfcalc,
            DIAGnsteps, DIAGtime, DRIFT, INCREMENT, INITIAL, PREVIOUS,
            RELDRIFT, RELINCREMENT, diag, diagptr */

/*
 *  diag.h: for "diag", a substructure of a "state".
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
/*-----------------------------------------------------------------------------
 *  diag, diagptr  --  contains a set of diagnostics of the global system,
 *                     such as the total amount of different forms of energy,
 *                     angular momentum, and possibly other global quantities
 *                     as well as recent and cumulative changes in those
 *                     quantities (providing important handles on the accuracy
 *                     of the integration).
 *                     NOTE: the diag_ncalc[3] is left in here for possible
 *                           future use; their function is currently taken 
 *                           over by their equivalents in "spec.h" . 
 *  note: this header file, too, should be structered using "Number_of_..."
 *        arrays as used in  ctrl.h  and  specs.h ; to be done.
 *-----------------------------------------------------------------------------
 */
typedef struct
    {
    int   diag_nsteps;       /* total number of integration steps performed  */
    int   diag_de_nsteps;    /* number steps since last energy check         */
    real  diag_time[3];      /* time of begin, previous and present check    */
    real  diag_etot[3];      /* total energy at each of these times          */
    real  diag_ekin[3];      /* kinetic energy at each of these times        */
    real  diag_epot[3];      /* potential energy at each of these times      */
    real  diag_emax_err;     /* maximum energy error through orbit history   */
    real  diag_min_inert;    /* minimum past moment of inertia of the system */
    real  diag_min_pair;     /* minimum pair distance through orbit history  */
#ifdef TREE
    int  diag_ncalc[3];      /* statistics for force calculations body/cells */
#endif
    } diag, *diagptr;

/*-----------------------------------------------------------------------------
 *  macros to extract individual (sub)components from a diag, which can be
 *  invoked as, for example: 
 *
 *  DIAGekin(ptr)[PREVIOUS]   for the previous value of the kinetic energy;
 *
 *  NOTE:  ptr should be a pointer of type "diagptr":
 *-----------------------------------------------------------------------------
 */                        
#define  INITIAL    0
#define  PREVIOUS   1
#define  CURRENT    2

#define  DIAGnsteps(ptr)     ((ptr)->diag_nsteps)           /* type: int     */
#define  DIAGde_nsteps(ptr)  ((ptr)->diag_de_nsteps)        /* type: int     */
#define  DIAGtime(ptr)       ((ptr)->diag_time)             /* type: realptr */
#define  DIAGetot(ptr)       ((ptr)->diag_etot)             /* type: realptr */
#define  DIAGekin(ptr)       ((ptr)->diag_ekin)             /* type: realptr */
#define  DIAGepot(ptr)       ((ptr)->diag_epot)             /* type: realptr */
#define  DIAGemax_err(ptr)   ((ptr)->diag_emax_err)         /* type: real    */
#define  DIAGmin_inert(ptr)  ((ptr)->diag_min_inert)        /* type: real    */
#define  DIAGmin_pair(ptr)   ((ptr)->diag_min_pair)         /* type: real    */
#ifdef TREE
#  define  DIAGnfcalc(ptr)   ((ptr)->diag_ncalc[0])         /* type: int     */
#  define  DIAGn2bcalc(ptr)  ((ptr)->diag_ncalc[1])         /* type: int     */
#  define  DIAGnbccalc(ptr)  ((ptr)->diag_ncalc[2])         /* type: int     */
/*
 *   DIAGnfcalc():	  number of n-on-1 force calculations
 *   DIAGn2bcalc():	  number of 2-body force calculations
 *   DIAGnbccalc():	  number of body-cell force calculations
 */
#endif

/*-----------------------------------------------------------------------------
 *  macros to extract absolute or relative changes from one of the arrays in
 *  a "diag" structure, which can be invoked as, for example:
 *
 *  INCREMENT(DIAGekin(ptr))  for the incremental change of the kinetic energy;
 *  RELDRIFT(DIAGepot(ptr))   for the relative total change of the pot. energy.
 *
 *  NOTE:  exp should be one of the above expressions, such as "DIAGekin(ptr)"
 *-----------------------------------------------------------------------------
 */
#define DRIFT(exp)      (((exp)[CURRENT]) - ((exp)[INITIAL]))  /* type: real */
#define INCREMENT(exp)  (((exp)[CURRENT]) - ((exp)[PREVIOUS])) /* type: real */
#define RELDRIFT(exp)       ((((exp)[CURRENT]) - ((exp)[INITIAL]))  /         \
                             ((exp)[INITIAL]))                 /* type: real */
#define RELINCREMENT(exp)   ((((exp)[CURRENT]) - ((exp)[PREVIOUS])) /         \
                             ((exp)[PREVIOUS]))                /* type: real */

/* endof: diag.h */
