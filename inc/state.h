/* state.h - Ctrls, Diags, Regsystem, Specs, System, state, stateptr */

/*
 *  state.h: for a state, the basic building block in newton0
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
#include  "syst.h"
#include  "spec.h"
#include  "ctrl.h"
#include  "diag.h"

/*-----------------------------------------------------------------------------
 *  state, stateptr  --  contains the information describing the state of the
 *                       gravitational many-body simulation at a point in time.
 *                       the four substructures contain (in that order)
 *                       basic, local, global and diagnostics information.
 *-----------------------------------------------------------------------------
 */
typedef struct
    {
    systptr  system;          /* points to basic system block, containing    */
                              /* the body array, and other variables:        */
                              /* the number of bodies, the time, etc. which  */
                              /* are needed for integrating eqs. of motion.  */
    specptr  specs;           /* points to local specifications block, which */
                              /* governs details of orbit integration and    */
                              /* collects force calculation diagnostics.     */
    ctrlptr  ctrls;           /* points to global control block, which       */
                              /* orchestrates orbit integration, scheduling  */
                              /* and directing of input/output, etc.         */
    diagptr  diags;           /* points to diagnostics block                 */
#ifdef REGULARIZATION
    systptr  regsystem;       /* points to the regularized system            */
#endif
    } state, *stateptr;

/*-----------------------------------------------------------------------------
 *  macros to extract individual components from a state:
 *  NOTE:  ptr should be a pointer of type "stateptr":
 *-----------------------------------------------------------------------------
 */
#define  System(ptr)        ((ptr)->system)              /* type: systptr    */
#define  Specs(ptr)         ((ptr)->specs)               /* type: specptr    */
#define  Ctrls(ptr)         ((ptr)->ctrls)               /* type: ctrlptr    */
#define  Diags(ptr)         ((ptr)->diags)               /* type: diagptr    */
#ifdef REGULARIZATION
#  define  Regsystem(ptr)   ((ptr)->regsystem)           /* type: systptr    */
#endif

/* endof: state.h */
