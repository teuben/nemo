/* differentiate.c - BIGNUMBER, d_system, deriv_system, 
                     init_deriv_system, n_square_interaction, 
                     push_system, reg_interaction, tree_interaction */

/*
 *  differentiate.c:  contains helper functions for  integrate.c 
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */

#include  "newton0.h"

static void n_square_interaction(systptr sys, specptr specs);
static void tree_interaction(systptr sys, specptr specs);
static void reg_interaction(systptr sys, specptr specs);


/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |   differentiate.c :                                                 |  */
/*  |                                                                     |  */
/*  |           This file contains the dynamical part of newton0,         |  */
/*  |       where Newton's gravitational equations are implemented        |  */
/*  |       to find the accelerations in terms of positions, from the     |  */
/*  |       set of second order differential equations                    |  */
/*  |                                                                     |  */
/*  |            2            __          m                               |  */
/*  |           d    -->     \             j        /--> --> \            |  */
/*  |         _____   r   =   >   G -------------- (  r - r   )           |  */
/*  |             2    i     /__     |--> --> | 3   \  j   i /            |  */
/*  |         (dt)            j      | r - r  |                           |  */
/*  |                                |  i   j |                           |  */
/*  |                                                                     |  */
/*  |       or generalizations thereof, using a softening method.         |  */
/*  |                                                                     |  */
/*  |           In addition to the procedures presented here, the         |  */
/*  |       files  systemalgebra.c  and  systemtranslate.c  contain       |  */
/*  |       helper functions which manipulate systems (by invoking        |  */
/*  |       lower-level helperfunctions to manipulate parts of            |  */
/*  |       systems, such as accelerations, velocities and positions).    |  */
/*  |       These helper functions are the ingredients for the            |  */
/*  |       recipes given in the file  integrate.c  .                     |  */
/*  |       Those recipes for higher order integration schemes            |  */
/*  |       naturally fall into several categories, each of which         |  */
/*  |       require a different set of helper functions, to               |  */
/*  |       provide an interface between the integrations schemes         |  */
/*  |       and the force calculation routines presented below.           |  */
/*  |                                                                     |  */
/*  |           In addition to the force calculations, the present        |  */
/*  |       file contains intermediate procedures which can be            |  */
/*  |       viewed either as syntactic sugar dressing up the general      |  */
/*  |       helper functions for specific usage (cf. d_system() ),        |  */
/*  |       or they can be viewed as providing an extra conceptual        |  */
/*  |       barrier between the recipes in  integrate.c  and the          |  */
/*  |       ingredients in  systemalgebra.c  and  systemconversion.c .    |  */
/*  |       Time will tell in which direction this file will evolve.      |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                              PART I                                 |  */
/*  |                                                                     |  */
/*  |       contains a functions to update the dynamical state            |  */
/*  |       of a system,using the positions to set the                    |  */
/*  |       accelerations and the potentials.                             |  */
/*  |           These functions use helperfunctions from the file         |  */
/*  |       systemalgebra.c , but not from the file  systemconversion.c . |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  init_deriv_system  --  initial computation of the system derivative,
 *                         before the first orbit calculations are started up;
 *                         see deriv_system() on the next page for a
 *                         description. init_deriv_system() is introduced as a
 *                         separate procedure because regularization requires
 *                         a different treatment at startup (3D energy
 *                         calculation is required to initialize 4D equations
 *                         of motion).
 *-----------------------------------------------------------------------------
 */
void  init_deriv_system(sys, specs)
systptr  sys;
specptr specs;
    {
#ifndef TREE
# ifndef REGULARIZATION
    n_square_interaction(sys, specs);
# endif
#endif
#ifdef TREE
    tree_interaction(sys, specs);
#endif
#ifdef REGULARIZATION                      /* the following line is the only */
    n_square_interaction(sys, specs);      /* difference with deriv_system() */
#endif
    }

/*-----------------------------------------------------------------------------
 *  deriv_system  --  computes the system derivative in the form of the
 *                    acceleration and potential for all particles in the
 *                    system, using their positions and Newtons law of gravity
 *                    modified by a softening method to dampen strong 
 *                    interactions during close encounters between particle
 *                    pairs (such encounters should be dampened since they are
 *                    unphysical if our particles represent only a small
 *                    sample of all points in a real system such as a galaxy).
 *                    accepts: sys: pointer to a nbody system in standard
 *                                  newton0 form;
 *                           specs: pointer to the specifications block,
 *                                  which contains that part of the control and
 *                                  diagnostics information which is
 *                                  intertwined with the force calculation.
 *                    effect: sets the potentials and accelerations in the
 *                            bodies in system "sys" to the computed values
 *                            (or dQ/ds and dP/ds etc. for regularized systems)
 *-----------------------------------------------------------------------------
 */
void  deriv_system(sys, specs)
systptr  sys;
specptr specs;
    {
#ifndef TREE
# ifndef REGULARIZATION
    n_square_interaction(sys, specs);
# endif
#endif
#ifdef TREE
    tree_interaction(sys, specs);
#endif
#ifdef REGULARIZATION
    reg_interaction(sys, specs);
#endif
    }

#define  BIGNUMBER  1.0e20

/*-----------------------------------------------------------------------------
 *  n_square_interaction  --  calculate a system derivative by explicitly
 *                            evaluating all interactions between all particles
 *                            see  deriv_system()  above.
 *-----------------------------------------------------------------------------
 */
local void  n_square_interaction(sys, specs)
systptr  sys;
specptr specs;
    {
    int  npart;                  /* number of bodies                         */
    real  r0_soft;               /* softening length                         */
    real  separ_ji[NDIM];        /* separation vector from i to j            */
    real  pot_helper;            /* auxiliary variable to store soft 1/R     */
    real  acc_helper;            /* auxiliary variable to store soft 1/R^3   */
    real  acc_i_helper;          /* auxiliary variable to store soft Mj/R^3  */
    real  acc_j_helper;          /* auxiliary variable to store soft -Mi/R^3 */
    real  acc_ji[NDIM];          /* acceleration vector from j on i          */
    real  acc_ij[NDIM];          /* acceleration vector from i on j          */
    real  inv_soft_pair_min;     /* minimum pair separation                  */
    bodyptr  bodies;             /* pointer to first body                    */
    bodyptr  body_i, body_j;     /* pairs to compute Newtonian gravity       */
#if 0
    proc  softener();            /* returns a pointer to the softening type  */
#endif
    proc  soften_potential;      /* computes softened potential and M/R^3 :  */
                                 /* a procedure, short for: real (*i_s)();   */


    npart = Nbody(sys);
    bodies = Bodies(sys);

    r0_soft = Softparam(specs);
    soften_potential = softener(Softfocus(specs));

    for (body_i = bodies; body_i - bodies < npart; body_i++)
	{
	CLRV(Acc(body_i));
	Pot(body_i) = 0.0;
	}

    inv_soft_pair_min = 1.0 / BIGNUMBER;

    for (body_i = bodies; body_i - bodies < npart - 1; body_i++)
        for (body_j = body_i + 1;  body_j - bodies < npart; body_j++)
            {
	    SUBV(separ_ji, Pos(body_j), Pos(body_i));
            (*soften_potential)(separ_ji, r0_soft, &pot_helper, &acc_helper);
            Pot(body_i) -= Mass(body_j) * pot_helper;
            Pot(body_j) -= Mass(body_i) * pot_helper;
            acc_i_helper = Mass(body_j) * acc_helper;
            acc_j_helper = -Mass(body_i) * acc_helper;
	    MULVS(acc_ji, separ_ji, acc_i_helper);
	    MULVS(acc_ij, separ_ji, acc_j_helper);
	    INCADDV(Acc(body_i), acc_ji);
	    INCADDV(Acc(body_j), acc_ij);

	    inv_soft_pair_min = MAX(inv_soft_pair_min, pot_helper);
	    }

    SPECmin_pair(specs) = MIN(SPECmin_pair(specs), 1.0/inv_soft_pair_min);
    }

#ifdef TREE

/*-----------------------------------------------------------------------------
 *  tree_interaction  --  calculate a system derivative by using the tree
 *                        approximation to evaluate the interactions between
 *                        the particles; see also  deriv_system()  above.
 *                        Litt.: J. Barnes and P. Hut, Nature 324, 446 (1986).
 *-----------------------------------------------------------------------------
 */
local void  tree_interaction(sys, specs)
systptr  sys;
specptr  specs;
    {
    mk_tree(sys);                                  /* in transformtree.c     */

    tree_tree_interaction(sys, specs);             /* in differentiatetree.c */
    }

#endif

#ifdef REGULARIZATION

/*-----------------------------------------------------------------------------
 *  reg_interaction  --  calculate a system derivative for a regularized system
 *                       note: only one choice for now; more to come
 *                       Litt.: D.C. Heggie, Celestial Mechanics 10, 217 (1974)
 *                        & S. Mikkola, Mon. Not. R. astr. Soc. 215, 171 (1985)
 *-----------------------------------------------------------------------------
 */
local void  reg_interaction(sys, specs)
systptr  sys;
specptr  specs;
    {
    heggie_mikkola_equations_of_motion(sys, specs);
    }

#endif

/*===========================================================================*/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                             PART II                                 |  */
/*  |                                                                     |  */
/*  |        contains functions which combine individual ingredients      |  */
/*  |        from the files  systemalgebra.c  and  systemconversion.c     |  */
/*  |        to define composite operations which are natural             |  */
/*  |        higher-order ingredients for use in procedures in the        |  */
/*  |        file  integrate.c .                                          |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  d_system  --  constructs the differential of the phase space information.
 *                accepts: sys: pointer to a nbody system in standard
 *                              newton0 form;
 *                       specs: pointer to the specifications block,
 *                              which contains that part of the control and
 *                              diagnostics information which is
 *                              intertwined with the force calculation.
 *                          ds: the pseudo-time differential, with which the
 *                              system derivative is multiplied to obtain
 *                              "d_sys"; for non-regularized system, ds is
 *                              simply the time differential dt, but for a
 *                              regularized system ds is the time-like
 *                              parameter introduced to regularize the
 *                              equations of motion.
 *-----------------------------------------------------------------------------
 */
psystptr  d_system(sys, specs, ds)
systptr  sys;
specptr  specs;
real  ds;
    {
    psystptr  d_sys;              /* will contain differential of system     */

    d_sys = sel_phasedot(sys);

    scale_psystem(d_sys, ds);

    return(d_sys);
    }

/*-----------------------------------------------------------------------------
 *  push_system  --  evolves the positions and velocities in a gravitational
 *                   many-body system "sys" by adding the incremental system
 *                   given by "d_sys".
 *                   accepts: sys: pointer to a nbody system in standard
 *                                 newton0 form;
 *                          specs: pointer to the specifications block,
 *                                 which contains that part of the control and
 *                                 diagnostics information which is
 *                                 intertwined with the force calculation.
 *-----------------------------------------------------------------------------
 */
void  push_system(sys, d_sys)
systptr  sys;
psystptr  d_sys;
    {
    annex_phase(sys, d_sys);                  /* in file  systemconversion.c */
    }

/*===========================================================================*/

/* endof: differentiate.c */
