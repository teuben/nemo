/* statealgebra.c - */

/*
 *  statealgebra.c:  for algebraic operations in newton0, other than on bodies
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */

#include  "newton0.h"

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |   statealgebra.c :                                                  |  */
/*  |                                                                     |  */
/*  |           This file contains procedures for operations on states,   |  */
/*  |           and on their substructures other than bodies and systems. |  */
/*  |                                                                     |  */
/*  |           THIS PART TO BE REWRITTEN                                 |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  mk_state  --  allocates memory for a many-body state in standard newton0
 *                form.
 *                returns: new_state: a pointer to the new state.
 *                NOTE: memory is allocated for the new "state" structure,
 *                      as well as for the generic "info", "control" and "diag"
 *                      substructures, which have a fixed size.
 *                      no memory is allocated for other substructures of a
 *                      state such as bodies and cells since their size is 
 *                      not known  a priori (also, these variable-size
 *                      structures all have there own "mk_whatever()"
 *                      memory allocation procedures which form the preferred
 *                      interface).
 *                      no memory is allocated for subsubstructures such as
 *                      the mass matrix in regularized coordinates, which is
 *                      a substructure of the "info" substructure of a "state",
 *                      for the same reasons.
 *-----------------------------------------------------------------------------
 */
stateptr  mk_state()
    {
    stateptr  new_state;

    new_state = (stateptr) malloc(sizeof(state));
    if (new_state == NULL)
	error("mk_state: not enough memory left for a new state\n");

    System(new_state) = NULL;               /* good form, and safer this way */
    Specs(new_state) = mk_specs();
    Ctrls(new_state) = mk_ctrls();
    Diags(new_state) = mk_diags();
#ifdef REGULARIZATION
    Regsystem(new_state) = NULL;            /* good form, and safer this way */
#endif

    return(new_state);
    }

/*-----------------------------------------------------------------------------
 *  rm_state  --  frees memory for a many-body state in standard newton0 form.
 *                accept: old_state: pointer to state which is to be removed.
 *                NOTE: see the limitations currently present when running
 *                      a tree or a regularized calculation, as described in
 *                      rm_system()  in the file  systemalgebra.c . BEWARE!
 *-----------------------------------------------------------------------------
 */
stateptr  rm_state(old_state)
stateptr  old_state;
    {
    rm_system(System(old_state));             /* in  systemalgebra.c */
    rm_specs(Specs(old_state));               /* see below */
    rm_ctrls(Ctrls(old_state));               /* see below */
    rm_diags(Diags(old_state));               /* see below */

    free(old_state);
    }

/*-----------------------------------------------------------------------------
 *  mk_specs  --  allocates memory for a spec structure in standard
 *                newton0 form.
 *                returns: new_specs: a pointer to the new spec structure
 *-----------------------------------------------------------------------------
 */
specptr  mk_specs()
    {
    specptr  new_specs;

    new_specs = (specptr) malloc(sizeof(spec));
    if (new_specs == NULL)
	error("mk_specs: not enough memory left for a new spec\n");

    return(new_specs);
    }

/*-----------------------------------------------------------------------------
 *  rm_specs  --  frees memory for a spec structure in standard newton0 form.
 *                accepts: old_specs: a pointer to an old spec structure
 *-----------------------------------------------------------------------------
 */
void  rm_specs(old_specs)
specptr  old_specs;
    {
    int  i;

    for (i = 0; i < Number_of_methods; i++)
	free(Methods(old_specs)[i]);

    free(old_specs);
    }

/*-----------------------------------------------------------------------------
 *  mk_ctrls  --  allocates memory for a control structure in standard
 *                   newton0 form.
 *                   returns: new_ctrls: a pointer to the new control 
 *                                          structure.
 *-----------------------------------------------------------------------------
 */
ctrlptr  mk_ctrls()
    {
    ctrlptr  new_ctrls;

    new_ctrls = (ctrlptr) malloc(sizeof(ctrl));
    if (new_ctrls == NULL)
	error("mk_ctrls: not enough memory left for a new control\n");

    return(new_ctrls);
    }

/*-----------------------------------------------------------------------------
 *  rm_ctrls  --  frees memory for a ctrl structure in standard newton0 form.
 *                accepts: old_ctrls: a pointer to an old ctrl structure
 *-----------------------------------------------------------------------------
 */
void  rm_ctrls(old_ctrls)
ctrlptr  old_ctrls;
    {
    int  i;

    for (i = 0; i < Number_of_filenames; i++)
	free(Filenames(old_ctrls)[i]);
    for (i = 0; i < Number_of_messages; i++)
	free(Messages(old_ctrls)[i]);

    free(old_ctrls);
    }

/*-----------------------------------------------------------------------------
 *  mk_diags  --  allocates memory for a diag structure in standard
 *                newton0 form.
 *                returns: new_diags: a pointer to the new diag structure
 *-----------------------------------------------------------------------------
 */
diagptr  mk_diags()
    {
    diagptr  new_diags;

    new_diags = (diagptr) malloc(sizeof(diag));
    if (new_diags == NULL)
	error("mk_diags: not enough memory left for a new diag\n");

    return(new_diags);
    }

/*-----------------------------------------------------------------------------
 *  rm_diags  --  frees memory for a diag structure in standard newton0 form.
 *                accepts: old_diags: a pointer to an old diag structure
 *-----------------------------------------------------------------------------
 */
void  rm_diags(old_diags)
diagptr  old_diags;
    {
    free(old_diags);
    }

#ifdef  FIND_SOME_TIME_TO_FIX_THIS

/*-----------------------------------------------------------------------------
 *  cp_state  --  makes a copy of a full state in standard newton0 form.
 *                accepts: old_state: pointer to a state;
 *                returns: new_state: pointer to a new copy of the state.
 *  PRELIMINARY  !!
 *  PRELIMINARY  !!
 *  PRELIMINARY  !!
 *  PRELIMINARY  !!
 *-----------------------------------------------------------------------------
 */
stateptr  cp_state(old_state)
stateptr  old_state;
    {
    stateptr  new_state;
#ifdef REGULARIZATION
    realptr  cp_massmatrix();
#endif
#ifdef TREE
    cellptr  cp_cells();
#endif

    new_state = mk_state();

    Nbody(new_state) = Nbody(old_state);
    Bodies(new_state) = cp_system(Bodies(old_state), Nbody(old_state));
    Ctrls(new_state) = cp_ctrls(Ctrls(old_system));
    Diags(new_state) = cp_diags(Diags(old_system));
#ifdef REGULARIZATION
    Nregbody(new_state) = Nregbody(old_state);
    Regbodies(new_state) = cp_system(Regbodies(old_state),
                                                         Nregbody(old_state));
    Massmatrix(new_state) = cp_massmatrix(Massmatrix(old_state),
                                                         Nbody(old_state));
#endif
#ifdef TREE
    Root(new_state) = cp_cells(Root(old_state),
                                       Nmaxcell(old_state), Ncell(old_state));
    Ncell(new_state) = Ncell(old_state);
    Nmaxcell(new_state) = Nmaxcell(old_state);
    SETV(Rootrmin(new_state), Rootrmin(old_state));
    Rootrsize(new_state) = Rootrsize(old_state);
#endif

    return(new_state);
    }

#endif

/* endof: statealgebra.c */
