/* differentiatetree.c - body_node_interaction, body_tree_interaction, 
                         celldivider, constant_opening_angle, 
                         tree_tree_interaction, treewalk */

/*
 *  differentiatetree.c: routines to compute gravity. 
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 *
 *        adapted from hackcode versions by Josh Barnes, 1985-7,
 *        with the quadrupole extension by Lars Hernquist, 1986.
 *      
 *      This code needs -DTREE to compile, and optionally -DQUADPOLE
 */

#include "newton0.h"

typedef bool (*bproc)();

static void treewalk(bodyptr bp, nodeptr np, specptr specs, proc workptr);
static void body_node_interaction(register bodyptr bp, register nodeptr np, real r0_soft, proc softening_method);
static bproc celldivider(string type_of_celldivision);
static bool constant_opening_angle(register bodyptr bp, register nodeptr np, real tolsq);

/*-----------------------------------------------------------------------------
 *  tree_tree_interaction  --  evaluate grav field for each particle in "sys".
 *-----------------------------------------------------------------------------
 */
void  tree_tree_interaction(sys, specs)
systptr  sys;
specptr  specs;
    {
    bodyptr  body_i;

    for (body_i = Bodies(sys); body_i - Bodies(sys) < Nbody(sys); body_i++)
	body_tree_interaction(body_i, Root(sys), specs);

    Nfcalc(specs) += Nbody(sys);             /* number of n-on-1 force calc. */
    }

/*-----------------------------------------------------------------------------
 *  body_tree_interaction  --  evaluate grav field at a given particle.
 *-----------------------------------------------------------------------------
 */
void  body_tree_interaction(bp, root, specs)
bodyptr  bp;
nodeptr  root;
specptr  specs;
    {

    Pot(bp) = 0.0;				       /* clear potential    */
    CLRV(Acc(bp));	   			       /* clear acceleration */
    treewalk(bp, root, specs, body_node_interaction); /* recursively compute */
    }

/*-----------------------------------------------------------------------------
 *  treewalk  --   walk the tree opening cells too close to a given point.
 *-----------------------------------------------------------------------------
 */
local void  treewalk(bp, np, specs, workptr)
bodyptr  bp;
nodeptr np;                     /* pointer into body-tree */
specptr  specs;
proc workptr;				/* routine to do calculation */
    {
    register nodeptr *npp;         /* npp for  node-pointer-pointer          */
    bproc  celldivision_method;    /* one particular celldivision method     */
    proc  softening_method;        /* computes softened potential and M/R^3. */

    celldivision_method = celldivider(Celldivisionmethod(specs));
    softening_method = softener(Softfocus(specs));
#ifdef QUADPOLE
    if (! streq("plummer", Softfocus(specs)))
	error("treewalk: quadrupole softening of \"%s\" not implemented\n",
              Softfocus(specs));
#endif

    if ((*celldivision_method)(bp, np, Tolsqparam(specs)))       /* open np? */
        {                                            /* loop over sub-cells  */
        for (npp = Subptr(np); npp - Subptr(np) < NSUB; npp++)
            if (*npp != NULL)                        /* does this one exist? */
                treewalk(bp, *npp, specs, workptr);          /*  then use it */
        } 
    else if (np != (nodeptr)bp)                      /* not to be skipped?   */
        {                                                    /*  then use it */
        (*workptr)(bp, np, Softparam(specs), softening_method);
	if (PartType(np) == BODY)
	    (N2bcalc(specs))++;			     /* count body-body int. */
	else
	    (Nbccalc(specs))++;			     /* count body-cell int. */
        }
    }

/*-----------------------------------------------------------------------------
 *  body_node_interaction  --  compute a single 2-body interaction.
 *-----------------------------------------------------------------------------
 */
local void  body_node_interaction(bp, np, r0_soft, softening_method)
register bodyptr bp;                        /* body whose force to calculate */
register nodeptr np;                        /* body or cell to interact with */
real  r0_soft;                              /* softening parameter           */
proc  softening_method;             /* computes softened potential and M/R^3 */
    {
    real  mor3;
    real  acc_body_node[NDIM], quaddr[NDIM];
    real  dr5inv, phiquad, drquaddr;
    real  pot_helper;            /* auxiliary variable to store soft 1/R     */
    real  acc_helper;            /* auxiliary variable to store soft 1/R^3   */
    real  separ[NDIM];  	 /* separation vector  body --> node .       */

    SUBV(separ, Pos(np), Pos(bp));
    (*softening_method)(separ, r0_soft, &pot_helper, &acc_helper);

    Pot(bp) -= Mass(np) * pot_helper;
    mor3 = Mass(np) * acc_helper;
    MULVS(acc_body_node, separ, mor3);
    INCADDV(Acc(bp), acc_body_node);             /* add to net acceleration  */
#ifdef QUADPOLE
    if(PartType(np) == CELL)                     /* if cell, add quad. term  */
        {
        dr5inv = acc_helper*pot_helper*pot_helper;        /*  separ ** (-5)  */
        MULMV(quaddr, Quad(np), separ);            /* form Q * separ         */
        DOTVP(drquaddr, separ, quaddr);            /* form separ * Q * separ */
        phiquad = -0.5 * dr5inv * drquaddr;        /* quad. part of poten.   */
        Pot(bp) += phiquad;                        /* increment potential    */
        phiquad = 5.0*phiquad*pot_helper*pot_helper;    /* save for acc.     */
        MULVS(acc_body_node, separ, phiquad);      /* components of acc.     */
        INCSUBV(Acc(bp), acc_body_node);           /* increment              */
        INCMULVS(quaddr, dr5inv);   
        INCSUBV(Acc(bp), quaddr);                  /* acceleration           */
	}
#endif
    }

/*-----------------------------------------------------------------------------
 *  celldivider  --  returns the pointer  ptr  to the function which decides
 *                   whether cells should be subdivided
 *-----------------------------------------------------------------------------
 */
local bproc  celldivider(type_of_celldivision)
string  type_of_celldivision;
    {
    if (streq("constant_theta", type_of_celldivision))
	return( constant_opening_angle );
    else
	error("celldivider: %s not implemented as a type of cell division\n",
                                                         type_of_celldivision);
    }

/*-----------------------------------------------------------------------------
 *  constant_opening_angle --  decide if a node should be opened,
 *                             using a constant opening angle criterion
 *-----------------------------------------------------------------------------
 */
local bool  constant_opening_angle(bp, np, tolsq)
register bodyptr bp;                     /* body of observer               */
register nodeptr np;                     /* body/cell to be tested         */
real tolsq;                              /* tolerance parameter squared    */
    {
    real  separ[NDIM];
    real  distsq;

    if (PartType(np) == BODY)                       /* at tip of tree?       */
        return (FALSE);                             /*   then cant subdivide */
    else
	{
        SUBV(separ, Pos(np), Pos(bp));              /* compute displacement  */
        DOTVP(distsq, separ, separ);                /* and find dist squared */

        return (Sizesq(np) > distsq * tolsq);       /* cell size too large ? */
	}
    }

/* endof: differentiatetree.c */
