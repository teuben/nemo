/* transformtree.c - cells, check_potsize, compute_multipoles, intcoord,
                     mk_cells, mk_empty_tree, mk_leaf, mk_tree, ncell,
                     new_cell, nmaxcell, potcorner, potsize, 
                     replace_body_by_cell, replace_root_by_cell,
                     report_monopoles, report_quadrupoles, repot_tree, 
                     stem_of_new_leaf, subindex */

/*
 *  transformtree.c: routines to create body-tree.
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 *
 *        adapted from hackcode versions by Josh Barnes, 1985-7,
 *        with the quadrupole extension by Lars Hernquist, 1986.
 */

#include "newton0.h"

local real  potsize;
local real  potcorner[NDIM];
local int  ncell;
local int  nmaxcell;
local cellptr cells = NULL;	/* cells are allocated from here */

static cellptr mk_empty_tree(systptr sys);
static void grow_new_leaves(systptr sys);
static void mk_leaf(bodyptr bp, systptr sys);
static void check_potsize(bodyptr bp, systptr sys);
static void repot_tree(bodyptr bp, systptr sys);
static nodeptr replace_root_by_cell(nodeptr old_root, real old_root_mid_pos[3 ]);
static nodeptr *stem_of_new_leaf(bodyptr bp, nodeptr *rootptr);
static nodeptr replace_body_by_cell(nodeptr np, int level);
static void compute_multipoles(register nodeptr upper);
static void report_monopoles(nodeptr upper, nodeptr lower);
static void report_quadrupoles(nodeptr upper, nodeptr lower);
static bool intcoord(int xp[3 ], real rp[3 ]);
static int subindex(int x[3 ], int level);
static cellptr mk_cells(int npart);
static cellptr new_cell(void);

/*-----------------------------------------------------------------------------
 *  mk_tree  --   initialize tree structure for hack force calculation.
 *-----------------------------------------------------------------------------
 */
void  mk_tree(sys)
systptr  sys;
    {
    if (cells == NULL)				/* first time through?   */
	cells = mk_empty_tree(sys);        	/*   allocate tree space */

    grow_new_leaves(sys);                          /* scatter the bodies */ 

    compute_multipoles(Root(sys));                 /* gather the poles   */
    }

/*-----------------------------------------------------------------------------
 *  mk_empty_tree  --   allocate memory for a tree structure
 *-----------------------------------------------------------------------------
 */
local cellptr  mk_empty_tree(sys)
systptr  sys;
    {
    real  nbody;               /* number of bodies in the system */
    real  cell_particle_ratio; /* maximum ratio of cell and particle numbers */
    cellptr  new_cells;

    nbody = Nbody(sys);
    cell_particle_ratio = 0.75*((double)(nbody+256))/((double)(nbody+32));
    Nmaxcell(sys) = cell_particle_ratio * nbody;   /* typ. need: 0.5 * nbody */

    new_cells = mk_cells(Nmaxcell(sys));	   /*   allocate cell space  */

    return (new_cells);
    }

/*-----------------------------------------------------------------------------
 *  grow_new_leaves  --  install all bodies of a system "sys" in the tree
 *-----------------------------------------------------------------------------
 */
local void  grow_new_leaves(sys)
systptr  sys;
    {
    bodyptr  bp;

    Root(sys) = NULL;				  /* deallocate current tree */
    potsize = sqrt(Potsizesq(sys));
    SETV(potcorner, Potmincorner(sys)); 
    nmaxcell = Nmaxcell(sys);
    ncell = 0;

    for (bp = Bodies(sys); bp - Bodies(sys) < Nbody(sys); bp++)
	mk_leaf(bp, sys);			  /*   load into body-tree */

    Ncell(sys) = ncell;
    Potsizesq(sys) = potsize * potsize;
    SETV(Potmincorner(sys), potcorner); 
    }

/*-----------------------------------------------------------------------------
 *  mk_leaf  --  insert a body into the tree.
 *-----------------------------------------------------------------------------
 */
local void  mk_leaf(bp, sys)
bodyptr bp;                      /* body to be loaded */
systptr  sys;
    {
    check_potsize(bp, sys);		            /* expand to include bp */
    *(stem_of_new_leaf(bp, &Root(sys))) = (nodeptr)bp;
                                         	    /* and insert into tree */
    }

/*-----------------------------------------------------------------------------
 *  check_potsize  --  check the size of the largest cell in the tree
 *-----------------------------------------------------------------------------
 */
local void  check_potsize(bp, sys)
bodyptr bp;                      /* body to be loaded */
systptr  sys;
    {
    int  xtmp[NDIM];

    if (! intcoord(xtmp, Pos(bp)))		/* expand box (rarely) */
	{
        repot_tree(bp, sys);
	check_potsize(bp, sys);
	}
    }

/*-----------------------------------------------------------------------------
 *  repot_tree  --   enlarge cubical "box", salvaging existing tree structure.
 *-----------------------------------------------------------------------------
 */
local void  repot_tree(bp, sys)
bodyptr  bp;                      /* body to be loaded */
systptr  sys;
    {
    int k;
    real rmid[NDIM];

    ADDVS(rmid, potcorner, 0.5 * potsize);        /*       find box midpoint */
    for (k = 0; k < NDIM; k++)                    /*    loop over dimensions */
        if (Pos(bp)[k] < rmid[k])                 /*      is bp left of mid? */
            potcorner[k] -= potsize;              /*          extend to left */
    potsize *= 2.0;                               /*    double length of box */
    Potsizesq(sys) *= 4.0;                        /*   remember it next time */
    if (Root(sys) != NULL)                        /* existing tree to repot? */
        Root(sys) = replace_root_by_cell(Root(sys), rmid);
    }

/*-----------------------------------------------------------------------------
 *  replace_root_by_cell  --  hang the old root from the new root cell and
 *                            return a node pointer to the new root cell
 *-----------------------------------------------------------------------------
 */
local nodeptr  replace_root_by_cell(old_root, old_root_mid_pos)
nodeptr  old_root;
real  old_root_mid_pos[NDIM];
    {
    int  intpos[NDIM];		  /* integerized coord. vector [0,INTBITMAX) */
    cellptr  new_root;

    new_root = new_cell();		               /* alloc a new cell   */
    Sizesq(new_root) = potsize * potsize;
    intcoord(intpos, old_root_mid_pos);                /* get integer coords */
    Subptr(new_root)[subindex(intpos, INTBITMAX>>1)] = old_root;
                                               /* repot old root in new cell */
    return((nodeptr) new_root);
    }

/*-----------------------------------------------------------------------------
 *  stem_of_new_leaf  --   descend tree and find the address for a new leaf
 *-----------------------------------------------------------------------------
 */
local nodeptr *stem_of_new_leaf(bp, rootptr)
bodyptr  bp;			/* body to load into tree */
nodeptr *rootptr;
    {
    int  level, xp[NDIM], xq[NDIM];
    nodeptr *npp;                            /* npp for node-pointer-pointer */
 
    intcoord(xp, Pos(bp));       		/* form integer coords */

    npp = rootptr;		         	/* start with tree root */
    for(level = INTBITMAX >> 1; level > 0; level >>=1)
        {                                       /* start with top bit */
        if (*npp == NULL)
            return(npp);	                /* found place to store bp */
	if (PartType(*npp) == BODY)		/*   reached a "leaf"? */
            *npp = replace_body_by_cell(*npp, level);
	npp = &Subptr(*npp)[subindex(xp, level)]; /* move down 1 level */
        }
    error("stem_of_new_leaf: ran out of bits ( level = 0 ) ...\n");
    }

/*-----------------------------------------------------------------------------
 *  replace_body_by_cell  --   descend one step and find address for a new cell
 *-----------------------------------------------------------------------------
 */

local nodeptr  replace_body_by_cell(np, level)
nodeptr  np;
int  level;
    {
     int  intpos[NDIM];		  /* integerized coord. vector [0,INTBITMAX) */
    real  cellsize;
    cellptr  cp;
 
    cp = new_cell();		                       /* alloc a new cell   */
    cellsize = potsize / (double)((INTBITMAX >> 1) / level);
    Sizesq(cp) = cellsize * cellsize;

    intcoord(intpos, Pos(np));	                       /* get integer coords */
    Subptr(cp)[subindex(intpos, level)] = np;          /* put body in cell   */

    return((nodeptr) cp);
    }

/*-----------------------------------------------------------------------------
 *  compute_multipoles  --  descend tree to compute multipoles moments of cells
 *-----------------------------------------------------------------------------
 */
local void  compute_multipoles(upper)
register nodeptr upper;                         /* pointer into body-tree */
    {
    register nodeptr *lowerptr;

    if (PartType(upper) == CELL)                     /* is this a cell?      */
        {
        Mass(upper) = 0.0;                             /* clear total mass   */
        CLRV(Pos(upper));			       /* and center of mass */

        for (lowerptr = Subptr(upper); lowerptr-Subptr(upper)<NSUB; lowerptr++)
            if (*lowerptr != NULL)                 
                report_monopoles(upper, *lowerptr);    /* loop over subcells */

        INCDIVVS(Pos(upper), Mass(upper));           /* rescale cms position */

#ifdef QUADPOLE
	CLRM(Quad(upper));			       /* clear quad. moment */
        for (lowerptr = Subptr(upper); lowerptr-Subptr(upper)<NSUB; lowerptr++)
            if (*lowerptr != NULL)                   
                report_quadrupoles(upper, *lowerptr);  /* loop over subcells */
#endif
        }
    }


/*-----------------------------------------------------------------------------
 *  report_monopoles  --  hand up the underlying mass and center-of-mass
 *                         information to the cell above.
 *-----------------------------------------------------------------------------
 */
local void  report_monopoles(upper, lower)
nodeptr  upper, lower;
    {
    real  leverage[NDIM];

    compute_multipoles(lower);
    Mass(upper) += Mass(lower);
    MULVS(leverage, Pos(lower), Mass(lower));
    INCADDV(Pos(upper), leverage);    /* Pos() still has to be normalized !! */
    }

#ifdef QUADPOLE

/*-----------------------------------------------------------------------------
 *  report_quadrupoles  --  hand up the underlying quadrupole information
 *                          to the cell above.
 *-----------------------------------------------------------------------------
 */
local void  report_quadrupoles(upper, lower)
nodeptr  upper, lower;
    {
    real  drsq;
    real  dr[NDIM];
    real  drdr[NDIM][NDIM], Idrsq[NDIM][NDIM], tmpm[NDIM][NDIM];

    SUBV(dr, Pos(lower), Pos(upper));    	/*       displacement vect. */
    OUTVP(drdr, dr, dr);                        /*       outer prod. of dr  */
    DOTVP(drsq, dr, dr);	        	/*       dot prod. dr * dr  */
    SETMI(Idrsq);		        	/*       init unit matrix   */
    MULMS(Idrsq, Idrsq, drsq);            	/*       scale by dr * dr   */
    MULMS(tmpm, drdr, 3.0);        		/*       scale drdr by 3 */
    SUBM(tmpm, tmpm, Idrsq);            	/*       form quad. moment  */
    MULMS(tmpm, tmpm, Mass(lower));     	/*       of cm of subnode,  */
    if (PartType(lower) == CELL)        	/*       if subnode is cell */
        ADDM(tmpm, tmpm, Quad(lower));  	/*         use its moment */
    ADDM(Quad(upper), Quad(upper), tmpm);	/*       add to q.m. of cell */
    }
    
#endif

/*-----------------------------------------------------------------------------
 *  intcoord  --  compute integerized coordinates.
 *                returns: TRUE unless rp was out of bounds.
 *-----------------------------------------------------------------------------
 */
local bool  intcoord(xp, rp)
int  xp[NDIM];			  /* integerized coord. vector [0,INTBITMAX) */
real  rp[NDIM];			  /* real (floating point) coordinate vector */
    {
    register int k;
    bool inb;
    double xsc;

    inb = TRUE;					     /* use to check bounds  */
    for (k = 0; k < NDIM; k++)     		     /* loop over dimensions */
        {
        xsc = (rp[k] - potcorner[k]) / potsize;      /* scale to [0,1)       */
        if (0.0 <= xsc && xsc < 1.0)                 /* within this range?   */
            xp[k] = floor(INTBITMAX * xsc);          /*   then integerize    */
        else                                         /* out of range         */
            inb = FALSE;                             /*   then remember it   */
        }
    return (inb);
    }

/*-----------------------------------------------------------------------------
 *  subindex  --  determine which subcell to select.
 *-----------------------------------------------------------------------------
 */
local int subindex(x, level)
int x[NDIM];		              /* integerized coordinates of particle */
int level;			      /* current level of tree               */
    {
    register int i, k;

    i = 0;                                         /* sum index in i         */
    for (k = 0; k < NDIM; k++)                     /* check each dimension   */
        if (x[k] & level)                          /*   if beyond midpoint   */
            i += NSUB >> (k + 1);                  /*     skip over subcells */

    return (i);
    }

/*-----------------------------------------------------------------------------
 *  mk_cells  --  allocates memory for cells used in the tree version of
 *                newton0.
 *                accepts: npart: the number of cells which can be accomodated.
 *                returns: new_cells: points to newly created array of cells.
 *-----------------------------------------------------------------------------
 */
local cellptr  mk_cells(npart)
int  npart;
    {
    cellptr  new_cells;

    new_cells = (cellptr) malloc((unsigned)npart * sizeof(cell));
    if (new_cells == NULL)
	error("mk_cells: not enough memory left for %d cells\n", npart);
    return(new_cells);
    }

/*-----------------------------------------------------------------------------
 *  new_cell  --   allocation routine for cells.
 *-----------------------------------------------------------------------------
 */
local  cellptr new_cell()
    {
    register cellptr c;
    register int i;

    if (ncell >= nmaxcell)
	error("new_cell: need more than %d cells; increase Nmaxcell\n",
	      nmaxcell);
    c = cells + ncell;
    ncell++;
    PartType(c) = CELL;
    for (i = 0; i < NSUB; i++)
	Subptr(c)[i] = NULL;
    return (c);
    }

/* endof: transformtree.c */
