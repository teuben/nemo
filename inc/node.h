/* node.h - CELL, INTBITMAX, NODE, Quad, Sizesq, Subptr, cell, cellptr, node, 
            nodeptr */

/*
 *  node.h: for nodes and cells in a tree containing bodies as their leaves
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
/*-----------------------------------------------------------------------------
 *  node, nodeptr  --  nodes are a generic form of particle type which includes
 *                     the classes of bodies and cells. Nodes contain the
 *                     minimum dynamic information for a particle, and are 
 *                     therefore equivalent to mbodies (see  body.h ).
 *                         Nodes are useful when performing operations which
 *                     apply equally to bodies and to cells; in these cases
 *                     nodes allow generic operations on particles.
 *                     An example is the definition of subptr in the defining
 *                     expression for a cell structure as being a nodeptr:
 *                     subnodes of cells can equally well be bodies as cells.
 *                         Another use of nodes is to express to the
 *                     human reader of a piece of code the intent of treating
 *                     generic particles, bodies as well as cells.
 *-----------------------------------------------------------------------------
 */
typedef struct
    {
    short  type;                 /* type of a node                      */
    real  phasespace[2][NDIM];   /* position and velocity vectors       */
    real  mass;                  /* mass of the node                    */
    } node, *nodeptr;

/*-----------------------------------------------------------------------------
 *  cell, cellptr  --  cells are a form of pseudo particles which summarize
 *                     the dynamical information of the cluster of particles
 *                     which they (hierarchically) point to through their
 *                     node pointers (i.e. their branches). In the present
 *                     implementation cells contain at the very least their
 *                     type, mass and position, and their set of node pointers.
 *                     In addition, they may contain higher order multipole
 *                     moments, starting with quadrupole moments.
 *                     Cells are situated at the center of mass of their
 *                     subservient particles, and hence the dipole moments
 *                     of their underlying mass distributions vanish
 *                     identically.
 *-----------------------------------------------------------------------------
 */
#define NSUB (1 << NDIM)   /* number of subcells per cell in NDIM dimensions */

typedef struct
    {
    short  type;                      /* type of a cell                      */
    real  phasespace[2][NDIM];        /* position and velocity vectors       */
    real  mass;                       /* mass of the cell                    */
    real  sizesq;                     /* square of cell size (between        */
#ifdef QUADPOLE                                      /* neighboring corners) */
    real  quad[NDIM][NDIM];           /* Cartesian quadrupole components     */
#endif
    nodeptr  subptr[NSUB];            /* pointers to descendants of the cell */
    } cell, *cellptr;

/*-----------------------------------------------------------------------------
 *  macros to extract individual (sub)components from a  cell or node:
 *  
 *  the following macros (defined in  body.h  for bodies and xbodies) can be 
 *  used without modification for cells and nodes as well:
 *
 *      #define  PartType(ptr)  ((ptr)->type)               [* type: short   *]
 *      #define  Pos(ptr)       ((ptr)->phasespace[0])      [* type: realptr *]
 *      #define  Vel(ptr)       ((ptr)->phasespace[1])      [* type: realptr *]
 *      #define  Mass(ptr)      ((ptr)->mass)               [* type: real    *]
 *      #define  Phase(ptr)     ((real *)(ptr)->phasespace) [* type: realptr *]
 *
 *  note: PartType addresses the type of particle, used here as a generic
 *        name for either a body or a cell.
 *
 *  no new macros are needed to describe nodes; cells, however, do need extra
 *  macros, as given below.
 *  NOTE:  ptr should be a pointer of type "cellptr":
 *-----------------------------------------------------------------------------
 */                        
#define  Sizesq(ptr)     (((cellptr) (ptr))->sizesq)       /* type: real     */
#define  Subptr(ptr)     (((cellptr) (ptr))->subptr)       /* type: *nodeptr */
#ifdef QUADPOLE
#  define  Quad(ptr)     (((cellptr) (ptr))->quad)         /* type: *realptr */
#endif
 
/*-----------------------------------------------------------------------------
 *  encoding for different types of particles:
 *  note: a node can carry a definite type, either a BODY or a CELL;
 *        however, a separate type NODE indicates that no choice has yet
 *        been made as to its identity. Therefore, a particle of PartType
 *        NODE necessarily is a node; but not every node has to carry the
 *        type NODE.
 *-----------------------------------------------------------------------------
 */
#define   CELL    020
#define   NODE    030

/*-----------------------------------------------------------------------------
 *  the following is a macro for the highest bit allowed in an integerized
 *  representation of position coordinates, which are used to hang particles
 *  in a tree (bit pattern if integer representation of position forms the
 *  address in an Eulerian tree).
 *-----------------------------------------------------------------------------
 */
#define INTBITMAX (1 << (8 * sizeof(int) - 2))

/* endof: node.h */
