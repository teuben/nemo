/*
 *  1D gridding utilities : linear, lookup table and external lookup function
 *        30-oct-93	Created				peter teuben
 *	  12-apr-95	no more ARGS in prototypes
 *        15-mar-05     proper real_proc for c++
 */
 
#define GRID_LINEAR  1
#define GRID_ARRAY   2
#define GRID_PROC    3

typedef struct grid {
    int mode;		/* one of the GRID_xxx modes */
    int n;              /* number of grid cells */

    bool up;		/* TRUE if gmin < gmax */
    real gmin, gmax;    /* requested min and max for gridder */
    real dg;            /* gridsize of one grid cell in linear mode */
    real *g;            /* monotonically increasing array of grid boundaries */

    real_proc f;        /* external function 'real f(real)' in case non-linear functions */
} Grid, *GridPtr;

extern void inil_grid (Grid *, int, real, real);
extern void inia_grid (Grid *, int, real *);
extern void inip_grid (Grid *, int, real, real, real_proc);
extern int  index_grid (Grid *, real);
extern real value_grid (Grid *, real);
