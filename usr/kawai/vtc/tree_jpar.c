#pragma global noalias
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "vtc.h"
#include "vtclocal.h"
#include "mdgp2defs.h"

typedef int Node;

#define NBODYMAX (2300000)
#define NCELLMAX (2300000)

/*
#define NBODYMAX (1100000)
#define NCELLMAX (1100000)
*/

/* body attrs */
typedef int BHindex;
static BHindex Body[NBODYMAX];

/* cell attrs */
#define NOCELL (-1)
typedef int Cell;
static double Cpos[3][NCELLMAX];
static double Ccmpos[3][NCELLMAX];
static double Ccmmass[NCELLMAX];
static double Csize[NCELLMAX];
static Cell Cparent[NCELLMAX];
static Cell Cchild[8][NCELLMAX];
static int Cchildnbody[8][NCELLMAX];
static int Cisleaf[NCELLMAX];
static int Cnbody[NCELLMAX];
static int Cndescendant[NCELLMAX];
static BHindex *Cbody[NCELLMAX];
static double Ccrit[NCELLMAX];
static Cell celltmp[NCELLMAX]; /* cells at the level in question */
static Cell celltmp_new[NCELLMAX]; /* cells at the next level */

static int cellp = -1; /* # of cells */

/* local funcs */
static int make_cell(void);
static void reset_cellbuf(int n, int npp);
static int sanity_check(Cell c0, int level, Nbodyinfo *nb);
#define NLISTMAX (80000)
static void
calculate_force_from_lists(int nlist, BHindex **nodelist,
			   int ni[NBOARDMAX], double (**xi)[3],
			   int nj[NBOARDMAX], double (**xj)[3], double **mj,
			   double (**a)[3], double **p,
			   Forceinfo *fi, Nbodyinfo *nb);

/* 'typical' params for tree algorithm */
void
vtc_get_default_tree_params(Forceinfo *tr)
{
    tr->node_div_crit = 8;
    tr->ncrit = 8192;
    tr->theta = 0.75;
    tr->eps = 0.02;
    tr->p = 1;
    tr->npp = 1;
    tr->full_dof = TRUE;
    tr->calculator = HOST;
}

static int
create_root_cell(Forceinfo *fi, Nbodyinfo *nb)
{
    int root;
    int i, k;
    double cs, val, posmax = 0.0;
    static double csold = 0.0;

    reset_cellbuf(nb->n, fi->npp);
    root = make_cell();
    for (k = 0; k < 3; k++) {
	Cpos[k][root] = 0.0;
    }
    cs = 1.0;
    posmax = 0.0;
    for (i = 0; i < nb->n; i++) {
	for (k = 0; k < 3; k++) {
	    val = fabs(nb->x[k][i]);
	    if (posmax < val) {
		posmax = val;
	    }
	}
    }
    while (posmax*2 > cs) {
	cs *= 2.0;
    }
    while (posmax*2 < cs*0.5) {
	cs *= 0.5;
    }
    if (cs != csold) {
	csold = cs;
	vtc_set_scale(cs);
    }
    Csize[root] = cs;
    Ccrit[root] = cs*cs/(fi->theta)/(fi->theta);
    return (root);
}

static void
assign_bodies(int n, /* # of bodies */
	      BHindex *body, /* bodies to be assiged */
	      Cell *bpid, /* cell the body currently assigned to */
	      int *bsubid, /* subid of the child of the cell to which
			      the body should be assigned */
	      double *bpos[3]) /* the body position */
{
#if 1 /* vector */
    int i, sid;
    double *bpos0 = bpos[0];
    double *cpos0 = Cpos[0];
    double *bpos1 = bpos[1];
    double *cpos1 = Cpos[1];
    double *bpos2 = bpos[2];
    double *cpos2 = Cpos[2];

    for (i = 0; i < n; i++) {
	sid = 0;
	if (bpos0[body[i]] > cpos0[bpid[i]]) {
	    sid += 1;
	}
	if (bpos1[body[i]] > cpos1[bpid[i]]) {
	    sid += 2;
	}
	if (bpos2[body[i]] > cpos2[bpid[i]]) {
	    sid += 4;
	}
	bsubid[i] = sid;
    }

#else /* scalar */
    int i, sid;

#pragma loop noalias
    for (i = 0; i < n; i++) {
	sid = 0;
	if (bpos[0][body[i]] > Cpos[0][bpid[i]]) {
	    sid += 1;
	}
	if (bpos[1][body[i]] > Cpos[1][bpid[i]]) {
	    sid += 2;
	}
	if (bpos[2][body[i]] > Cpos[2][bpid[i]]) {
	    sid += 4;
	}
	bsubid[i] = sid;
    }
#endif
}

static void
count_particles(int n, /* # of bodies */
		Cell *bpid, /* cell the body currently assigned to */
		int *bsubid) /* subid of the child of the cell to which
				the body should be assigned */
{
    int i, s;

    for (i = 0; i < n; i++) {
	Cchildnbody[bsubid[i]][bpid[i]]++;
    }
}

static int /* # of non-leaf cells newly opened */
open_cells(int ndiv, /* min # of particle to process cell division */
	   int ncell0, Cell (*cell0)[NCELLMAX], /* cells in question */
	   Cell (*cell1)[NCELLMAX], /* cells newly opened */
	   BHindex **bodypp, /* pointer to the max index of Body[] in use */
	   double theta2) 
{
    int i, s, k, cn1;
    int ncell1 = 0;
    Cell c0, c1;
    double newsize = 0.5*Csize[(*cell0)[0]];
    double offset[3];

    for(s = 0; s < 8; s++) {
	for (k = 0; k < 3; k++) {
	    if ((s>>k)%2) {
		offset[k] = 0.5*newsize;
	    }
	    else {
		offset[k] = -0.5*newsize;
	    }
	}
	for (i = 0; i < ncell0; i++) {
	    c0 = (*cell0)[i];
	    /* child s of cell c0 is in question */
	    cn1 = Cchildnbody[s][c0];
	    if (0 == cn1) {
		continue; /* empty cell. do not process anymore */
	    }
	    /* open this cell */
	    Cchild[s][c0] = c1 = make_cell();
	    Cparent[c1] = c0;
	    Csize[c1] = newsize;
	    Ccrit[c1] = newsize*newsize/theta2;
	    for (k = 0; k < 3; k++) {
		Cpos[k][c1] = Cpos[k][c0] + offset[k];
	    }
	    if (cn1 > ndiv) {
		(*cell1)[ncell1] = c1;
		ncell1++;
	    }
	    else {
		Cisleaf[c1] = TRUE;
		Cbody[c1] = *bodypp;
		(*bodypp) += cn1;
	    }
	}
    }
    return (ncell1);
}

#if 1 /* partly vectorized */
int /* # of bodies to be processed at the next level */
add_particles(int n0, /* # of bodies to be processed */
	      BHindex *body0, /* bodies to be processed */
	      BHindex *body1, /* bodies left to be processed at the next level */
	      Cell *bpid, /* cell the body currently assigned to */
	      int *bsubid) /* subid of the child of the cell to which
			      the body should be added */
{
    int i, n1 = 0;
    Cell c0, c1;
    static int notyet[NBODYMAX];

    for (i = 0; i < n0; i++) {
	notyet[i] = TRUE;
    }

    for (i = 0; i < n0; i++) {
	c0 = bpid[i];
	c1 = Cchild[bsubid[i]][c0];
	if (Cisleaf[c1]) {
	    Cbody[c1][Cnbody[c1]] = body0[i];
	    Cnbody[c1]++;
	    notyet[i] = FALSE;
	}
    }
#pragma loop novrec
    for (i = 0; i < n0; i++) {
	if (notyet[i]) {
	    c0 = bpid[i];
	    c1 = Cchild[bsubid[i]][c0];
	    body1[n1] = body0[i]; /* prepare bodies for the next level */
	    bpid[n1] = c1;
	    n1++;
	}
    }
    return(n1);
}

#else /* scalar */

static int /* # of bodies to be processed at the next level */
add_particles(int n0, /* # of bodies to be processed */
	      BHindex *body0, /* bodies to be processed */
	      BHindex *body1, /* bodies left to be processed at the next level */
	      Cell *bpid, /* cell the body currently assigned to */
	      int *bsubid) /* subid of the child of the cell to which
			      the body should be added */
{
    int i, n1 = 0;
    Cell c0, c1;

    for (i = 0; i < n0; i++) {
	c0 = bpid[i];
	c1 = Cchild[bsubid[i]][c0];
	if (Cisleaf[c1]) {
	    Cbody[c1][Cnbody[c1]] = body0[i];
	    Cnbody[c1]++;
	}
	else {
	    body1[n1] = body0[i]; /* prepare bodies for the next level */
	    bpid[n1] = c1;
	    n1++;
	}
    }
    return(n1);
}

#endif

static Cell /* root cell */
create_tree(Forceinfo *fi, Nbodyinfo *nb)
{
    int root;
    int i;
    int ncelltmp; /* # of cells at the level in question */
    int ncelltmp_new; /* # of cells processed at the next level */
    int nbody; /* # of particles at the level in question */
    int nbody_new; /* # of particles to be processed at the next level */
    double theta2;
    static BHindex body[NBODYMAX]; /* particles at the level in question */
    static BHindex body_new[NBODYMAX]; /* particles to be processed at the next level */
    static Cell bpid[NCELLMAX]; /* particle's parent cell */
    static int bsubid[NBODYMAX]; /* bpid's child cell which the particle belongs to */
    static BHindex *bodyp;

    root = create_root_cell(fi, nb);

    if (NBODYMAX < nb->n) {
	fprintf(stderr, "too large n: %d NMAX: %d\n",
		nb->n, NBODYMAX);
	exit(1);
    }
    ncelltmp = 1;
    celltmp[0] = root;
    nbody = nb->n;
    for (i = 0; i < nbody; i++) {
	bpid[i] = root;
	body[i] = i;
    }
    bodyp = Body;
    theta2 = (fi->theta)*(fi->theta);
    while (ncelltmp > 0) {
	assign_bodies(nbody, body, bpid, bsubid, nb->x);
	count_particles(nbody, bpid, bsubid);

	ncelltmp_new = open_cells(fi->node_div_crit, ncelltmp, &celltmp,
				  &celltmp_new, &bodyp, theta2);
	nbody_new = add_particles(nbody, body, body_new, bpid, bsubid);
	ncelltmp = ncelltmp_new;
	for (i = 0; i < ncelltmp; i++) {
	    celltmp[i] = celltmp_new[i];
	}
	nbody = nbody_new;
	for (i = 0; i < nbody; i++) {
	    body[i] = body_new[i];
	}
    }
    return (root);
}


#if 1

#define LEVELMAX (32)
static void
m2m_1storder(Forceinfo *fi, Nbodyinfo *nb)
{
    int c, k, s, lv;
    double m1;
    int level;
    Cell first[LEVELMAX];
    int *cb;
    double *cmp0, *cmp1, *cmp2;
    double *nbm, *nbx0, *nbx1, *nbx2;

    for (c = 0; c <= cellp; c++) {
	Ccmmass[c] = 0.0;
	Ccmpos[0][c] = 0.0;
	Ccmpos[1][c] = 0.0;
	Ccmpos[2][c] = 0.0;
    }
    first[0] = 0;
    level = 0;
    Cndescendant[0] = 0;
    for (c = 1; c <= cellp; c++) {
	if (Csize[c] < Csize[c-1]*0.8) { /* go down the tree */
	    level++;
	    first[level] = c;
	}
	Cndescendant[c] = 0;
    }
    first[level+1] = c;

    cmp0 = Ccmpos[0];
    cmp1 = Ccmpos[1];
    cmp2 = Ccmpos[2];
    nbm = nb->m;
    nbx0 = nb->x[0];
    nbx1 = nb->x[1];
    nbx2 = nb->x[2];
    for (lv = level; lv >= 0; lv--) {
	/* leaf cells */
#pragma loop novrec
	for (c = first[lv]; c < first[lv+1]; c++) {
	    int i;
	    cb = Cbody[c];
	    if (!Cisleaf[c]) {
		continue;
	    }
	    for (i = 0; i < Cnbody[c]; i++) {
		BHindex b;
		b = cb[i];
		m1 = nbm[b];
		cmp0[c] += nbx0[b]*m1;
		cmp1[c] += nbx1[b]*m1;
		cmp2[c] += nbx2[b]*m1;
		Ccmmass[c] += m1;
	    }
	    Cndescendant[c] = Cnbody[c];
	}

	/* non-leaf cells */
	for (s = 0; s < 8; s++) {
	    int *ccs;
	    ccs = Cchild[s];
#pragma loop novrec
	    for (c = first[lv]; c < first[lv+1]; c++) {
		Cell c1;
		c1 = ccs[c];
		if (Cisleaf[c]) {
		    continue;
		}
		if (NOCELL == c1) {
		    continue;
		}
		m1 = Ccmmass[c1];
		cmp0[c] += cmp0[c1]*m1;
		cmp1[c] += cmp1[c1]*m1;
		cmp2[c] += cmp2[c1]*m1;
		Ccmmass[c] += m1;
		Cndescendant[c] += Cndescendant[c1];
	    }
	}

	for (c = first[lv]; c < first[lv+1]; c++) {
	    cmp0[c] /= Ccmmass[c];
	    cmp1[c] /= Ccmmass[c];
	    cmp2[c] /= Ccmmass[c];
	}
    }
}

#else

static void
m2m_1storder(Cell c0, Forceinfo *fi, Nbodyinfo *nb)
{
    int i, k, s;
    Cell c1;
    BHindex b;
    double m1;

    for (k = 0; k < 3; k++) {
	Ccmpos[k][c0] = 0.0;
    }
    Ccmmass[c0] = 0.0;
    Cndescendant[c0] = 0;
    if (Cisleaf[c0]) {
	for (i = 0; i < Cnbody[c0]; i++) {
	    b = Cbody[c0][i];
	    m1 = nb->m[b];
	    for (k = 0; k < 3; k++) {
		Ccmpos[k][c0] += (nb->x[k][b])*m1;
	    }
	    Ccmmass[c0] += m1;
	}
	Cndescendant[c0] = Cnbody[c0];
    }
    else {
	for (s = 0; s < 8; s++) {
	    c1 = Cchild[s][c0];
	    if (NOCELL == c1) {
		continue;
	    }
	    m2m_1storder(c1, fi, nb);
	    m1 = Ccmmass[c1];
	    for (k = 0; k < 3; k++) {
		Ccmpos[k][c0] += Ccmpos[k][c1]*m1;
	    }
	    Ccmmass[c0] += m1;
	    Cndescendant[c0] += Cndescendant[c1];
	}
    }
    for (k = 0; k < 3; k++) {
	Ccmpos[k][c0] /= Ccmmass[c0];
    }
}

#endif

static void
m2m_2ndorder(Cell root, Forceinfo *fi, Nbodyinfo *nb)
{
    fprintf(stderr, "m2m_2ndorder not implemented yet\n");
}

static void
m2m_anyorder(Cell root, Forceinfo *fi, Nbodyinfo *nb)
{
    fprintf(stderr, "m2m_anyorder not implemented yet\n");
}

static void
process_m2m(Cell root, Forceinfo *fi, Nbodyinfo *nb)
{
    if (!fi->full_dof || fi->p > 2) {
	m2m_anyorder(root, fi, nb);
    }
    else if (1 == fi->p) {
	m2m_1storder(fi, nb);
    }
    else if (2 == fi->p) {
	m2m_2ndorder(root, fi, nb);
    }
    else {
	fprintf(stderr, "process_m2m: order %d invalid\n", fi->p);
	exit(1);
    }
}

static int /* length of the list */
get_descendant_list(Cell c0, Forceinfo *fi, Nbodyinfo *nb, BHindex *list)
{
    Cell c1;
    int i, s;
    int len0, len1;

    len0 = 0;
    if (Cisleaf[c0]) {
	len0 = Cnbody[c0];
	for (i = 0; i < len0; i++) {
	    list[i] = Cbody[c0][i];
	}
    }
    else {
	for (s = 0; s < 8; s++) {
	    c1 = Cchild[s][c0];
	    if (NOCELL == c1) {
		continue;
	    }
	    len0 += get_descendant_list(c1, fi, nb, list+len0);
	}
    }
    return (len0);
}

static double /* (separation)^2 */
get_separation(double isize, double ipos[3], Cell jcell)
{
    int k;
    double dx, dy, dz, dr2;
    double xmin, jsize;

    jsize = Csize[jcell];
    xmin = (isize)*0.5; /* i & j are guaranteed not to be overlaped */

/*  more strict criterion
    xmin = (isize+jsize)*0.5;
*/
    dx = fabs(ipos[0]-Cpos[0][jcell]);
    dy = fabs(ipos[1]-Cpos[1][jcell]);
    dz = fabs(ipos[2]-Cpos[2][jcell]);
    dx -= xmin;
    dy -= xmin;
    dz -= xmin;
    if (dx < 0.0) {
	dx = 0.0;
    }
    if (dy < 0.0) {
	dy = 0.0;
    }
    if (dz < 0.0) {
	dz = 0.0;
    }
    dr2 = dx*dx+dy*dy+dz*dz;
    return (dr2);
}

#if 1
static int /* length of the list excluding i-particles */
get_interaction_list(Cell icell, Cell root, double theta2, Nbodyinfo *nb,
		     Node *list, int *iscell, int dummy)
{
    static int qualified[NCELLMAX];
    static double dr2[NCELLMAX];
    int i, j, k, s;
    Cell jcell0, jcell1;
    int nlist;
    int ncelltmp; /* # of cells at the level in question */
    int ncelltmp_new; /* # of cells processed at the next level */
    double isize, ipos[3];
    double xmin;

    isize = Csize[icell];
    xmin = (isize)*0.5;
    for (k = 0; k < 3; k++) {
	ipos[k] = Cpos[k][icell];
    }
    nlist = 0;
    ncelltmp = ncelltmp_new = 1;
    celltmp[0] = root;

    while (ncelltmp > 0) {
	double *cp0, *cp1, *cp2;
	ncelltmp = ncelltmp_new;
	ncelltmp_new = 0;

	cp0 = Cpos[0];
	cp1 = Cpos[1];
	cp2 = Cpos[2];

	for (i = 0; i < ncelltmp; i++) {
	    double dx, r;

	    r = 0.0;
	    dx = ipos[0]-cp0[celltmp[i]];
	    dx = dx > 0.0 ? dx:-dx;
	    dx -= xmin;
	    dx = dx > 0.0 ? dx:0.0;
	    r += dx*dx;

	    dx = ipos[1]-cp1[celltmp[i]];
	    dx = dx > 0.0 ? dx:-dx;
	    dx -= xmin;
	    dx = dx > 0.0 ? dx:0.0;
	    r += dx*dx;

	    dx = ipos[2]-cp2[celltmp[i]];
	    dx = dx > 0.0 ? dx:-dx;
	    dx -= xmin;
	    dx = dx > 0.0 ? dx:0.0;
	    r += dx*dx;

	    dr2[i] = r;
	}
	for (i = 0; i < ncelltmp; i++) {
	    qualified[i] = (dr2[i] * theta2 > (Csize[celltmp[i]])*(Csize[celltmp[i]])) ? TRUE:FALSE;
	    if (icell == celltmp[i]) {
		qualified[i] = -1; /* neither TRUE nor FALSE */
	    }
	}

	/* apply multipole expansion */
	for (i = 0; i < ncelltmp; i++) {
	    if (qualified[i] == TRUE) {
		list[nlist] = celltmp[i];
		iscell[nlist] = TRUE;
		nlist++;
		qualified[i] = -1;
	    }
	}
	for (i = 0; i < ncelltmp; i++) {
	    if ((qualified[i] == FALSE) && Cisleaf[celltmp[i]]) {
	      qualified[i] = TRUE;
	    }
	}
	/* leaf cell. handle particle force directly */
	for (i = 0; i < ncelltmp; i++) {
	    if (qualified[i] == TRUE) {
		BHindex *cb = Cbody[celltmp[i]];
		for (j = 0; j < Cnbody[celltmp[i]]; j++, nlist++) {
		    list[nlist] = cb[j];
		    iscell[nlist] = FALSE;
		}
		Cfprintf(stderr, "leaf: %d\n", Cnbody[celltmp[i]]);
		qualified[i] = -1;
	    }
	}
	for (s = 0; s < 8; s++) {
	    for (i = 0; i < ncelltmp; i++) {
		if (qualified[i] == FALSE) { /* descend the tree */
		    jcell1 = Cchild[s][celltmp[i]];
		    if (NOCELL != jcell1) {
			celltmp_new[ncelltmp_new] = jcell1;
			ncelltmp_new++;
		    }
		}
	    }
	}
	for (i = 0; i < ncelltmp_new; i++) {
	    celltmp[i] = celltmp_new[i];
	}
    } /* while */
    return (nlist);
}
#else
static int /* length of the list excluding i-particles */
get_interaction_list(Cell icell, Cell jcell0, double theta2, Nbodyinfo *nb,
		     BHindex *list, int *iscell, int level)
{
    int i, s;
    int len0 = 0;
    Cell jcell1;
    static double isize, ipos[3];

    if (level == 0) {
	int k;
	isize = Csize[icell];
	for (k = 0; k < 3; k++) {
	    ipos[k] = Cpos[k][icell];
	}
    }
    if (icell == jcell0) { /* self interaction is already taken into account */
	len0 = 0;
    }
#if 1
    else if (get_separation(isize, ipos, jcell0) * theta2 > (Csize[jcell0])*(Csize[jcell0])) { 
#else
    else if (get_separation(isize, ipos, jcell0) > Ccrit[jcell0]) {/* apply multipole expansion */
#endif
	list[0] = jcell0;
	iscell[0] = TRUE;
	len0 = 1;
    }
    else { /* not well separated */
	if (Cisleaf[jcell0]) { /* handle particle force directly */
	    BHindex *cb = Cbody[jcell0];
	    len0 = Cnbody[jcell0];
	    for (i = 0; i < len0; i++) {
		list[i] = cb[i];
		iscell[i] = FALSE;
	    }
	}
	else { /* descend the tree */
	    for (s = 0; s < 8; s++) {
		jcell1 = Cchild[s][jcell0];
		if (NOCELL == jcell1) {
		    continue;
		}
		len0 += get_interaction_list(icell, jcell1, theta2, nb,
					     list+len0, iscell+len0, level++);
	    }
	}
    }
    return (len0);
}

#endif

static void
traverse_tree(Cell c0, Cell root, Forceinfo *fi, Nbodyinfo *nb)
{
    Cell c1;
    int node; /* index to a cell or a particle */
    int i, j, k, s, b;
    double *cmp0, *cmp1, *cmp2;
    double *nbm, *nbx0, *nbx1, *nbx2;
    BHindex *nl;
    double *ml, (*xl)[3];
    static int ilen[NBOARDMAX], jlen[NBOARDMAX];
    static int nwalk;
    static double pptot;
    static int iscell[NLISTMAX]; /* node is a cell or a particle */
    static BHindex *nodelist[NBOARDMAX];
    static double *mlist[NBOARDMAX];
    static double (*xlist[NBOARDMAX])[3];
    static double (*alist[NBOARDMAX])[3];
    static double *plist[NBOARDMAX];
    static int nisum;
    static int nlist;
    static int nboard;
    static int firstcall = 1;
    static int verybig;
    static double maxratio;

    if (firstcall) {
	firstcall = 0;
	nboard = vtc_get_grape();
	fprintf(stderr, "nboard: %d\n", nboard);
	if (nboard > NBOARDMAX) {
	    fprintf(stderr, "too large nboard.\n", nboard);
	    exit(1);
	}
	for (b = 0; b < nboard; b++) {
	    nodelist[b] = (BHindex *)malloc(sizeof(BHindex)*NLISTMAX);
	    mlist[b] = (double *)malloc(sizeof(double)*NLISTMAX);
	    xlist[b] = (double (*)[3])malloc(sizeof(double)*3*NLISTMAX);
	    alist[b] = (double (*)[3])malloc(sizeof(double)*3*NLISTMAX);
	    plist[b] = (double *)malloc(sizeof(double)*NLISTMAX);
	}
    }
    if (c0 == root) {
	nlist = 0;
	nisum = 0;
	pptot = 0.0;
	nwalk = 0;
	verybig = 0;
	maxratio = 0.0;
    }
    if (Cndescendant[c0] > fi->ncrit && !Cisleaf[c0]) {
	for (s = 0; s < 8; s++) {
	    c1 = Cchild[s][c0];
	    if (NOCELL == c1) {
		continue;
	    }
	    traverse_tree(c1, root, fi, nb);
	}
    }
    else {
	ilen[nlist] = get_descendant_list(c0, fi, nb, nodelist[nlist]);
	if (ilen[nlist] != Cndescendant[c0]) {
	    fprintf(stderr, "ilen[%d] (%d) differs from Cndescendant[%d] (%d)\n",
		    ilen[nlist], nlist, c0, Cndescendant[c0]);
	    exit(1);
	}
	for (i = 0; i < ilen[nlist]; i++) {
	    iscell[i] = FALSE;
	}
	nisum += ilen[nlist];
	Cprintf("ilen0[%d]: %d\n", nlist, ilen[nlist]);
	jlen[nlist] = ilen[nlist] + get_interaction_list(c0, root, (fi->theta)*(fi->theta), nb,
							 nodelist[nlist]+ilen[nlist], iscell+ilen[nlist], 0);
	if (NLISTMAX < jlen[nlist]) {
	    fprintf(stderr, "list too long. jlen[%d]: %d\n", nlist, jlen[nlist]);
	    exit(1);
	}

	cmp0 = Ccmpos[0];
	cmp1 = Ccmpos[1];
	cmp2 = Ccmpos[2];
	nbm = nb->m;
	nbx0 = nb->x[0];
	nbx1 = nb->x[1];
	nbx2 = nb->x[2];
	nl = nodelist[nlist];
	ml = mlist[nlist];
	xl = xlist[nlist];
	for (j = 0; j < jlen[nlist]; j++) {
	    node = nl[j];
	    if (iscell[j]) { /* node is an index to a cell */
		ml[j] = Ccmmass[node];
		xl[j][0] = cmp0[node];
		xl[j][1] = cmp1[node];
		xl[j][2] = cmp2[node];
	    }
	    else { /* node is an index to a particle */
		ml[j] = nbm[node];
		xl[j][0] = nbx0[node];
		xl[j][1] = nbx1[node];
		xl[j][2] = nbx2[node];
	    }
	}
	nwalk++;
	pptot += ilen[nlist]*jlen[nlist];
	nlist++;
    }
    if (nlist == nboard || (c0 == root && nlist != 0)) {
	int b, id, k;
#if 0
	int big = 0, small = NLISTMAX;

	for (b = 0; b < nlist; b++) {
	    if (big < jlen[b]) {
		big = jlen[b];
	    }
	    if (small > jlen[b]) {
		small = jlen[b];
	    }
	}
	if (big/small > 0.5) {
	    verybig++;
	}
	maxratio += big/small;
#endif
	calculate_force_from_lists(nlist, nodelist, ilen, xlist,
				   jlen, xlist, mlist,
				   alist, plist,
				   fi, nb);
	nlist = 0;
    }
    if (c0 == root) {
	fprintf(stderr, "nisum: %d\n", nisum);
	fi->ninteraction = pptot;
	fi->nwalk = nwalk;
	fprintf(stderr, "verybig: %d  nwalk: %d  verybig/nwalk: %f big/small avg: %f\n",
		verybig, nwalk, (double)verybig/nwalk, maxratio/nwalk);
    }
}

/* calculate force from nlist interaction lists to nlist acceleration boxes */
static void
calculate_force_from_lists(int nlist, BHindex **nodelist,
			   int ni[NBOARDMAX], double (**xi)[3],
			   int nj[NBOARDMAX], double (**xj)[3], double **mj,
			   double (**alist)[3], double **plist,
			   Forceinfo *fi, Nbodyinfo *nb)
{
    int b, i, k;
    int node; /* index to a cell or a particle */
    double epsinv;

    if (fi->eps != 0.0) { 
	epsinv = 1.0/fi->eps;
    }
    else {
	epsinv = 0.0;
    }
    (vtc_force_calculatorMB[fi->calculator])(nlist, ni, xi, nj, xj, mj, fi->eps,
					     alist, plist);
    for (b = 0; b < nlist; b++) { /* for each board */
	if (fi->calculator != GRAPE_POTENTIALONLY) {
	    for (i = 0; i < ni[b]; i++) {
		node = nodelist[b][i];
		for (k = 0; k < 3; k++) {
#if 1
		    nb->a[k][node] = alist[b][i][k];
#else /* returns no force */
		    nb->a[k][node] = 0.0;
#endif
		}
	    }
	}
	if (fi->calculator != GRAPE_FORCEONLY) {
	    for (i = 0; i < ni[b]; i++) {
		node = nodelist[b][i];
		nb->p[node] = plist[b][i] + nb->m[node] * epsinv;

		/* MD2 does not calc self interaction...but on vpp5k
		 * this correction seems to be necessary.
		 * I don't know why...
		 nb->p[node] = plist[b][i];
		 nb->p[node] = plist[b][i] + nb->m[node] * epsinv;
		*/
	    }
	}
    }
#if 0
    {
	static int err = 0;

	for (i = 0; i < nb->n; i++) {
	    for (k = 0; k < 3; k++) {
		if (fabs(nb->a[k][i]) > 10.0 && !err) {
		    fprintf(stderr, "\n\n!!! nb->a[%d][%d]: %f\n\n\n",
			    k, i, nb->a[k][i]);
		    err = 1;
		}
	    }
	}
    }
#endif
}

static void
print_pos(Nbodyinfo *nb)
{
    int i;
    for (i = 0; i <2; i++) {
	fprintf(stderr, "x: %6.5f %6.5f %6.5f a: %6.5f %6.5f %6.5f\n",
		nb->x[0][i], nb->x[1][i], nb->x[2][i],
		nb->a[0][i], nb->a[1][i], nb->a[2][i]);
    }
    fprintf(stderr, "\n");
    exit(1);
}

static Cell
make_cell(void)
{
    int s;

    cellp++;
    if (cellp >= NCELLMAX) {
	fprintf(stderr, "cell buf overflown. NCELLMAX: %d\n",
		NCELLMAX);
	exit(1);
    }
    Cisleaf[cellp] = FALSE;
    Cnbody[cellp] = 0;
    Cbody[cellp] = NULL;
    for (s = 0; s < 8; s++) {
	Cchildnbody[s][cellp] = 0;
	Cchild[s][cellp] = NOCELL;
    }
    return (cellp);
}

static void
reset_cellbuf(int n, int npp)
{
    cellp = -1;
}

void vtc_get_force_tree(Forceinfo *fi, Nbodyinfo *nb)
{
    Cell root;
    int i;

    fflush(stderr);
    vtc_print_cputime("create_tree start at");
    root = create_tree(fi, nb);
    fflush(stderr);
    PR(cellp, d);
    vtc_print_cputime("process_m2m start at");
    process_m2m(root, fi, nb);
#if 0
    if (sanity_check(root, 0, nb)) {
	fprintf(stderr, "sanity_check failed\n");
	exit(1);
    }
    else {
	fprintf(stderr, "sanity_check passed\n");
	exit(1);
    }
#endif
    vtc_print_cputime("<<< traverse_tree start at");
    traverse_tree(root, root, fi, nb);
}

#if 0
#define Sfprintf
#else
#define Sfprintf fprintf
#endif

static void
ps(FILE *fp, int len)
{
    int i;
    for (i = 0; i < len; i++) {
	Sfprintf(fp, "  ");
    }
}

static int
inbox(Cell c, BHindex b, Nbodyinfo *nb)
{
    int k;
    int isin = 1;
    double len = Csize[c]*0.5;
    double **x = nb->x;

    for (k = 0; k < 3; k++) {
	Cfprintf(stderr, ">>> %6.5f < %6.5f < %6.5f ?\n",
		 Cpos[k][c]-len, x[k][b], Cpos[k][c]+len);
	if (Cpos[k][c]+len < x[k][b]) {
	    isin = 0;
	}
	if (Cpos[k][c]-len > x[k][b]) {
	    isin = 0;
	}
    }
    return (isin);
}

static int
sanity_check(Cell c0, int level, Nbodyinfo *nb)
{
    int i, k, s;
    int err = 0;
    static ntot = 0;

    if (level == 0) {
	ntot = 0;
    }
    ps(stderr, level);
    Sfprintf(stderr, "cell%d size: %f pos: %6.5f %6.5f %6.5f\n",
	    c0, Csize[c0], Cpos[0][c0], Cpos[1][c0], Cpos[2][c0]);
    ps(stderr, level);
    Sfprintf(stderr, "cmmass: %f cmpos: %6.5f %6.5f %6.5f ndescendant: %d\n",
	    Ccmmass[c0], Ccmpos[0][c0], Ccmpos[1][c0], Ccmpos[2][c0], Cndescendant[c0]);
    if (Cisleaf[c0]) {
	/* check if all particles are in the cell */
	BHindex * bp0 = Cbody[c0];
	int cn0 = Cnbody[c0];

	ps(stderr, level);
	Sfprintf(stderr, "is a leaf (Cnbody: %d)\n", cn0);
	for (i = 0; i < cn0; i++) {
	    BHindex b = Cbody[c0][i];
	    ps(stderr, level);
	    Sfprintf(stderr, "b: %d  x: %6.5f %6.5f %6.5f",
		    b, nb->x[0][b], nb->x[1][b], nb->x[2][b]);
	    if (!inbox(c0, b, nb)) {
		fprintf(stderr, " NG");
		err++;
	    }
	    Sfprintf(stderr, "\n");
	}
	ntot += cn0;
	ps(stderr, level);
        Sfprintf(stderr, "---- ntot: %d\n", ntot);
    }
    else {
	ps(stderr, level);
	Sfprintf(stderr, "is not a leaf");
	for (s = 0; s < 8; s++) {
	    Cell c1 = Cchild[s][c0];
	    if (NOCELL == c1) {
		continue;
	    }
	    if (Csize[c0]*0.5 != Csize[c1]) {
		fprintf(stderr, " NG (Csize)");
		err++;
	    }
	    for (k = 0; k < 3; k++) {
		if (fabs(Cpos[k][c0]-Cpos[k][c1]) != Csize[c1]*0.5) {
		    fprintf(stderr, " NG (Cpos)");
		    err++;
		}
	    }
	    Sfprintf(stderr, "\n");
	    err += sanity_check(c1, level+1, nb);
	}
    }
    if (level == 0) {
	if (nb->n != ntot) {
	    fprintf(stderr, "ntot == %d != n == %d\n", ntot, nb->n);
	    exit(1);
	}
	else {
	    fprintf(stderr, "ntot == n == %d\n", ntot);
	}
    }
    return err;
}
