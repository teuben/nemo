#pragma global noalias
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "vtc.h"
#include "vtclocal.h"

/*
 * local defs
 */
static int cellp = -1; /* # of cells */
static int cellpmax = 0; /* size of cell buffers */
static int usep2m2; /* FALSE:treecode    TRUE:P2M2 treecode */
static int negativemass; /* if TRUE, partices may have negative mass */
#if DYNMEM /* (dynamical memory allocation) */
static Body *Bindex = NULL;
static Mortonkey *Bkey = NULL;
static double (*Cpos)[3] = NULL;
static double (*Ccmpos)[3] = NULL;
static double *Ccmmass = NULL;
static double (*Ccmpos2)[3] = NULL; /* for negative mass */
static double *Ccmmass2 = NULL; /* for negative mass */
static double *Csize = NULL;
static Cell (*Cchild)[8] = NULL;
static int *Cisleaf = NULL;
static int *Cnbody = NULL;
static int *Cndescendant = NULL;
static int *Cnpositive = NULL;
static int *Cnnegative = NULL;
static Body **Cbody = NULL;
static double *Ccrit = NULL;
static double (*Cpppos)[3] = NULL;
static double *Cppmass = NULL;
#if BFLOAD
static int (*Cchildnbody)[8] = NULL;
static Cell *celltmp = NULL; /* cells in the level in question */
static Cell *celltmp_new = NULL; /* cells in the next level */
#endif /* BFLOAD */
#else /* !DYNMEM (static memory allocation) */
#define NCELLMAX (2000000)
#define NBODYMAX (2000000)
static Body Bindex[NBODYMAX];
static Mortonkey Bkey[NBODYMAX];
static double Cpos[NCELLMAX][3];
static double Ccmpos[NCELLMAX][3];
static double Ccmmass[NCELLMAX];
static double Ccmpos2[NCELLMAX][3]; /* for negative mass */
static double Ccmmass2[NCELLMAX]; /* for negative mass */
static double Csize[NCELLMAX];
static Cell Cchild[NCELLMAX][8];
static int Cisleaf[NCELLMAX];
static int Cnbody[NCELLMAX];
static int Cndescendant[NCELLMAX];
static int Cnpositive[NCELLMAX];
static int Cnnegative[NCELLMAX];
static Body *Cbody[NCELLMAX];
static double Ccrit[NCELLMAX];
static double Cpppos[NCELLMAX][3];
static double Cppmass[NCELLMAX];
#if BFLOAD
static int Cchildnbody[NCELLMAX][8];
static Cell celltmp[NCELLMAX];
static Cell celltmp_new[NCELLMAX];
#else /* !BFLOAD */
static void create_tree_dfl(Cell c0, int n0, Mortonkey key0, int keylen, int ndiv, double theta2, Nbodyinfo *nb);
#endif /* BFLOAD */
#endif /* DYNMEM */

/* prototype defs of local funcs */
static int create_root_cell(Forceinfo *fi, Nbodyinfo *nb);
static void create_key(Body *a, Mortonkey *key, Nbodyinfo *nb, double rsize);
static Cell create_tree(Forceinfo *fi, Nbodyinfo *nb);
static void m2m_1storder(Cell c0, Forceinfo *fi, Nbodyinfo *nb);
static void m2m_1storder_md(Cell c0, Forceinfo *fi, Nbodyinfo *nb);
static void printmatrix(double a[3][3]);
static int calc_pppos_2ndorder(double q[3][3], double cmpos[3], double cmmass, double pppos[3][3], int c0);
static void err_2ndorder(Nbodyinfo *nb, int sign, int c0, double q[3][3],
			 double xx, double yy, double zz, double xy, double yz, double zx);
static void copy_pppos_2ndorder(Nbodyinfo *nb, Cell c0, int off, int nparticle, int sign, int npp);
static void m2m_2ndorder(Cell c0, Forceinfo *fi, Nbodyinfo *nb);
static void m2m_2ndorder_md(Cell c0, Forceinfo *fi, Nbodyinfo *nb);
static void m2m_anyorder(Cell c0, Forceinfo *fi, Nbodyinfo *nb);
static void m2m_anyorder_md(Cell c0, Forceinfo *fi, Nbodyinfo *nb);
static void process_m2m(Cell root, Forceinfo *fi, Nbodyinfo *nb);
static int is_overlapped(Cell icell, Cell jcell);
static double get_separation(double isize, double ipos[3], Cell jcell);
static int get_interaction_list(Cell icell, Cell jcell0, double theta2, Nbodyinfo *nb,
				Body *list, int *iscell, int npp);
static void traverse_tree(Cell c0, Cell root, Forceinfo *fi, Nbodyinfo *nb);
static void reallocate_cellbuf(int npp);
static int make_cell(void);
static void reset_cellbuf(int n, int npp);
static void viewtree_test(Nbodyinfo *nb);
static void viewlist_test(Nbodyinfo *nb, int ni, int nj, Cell c0, Body *nodelist, int *iscell);
static void ps(FILE *fp, int len);
static int inbox(Cell c, Body b, Nbodyinfo *nb);
static int sanity_check(Cell c0, int level, Nbodyinfo *nb);
static int pack_interaction_list(Nbodyinfo *nb, int jlen, int npp, int usenegativemass,
				 Node *nodelist, int *iscell, double *mlist, double (*xlist)[3]);
#if BFLOAD
static void assign_bodies(int n, Body *body, Cell *bpid, int *bsubid, double (*bpos)[3]);
static void count_particles(int n, Cell *bpid, int *bsubid);
#if DYNMEM
static int open_cells(int ndiv, int ncell0, Cell **cell0, Cell **cell1, Body **bodypp, double theta2);
#else /* !DYNMEM */
static int open_cells(int ndiv, int ncell0, Cell (*cell0)[NCELLMAX], Cell (*cell1)[NCELLMAX], Body **bodypp, double theta2);
#endif /* DYNMEM */
static int add_particles(int n0, Body *body0, Body *body1, Cell *bpid, int *bsubid);
static int get_descendant_list(Cell c0, Forceinfo *fi, Nbodyinfo *nb, Body *list);
#else /* !BFLOAD */
static void create_tree_dfl(Cell c0, int n0, Mortonkey key0, int keylen, int ndiv, double theta2, Nbodyinfo *nb);
#endif /* BFLOAD */

/*
 * global functions
 */
/* 'typical' params for tree algorithm on GRAPE */
void
vtc_get_default_tree_params(Forceinfo *tr)
{
    tr->theta = 0.75;
    tr->eps = 0.02;
    tr->ncrit = 8192;
    tr->node_div_crit = 8;
    tr->p = 1;
    tr->negativemass = FALSE;
    tr->pprad = 0.2;
    tr->npp = 1;
    tr->pppos = NULL;
    tr->full_dof = TRUE;
    tr->calculator = HOST;
    tr->test_id = -1;
}

void vtc_get_force_tree(Forceinfo *fi, Nbodyinfo *nb)
{
    Cell root;
    int i, k;
    static int firstcall = 1;

    if (firstcall) {
	firstcall = 0;
	negativemass = fi->negativemass;
	if (fi->p == 1 && fi->full_dof) {
	    usep2m2 = FALSE;
	}
	else {
	    usep2m2 = TRUE;
	    load_design(fi->p, fi->pprad, fi->full_dof, fi->negativemass,
			&(fi->npp), &(fi->pppos));
	}
#if 0
	fprintf(stderr, "p: %d pprad: %f full_dof: %d npp: %d\n",
		fi->p, fi->pprad, fi->full_dof, fi->npp);
	if (fi->p > 2 || !(fi->full_dof)) {
	    for (i = 0; i < fi->npp; i++) {
		fprintf(stderr, "pp%03d: % 20.18f % 20.18f % 20.18f\n",
			i, fi->pppos[i][0], fi->pppos[i][1], fi->pppos[i][2]);
	    }
	}
#endif
    }

    vtc_print_cputime("create_tree start at");
    root = create_tree(fi, nb);
/*
    PR(cellp, d); PRL(cellpmax, d);
    */
    vtc_print_cputime("process_m2m start at");
    process_m2m(root, fi, nb);
#if 0
    if (sanity_check(root, 0, nb)) {
	fprintf(stderr, "sanity_check failed\n");
	exit(1);
    }
    else {
	fprintf(stderr, "sanity_check passed\n");
    }
#endif
    if (fi->test_id == 1) {
	viewtree_test(nb);
    }
    vtc_print_cputime("traverse_tree start at");
    traverse_tree(root, root, fi, nb);
}

/*
 * local functions
 */
static int
create_root_cell(Forceinfo *fi, Nbodyinfo *nb)
{
    int root;
    int i, k;
    double cs, val, posmax = 0.0;
    static double csold = 0.0;
    static double mmin;
    static int cnt = 0;
    static int firstcall = TRUE;

    if (firstcall) { /* necessary only for GRAPE-3&5 */
	firstcall = FALSE;
	mmin = nb->m[0];
	for (i = 0; i < nb->n; i++) {
	    if (nb->m[i] < mmin) {
		mmin = nb->m[i];
	    }
	}
	mmin /= 3.0; /* to resolve pseudoparticle mass in quadrupole expansion */
    }
    reset_cellbuf(nb->n, fi->npp);
    root = make_cell();
    cs = 1.0;
    posmax = 0.0;
    for (i = 0; i < nb->n; i++) {
	for (k = 0; k < 3; k++) {
	    val = fabs(nb->x[i][k]);
	    if (posmax < val) {
		posmax = val;
	    }
	}
    }
    while (posmax*2 > cs) {
	cs *= 2.0;
    }
    while (posmax*2 <= cs*0.5) {
	cs *= 0.5;
    }
    if (cs != csold) {/* necessary only for GRAPE-3,5 and MDGRAPE-2 */
	csold = cs;
	vtc_set_scale(cs*2.0, mmin);
	fprintf(stderr, "rescaled. cnt: %d xmax: %f mmin: %f\n",
		cnt, cs*2.0, mmin);
    }
    Cbody[root] = Bindex;
    Csize[root] = cs;
    Ccrit[root] = cs*cs/(fi->theta)/(fi->theta);
    for (k = 0; k < 3; k++) {
	Cpos[root][k] = 0.0;
    }
    Cnbody[root] = nb->n;
    cnt++;
    return (root);
}

#if (MORTONKEY_LENGTH == 16) /* hand tuned for key length 16 */
static void
create_key(Body *a, Mortonkey *key, Nbodyinfo *nb, double rsize)
{
    int i, k, b, b3;
    Mortonkey srcval, dstval;
    double rscale = 2.0/rsize*MORTONKEY_OFFSET;
    int n = nb->n;
    double *nbx0, *nbx1, *nbx2;

    for (i = 0; i < n; i++) {
	key[i] = 0ll;
    }
    nbx0 = (double *)&(nb->x[0][0]);
    nbx1 = (double *)&(nb->x[0][1]);
    nbx2 = (double *)&(nb->x[0][2]);
    for (i = 0; i < n; i++) {
	/* k = 0 */
	dstval = 0;
	srcval = (Mortonkey)(nbx0[a[i]*3] * rscale+MORTONKEY_OFFSET);
	if (srcval & 0x0001) dstval |= 0x0000000000000001ll;
	if (srcval & 0x0002) dstval |= 0x0000000000000008ll;
	if (srcval & 0x0004) dstval |= 0x0000000000000040ll;
	if (srcval & 0x0008) dstval |= 0x0000000000000200ll;
	if (srcval & 0x0010) dstval |= 0x0000000000001000ll;
	if (srcval & 0x0020) dstval |= 0x0000000000008000ll;
	if (srcval & 0x0040) dstval |= 0x0000000000040000ll;
	if (srcval & 0x0080) dstval |= 0x0000000000200000ll;
	if (srcval & 0x0100) dstval |= 0x0000000001000000ll;
	if (srcval & 0x0200) dstval |= 0x0000000008000000ll;
	if (srcval & 0x0400) dstval |= 0x0000000040000000ll;
	if (srcval & 0x0800) dstval |= 0x0000000200000000ll;
	if (srcval & 0x1000) dstval |= 0x0000001000000000ll;
	if (srcval & 0x2000) dstval |= 0x0000008000000000ll;
	if (srcval & 0x4000) dstval |= 0x0000040000000000ll;
	if (srcval & 0x8000) dstval |= 0x0000200000000000ll;
	key[i] |= dstval<<2;

	/* k = 1 */
	dstval = 0;
	srcval = (Mortonkey)(nbx1[a[i]*3] * rscale+MORTONKEY_OFFSET);
	if (srcval & 0x0001) dstval |= 0x0000000000000001ll;
	if (srcval & 0x0002) dstval |= 0x0000000000000008ll;
	if (srcval & 0x0004) dstval |= 0x0000000000000040ll;
	if (srcval & 0x0008) dstval |= 0x0000000000000200ll;
	if (srcval & 0x0010) dstval |= 0x0000000000001000ll;
	if (srcval & 0x0020) dstval |= 0x0000000000008000ll;
	if (srcval & 0x0040) dstval |= 0x0000000000040000ll;
	if (srcval & 0x0080) dstval |= 0x0000000000200000ll;
	if (srcval & 0x0100) dstval |= 0x0000000001000000ll;
	if (srcval & 0x0200) dstval |= 0x0000000008000000ll;
	if (srcval & 0x0400) dstval |= 0x0000000040000000ll;
	if (srcval & 0x0800) dstval |= 0x0000000200000000ll;
	if (srcval & 0x1000) dstval |= 0x0000001000000000ll;
	if (srcval & 0x2000) dstval |= 0x0000008000000000ll;
	if (srcval & 0x4000) dstval |= 0x0000040000000000ll;
	if (srcval & 0x8000) dstval |= 0x0000200000000000ll;
	key[i] |= dstval<<1;

	/* k = 2 */
	dstval = 0;
	srcval = (Mortonkey)(nbx2[a[i]*3] * rscale+MORTONKEY_OFFSET);
	if (srcval & 0x0001) dstval |= 0x0000000000000001ll;
	if (srcval & 0x0002) dstval |= 0x0000000000000008ll;
	if (srcval & 0x0004) dstval |= 0x0000000000000040ll;
	if (srcval & 0x0008) dstval |= 0x0000000000000200ll;
	if (srcval & 0x0010) dstval |= 0x0000000000001000ll;
	if (srcval & 0x0020) dstval |= 0x0000000000008000ll;
	if (srcval & 0x0040) dstval |= 0x0000000000040000ll;
	if (srcval & 0x0080) dstval |= 0x0000000000200000ll;
	if (srcval & 0x0100) dstval |= 0x0000000001000000ll;
	if (srcval & 0x0200) dstval |= 0x0000000008000000ll;
	if (srcval & 0x0400) dstval |= 0x0000000040000000ll;
	if (srcval & 0x0800) dstval |= 0x0000000200000000ll;
	if (srcval & 0x1000) dstval |= 0x0000001000000000ll;
	if (srcval & 0x2000) dstval |= 0x0000008000000000ll;
	if (srcval & 0x4000) dstval |= 0x0000040000000000ll;
	if (srcval & 0x8000) dstval |= 0x0000200000000000ll;
	key[i] |= dstval;
    }
}

#else /* MORTONKEY_LENGTH != 16 */

static void
create_key(Body *a, Mortonkey *key, Nbodyinfo *nb, double rsize)
{
    int i, k, b, b3;
    Mortonkey srcval, dstval;
    double rscale = 1.0/rsize*MORTONKEY_OFFSET;
    int n = nb->n;

    for (i = 0; i < n; i++) {
	dstval = 0;
	for (k = 0; k < 3; k++) {
	    srcval = (Mortonkey)(nb->x[a[i]][k] * rscale + MORTONKEY_OFFSET);
	    /*
	    fprintf(stderr, "x: %f srcval: %016llx\n", nb->x[a[i]][k], srcval);
	    */
	    /* now dstval is in the range [0, 2*MORTONKEY_OFFSET] */
	    for (b = 0, b3 = (2-k); b < MORTONKEY_LENGTH; b++, b3 += 3) {
		if (srcval & (1<<b)) {
		    dstval |= ((Mortonkey)1)<<b3;
		}
	    }
	}
	key[i] = dstval;
	/*
	fprintf(stderr, "key[%d]: %016llx\n\n", i, key[i]);
	*/
    }
}
#endif /* MORTONKEY_LENGTH == 16 */

#if BFLOAD
#include "treebfl.c"
#else /* !BFLOAD */

/* create tree structure using depth-first algorithm */
static Cell /* root cell */
create_tree(Forceinfo *fi, Nbodyinfo *nb)
{
    int i;
    int n = nb->n;
    Cell root;
    static int bufsize;

    if (bufsize < n) {
#if DYNMEM
	Bindex = (Body *)realloc(Bindex, sizeof(Body)*n+10);
	Bkey = (Mortonkey *)realloc(Bkey, sizeof(Mortonkey)*n+10);
	if (NULL == Bindex || NULL == Bkey) {
	    perror("create_tree");
	    exit(1);
	}
#else /* !DYNMEM */
	if (n > NBODYMAX) {
	    fprintf(stderr, "# of reqired body (%d) exceeds NBODYMAX (%d)\n",
		    nb->n, NBODYMAX);
	    exit(1);
	}
#endif /* DYNMEM */
	bufsize = n;
    }
    for (i = 0; i < n; i++) {
	Bindex[i] = i;
    }
    root = create_root_cell(fi, nb);
    create_key(Bindex, Bkey, nb, Csize[root]);
    sort_body(Bindex, Bkey, n);
    create_tree_dfl(root, n, 0, MORTONKEY_LENGTH,
		    fi->node_div_crit, (fi->theta)*(fi->theta), nb);
    return(root);
}

static void
create_tree_dfl(Cell c0, /* cell in question */
		int n0, /* # of bodies in c0 */
		Mortonkey key0, /* c0's key */
		int keylen, /* min non-zero c0's key */
		int ndiv, /* min # of particle to process cell division */
		double theta2, Nbodyinfo *nb)
{
    int s, k, n1, btop, bend, btopsave;
    Cell c1;
    Mortonkey keyscale, key1;
    double size1;

    if (n0 <= ndiv || keylen == 0) {
	return;
    }
    keyscale = ((Mortonkey)1)<<((keylen-1)*3);
    size1 = Csize[c0]*0.5;
    btop = Cbody[c0]-Bindex; /* top of bodies in question */
    for (s = 0; s < 8; s++) {
	key1 = key0 + keyscale * s;
	if (Bkey[btop]-key1 >= keyscale) {
	    continue; /* no body belongs to child-cell s */
	}
	Cisleaf[c0] = FALSE;
	Cchild[c0][s] = c1 = make_cell();
	Cbody[c1] = Bindex+btop;
	Csize[c1] = size1;
	Ccrit[c1] = size1*size1/theta2;
	for (k = 0; k < 3; k++) {
	    int mask = 1<<(2-k);
	    Cpos[c1][k] = Cpos[c0][k] + ((s&mask)/mask-0.5)*size1;
	}
	/* count # of particles to be assigned to s */
	btopsave = btop;
	bend = btop + n0;
	while (bend - btop > 1) {
	    int b;
	    b = (bend+btop)/2;
	    if (Bkey[b]-key1 >= keyscale) {
		bend = b;
	    }
	    else {
		btop = b;
	    }
	}
	btop = bend;
	Cnbody[c1] = n1 = bend-btopsave;
	create_tree_dfl(c1, n1, key1, keylen-1, ndiv, theta2, nb);
	n0 -= n1;
#if 0
	{
	    int i;
	    
	    fprintf(stderr, "\nc1: %d n1: %d key: %016llx size: %f pos: %f %f %f\n",
		    c1, n1, key1, size1, Cpos[c1][0], Cpos[c1][1], Cpos[c1][2]);
	    for (i = 0; i < n1; i++) {
		Body b = Cbody[c1][i];
		fprintf(stderr, "b: %d key: %016llx x: %6.5f %6.5f %6.5f",
			b, Bkey[Cbody[c1]+i-Bindex], nb->x[b][0], nb->x[b][1], nb->x[b][2]);
		if (!inbox(c1, b, nb)) {
		    fprintf(stderr, " NG");
		}
		fprintf(stderr, "\n");
	    }
	}
#endif
	if (n0 <= 0) {
	    return;
	}
    }
}

static void
m2m_1storder(Cell c0, Forceinfo *fi, Nbodyinfo *nb)
{
    int i, k, s;
    Cell c1;
    Body b;
    double m1;

    for (k = 0; k < 3; k++) {
	Ccmpos[c0][k] = 0.0;
    }
    Ccmmass[c0] = 0.0;
    Cndescendant[c0] = 0;
    if (Cisleaf[c0]) {
	for (i = 0; i < Cnbody[c0]; i++) {
	    b = Cbody[c0][i];
	    m1 = nb->m[b];
	    for (k = 0; k < 3; k++) {
		Ccmpos[c0][k] += (nb->x[b][k])*m1;
	    }
	    Ccmmass[c0] += m1;	}
	Cndescendant[c0] = Cnbody[c0];
    }
    else {
	for (s = 0; s < 8; s++) {
	    c1 = Cchild[c0][s];
	    if (NOCELL == c1) {
		continue;
	    }
	    m2m_1storder(c1, fi, nb);
	    m1 = Ccmmass[c1];
	    for (k = 0; k < 3; k++) {
		Ccmpos[c0][k] += Ccmpos[c1][k]*m1;
	    }
	    Ccmmass[c0] += m1;
	    Cndescendant[c0] += Cndescendant[c1];
	}
    }
    for (k = 0; k < 3; k++) {
	Ccmpos[c0][k] /= Ccmmass[c0];
    }
}

#endif /* BFLOAD */

#include "highorder.c"

static void
process_m2m(Cell root, Forceinfo *fi, Nbodyinfo *nb)
{
    if (fi->negativemass) { /* particle may have negative mass */
	if (fi->p > 2 || !(fi->full_dof)) {
	    m2m_1storder_md(root, fi, nb);
	    m2m_anyorder_md(root, fi, nb);
	}
	else if (1 == fi->p) {
	    m2m_1storder_md(root, fi, nb);
	}
	else if (2 == fi->p) {
	    m2m_1storder_md(root, fi, nb);
	    m2m_2ndorder_md(root, fi, nb);
	}
	else {
	    fprintf(stderr, "process_m2m: order %d invalid\n", fi->p);
	    exit(1);
	}
    }
    else { /* particle mass is always positive */
	if (fi->p > 2 || !(fi->full_dof)) {
	    m2m_1storder(root, fi, nb);
	    m2m_anyorder(root, fi, nb);
	}
	else if (1 == fi->p) {
	    m2m_1storder(root, fi, nb);
	}
	else if (2 == fi->p) {
	    m2m_1storder(root, fi, nb);
	    m2m_2ndorder(root, fi, nb);
	}
	else {
	    fprintf(stderr, "process_m2m: order %d invalid\n", fi->p);
	    exit(1);
	}
    }
}

/* currently not used */
static int
is_overlapped(Cell icell, Cell jcell)
{
    int k;
    double dx[3];
    double xmin = (Csize[icell]+Csize[jcell])*0.4999999999999999;

    for (k = 0; k < 3; k++) {
	dx[k] = Cpos[icell][k] - Cpos[jcell][k];
    }
    if ((fabs(dx[0]) > xmin) || (fabs(dx[1]) > xmin) || (fabs(dx[2]) > xmin)) {
	return (FALSE);
    }
    fprintf(stderr, "---- overlaped ----\n");
    return (TRUE);
}

#if VECTORIZED
#include "treev.c"
#else /* scalar */
static int /* # of nodes packed */
pack_interaction_list(Nbodyinfo *nb, int jlen, int npp, int usenegativemass,
		      Node *nodelist, int *iscell, double *mlist, double (*xlist)[3])
{
    int i, j;
    int pplen = 0;
    Node node;

    if (usep2m2) {
	for (j = 0; j < jlen; j++) {
	    node = nodelist[j];
	    if (iscell[j]) { /* node is an index to a cell */
		for (i = 0; i < npp; i++) {
		    mlist[pplen] = Cppmass[npp*node+i];
		    if (mlist[pplen] == 0.0) {
			continue;
		    }
		    xlist[pplen][0] = Cpppos[npp*node+i][0];
		    xlist[pplen][1] = Cpppos[npp*node+i][1];
		    xlist[pplen][2] = Cpppos[npp*node+i][2];
		    pplen++;
		}
	    }
	    else { /* node is an index to a particle */
		mlist[pplen] = nb->m[node];
		xlist[pplen][0] = nb->x[node][0];
		xlist[pplen][1] = nb->x[node][1];
		xlist[pplen][2] = nb->x[node][2];
		pplen++;
	    }
	}
    }
    else if (usenegativemass) { /* no P2M2 */
	for (j = 0; j < jlen; j++) {
	    node = nodelist[j];
	    if (iscell[j]) { /* node is an index to a cell */
		mlist[pplen] = Ccmmass[node];
		if (mlist[pplen] != 0.0) {
		    xlist[pplen][0] = Ccmpos[node][0];
		    xlist[pplen][1] = Ccmpos[node][1];
		    xlist[pplen][2] = Ccmpos[node][2];
		    pplen++;
		}
		mlist[pplen] = Ccmmass2[node];
		if (mlist[pplen] != 0.0) {
		    xlist[pplen][0] = Ccmpos2[node][0];
		    xlist[pplen][1] = Ccmpos2[node][1];
		    xlist[pplen][2] = Ccmpos2[node][2];
		    pplen++;
		}
	    }
	    else { /* node is an index to a particle */
		mlist[pplen] = nb->m[node];
		xlist[pplen][0] = nb->x[node][0];
		xlist[pplen][1] = nb->x[node][1];
		xlist[pplen][2] = nb->x[node][2];
		pplen++;
	    }
	}
    }
    else { /* no P2M2, no negative mass particle */
	for (j = 0; j < jlen; j++) {
	    node = nodelist[j];

	    if (iscell[j]) { /* node is an index to a cell */
		mlist[pplen] = Ccmmass[node];
		xlist[pplen][0] = Ccmpos[node][0];
		xlist[pplen][1] = Ccmpos[node][1];
		xlist[pplen][2] = Ccmpos[node][2];
		pplen++;
	    }
	    else { /* node is an index to a particle */
		mlist[pplen] = nb->m[node];
		xlist[pplen][0] = nb->x[node][0];
		xlist[pplen][1] = nb->x[node][1];
		xlist[pplen][2] = nb->x[node][2];
		pplen++;
	    }
	}
    }
    return (pplen);
}

static int /* length of the list excluding i-particles */
get_interaction_list(Cell icell, Cell jcell0, double theta2, Nbodyinfo *nb,
		     Body *list, int *iscell, int npp)
{
    int i, s;
    int len0 = 0;
    Cell jcell1;

    if (icell == jcell0) { /* self interaction is already taken into account */
	len0 = 0;
    }
    else if (get_separation(Csize[icell], Cpos[icell], jcell0) * theta2 > (Csize[jcell0])*(Csize[jcell0])) {
	/*
	  else if (get_separation(Csize[icell], Cpos[icell], jcell0) > Ccrit[jcell0]) {
	*/
	/* apply multipole expansion */
	list[0] = jcell0;
	iscell[0] = TRUE;
	len0 = 1;
    }
    else { /* not well separated */
	if (Cisleaf[jcell0] || Cndescendant[jcell0] <= npp) { /* handle particle force directly */
	    Body *cb = Cbody[jcell0];
	    len0 = Cnbody[jcell0];
	    for (i = 0; i < len0; i++) {
		list[i] = cb[i];
		iscell[i] = FALSE;
	    }
	}
	else { /* descend the tree */
	    for (s = 0; s < 8; s++) {
		jcell1 = Cchild[jcell0][s];
		if (NOCELL == jcell1) {
		    continue;
		}
		len0 += get_interaction_list(icell, jcell1, theta2, nb,
					     list+len0, iscell+len0, npp);
	    }
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
#if 1
    xmin = (isize)*0.5; /* i & j are guaranteed not to be overlaped */
#else
    /* you may want to use more strict criterion */
    xmin = (isize+jsize)*0.5;
#endif

    if (negativemass) { /* use distance from geometric center of jcell */
	dx = fabs(ipos[0]-Cpos[jcell][0]);
	dy = fabs(ipos[1]-Cpos[jcell][1]);
	dz = fabs(ipos[2]-Cpos[jcell][2]);
    }
    else {  /* use distance from center of mass of jcell.
	       this works more better for lower order expansion */
       dx = fabs(ipos[0]-Ccmpos[jcell][0]);
       dy = fabs(ipos[1]-Ccmpos[jcell][1]);
       dz = fabs(ipos[2]-Ccmpos[jcell][2]);
    }
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

#endif /* VECTORIZED */

static void
traverse_tree(Cell c0, Cell root, Forceinfo *fi, Nbodyinfo *nb)
{
    Cell c1;
    int node; /* index to a cell or a particle */
    int i, j, k, s;
    int ilen, jlen, pplen;
    int npp = fi->npp;
    double epsinv;
    double *nba0, *nba1, *nba2, *nbm, *nbp;
    static int pplenmax = 0;
    static double amax = 0.0;
    static int nwalk;
    static double pptot;
    static int iscell[NJMAX]; /* node is a cell or a particle */
    static Body nodelist[NJMAX];
    static double mlist[NJMAX];
    static double xlist[NJMAX][3];
    static double plist[NJMAX];
    static double alist[NJMAX][3];

    if (c0 == root) {
	pptot = 0.0;
	nwalk = 0;
    }
    if (Cndescendant[c0] > fi->ncrit && !Cisleaf[c0]) {
	for (s = 0; s < 8; s++) {
	    c1 = Cchild[c0][s];
	    if (NOCELL == c1) {
		continue;
	    }
	    traverse_tree(c1, root, fi, nb);
	}
    }
    else {

#if BFLOAD
	ilen = get_descendant_list(c0, fi, nb, nodelist);
	if (ilen != Cndescendant[c0]) {
	    fprintf(stderr, "ilen (%d) differs from Cndescendant[%d] (%d)\n",
		    ilen, c0, Cndescendant[c0]);
	    exit(1);
	}
#else /* !BFLOAD */
	ilen = Cndescendant[c0];
	for (i = 0; i < ilen; i++) {
	    nodelist[i] = Cbody[c0][i];
	}
#endif /* BFLOAD */
	for (i = 0; i < ilen; i++) {
	    iscell[i] = FALSE;
	}
	jlen = ilen + get_interaction_list(c0, root, (fi->theta)*(fi->theta), nb,
					   nodelist+ilen, iscell+ilen, npp);
	pplen = pack_interaction_list(nb, jlen, npp, fi->negativemass, nodelist, iscell, mlist, xlist);
	if (pplen > NJMAX) {
	    fprintf(stderr, "too long interaction list (%d)\n", pplen);
	    exit(1);
	}
	if (pplen > pplenmax) {
	    pplenmax = pplen;
	}

	if (fi->test_id == 2) {
	    viewlist_test(nb, ilen, pplen, c0, nodelist, iscell);
	}
	(vtc_force_calculator[fi->calculator])(ilen, xlist,
					       pplen, xlist, mlist, fi->eps,
					       alist, plist);
	nwalk++;
	pptot += ilen*pplen;

	if (fi->eps != 0.0) { 
	    epsinv = 1.0/fi->eps;
	}
	else {
	    epsinv = 0.0;
	}
        if (fi->calculator != GRAPE_POTENTIALONLY &&
	    fi->calculator != HOST_POTENTIALONLY) {
	    nba0 = (double *)&(nb->a[0][0]);
	    nba1 = (double *)&(nb->a[0][1]);
	    nba2 = (double *)&(nb->a[0][2]);
	    for (i = 0; i < ilen; i++) {
		node = nodelist[i];
		nba0[node*3] = alist[i][0];
		nba1[node*3] = alist[i][1];
		nba2[node*3] = alist[i][2];
	    }
	}
        if (fi->calculator != GRAPE_FORCEONLY &&
	    fi->calculator != HOST_FORCEONLY) {
	    nbp = nb->p;
	    nbm = nb->m;
	    for (i = 0; i < ilen; i++) {
		node = nodelist[i];
		nbp[node] = plist[i] + nbm[node] * epsinv;
		/* MD2 does not calc self interaction
		   nbp[node] = plist[i];
		*/
		/* but on vpp5k this correction seems to be necessary
		 * I don't know why... */
	    }
	}
    }
    if (c0 == root) {
	fi->ninteraction = pptot;
	fi->nwalk = nwalk;
	fprintf(stderr, "pplenmax: %d amax: %f\n", pplenmax, amax);
    }
}

static void
reallocate_cellbuf(int npp)
{
#if DYNMEM
    size_t dsize = sizeof(double)*cellpmax;
    size_t isize = sizeof(int)*cellpmax;
    size_t csize = sizeof(Cell)*cellpmax;
    size_t bsize = sizeof(Body*)*cellpmax;
    static int cnt = 0;

#if BFLOAD
    celltmp = (Cell *)realloc(celltmp, csize);
    celltmp_new = (Cell *)realloc(celltmp_new, csize);
    Cchildnbody = (int (*)[8])realloc(Cchildnbody, isize*8);
#endif /* BFLOAD */
    Cnbody = (int *)realloc(Cnbody, isize);
    Cndescendant = (int *)realloc(Cndescendant, isize);
    Cbody = (Body **)realloc(Cbody, bsize);
    Cchild = (Cell (*)[8])realloc(Cchild, csize*8);
    Cisleaf = (int *)realloc(Cisleaf, isize);
    Cpos = (double (*)[3])realloc(Cpos, dsize*3);
    Ccmpos = (double (*)[3])realloc(Ccmpos, dsize*3);
    Ccmmass = (double *)realloc(Ccmmass, dsize);
    if (negativemass) {
	Ccmpos2 = (double (*)[3])realloc(Ccmpos2, dsize*3);
	Ccmmass2 = (double *)realloc(Ccmmass2, dsize);
	Cnpositive = (int *)realloc(Cnpositive, isize);
	Cnnegative = (int *)realloc(Cnnegative, isize);
    }
    Ccrit = (double *)realloc(Ccrit, dsize);
    Csize = (double *)realloc(Csize, dsize);
    if (usep2m2) {
	static int npp0 = 0;
	if (npp > npp0) {
	    npp0 = npp;
	}
	Cpppos = (double (*)[3])realloc(Cpppos, dsize*npp0*3);
	Cppmass = (double *)realloc(Cppmass, dsize*npp0);
    }
    if (NULL == Csize) {
	perror("reallocate_cellbuf");
	exit(1);
    }
    cnt++;
    fprintf(stderr, "reallocate_cellbuf done. cnt: %d cellpmax: %d\n",
	    cnt, cellpmax);
#else /* !DYNMEM */
    if (cellpmax > NCELLMAX) {
	fprintf(stderr, "# of reqired cell (%d) exceeds NCELLMAX (%d)\n",
		cellpmax, NCELLMAX);
	exit(1);
    }
#endif /* DYNMEM */
}

static Cell
make_cell(void)
{
    int s;

    cellp++;
    if (cellp >= cellpmax) {
	cellpmax = (cellpmax+1)*1.5;
	reallocate_cellbuf(0);
    }
    Cisleaf[cellp] = TRUE;
    Cnbody[cellp] = 0;
    Cbody[cellp] = NULL;
    for (s = 0; s < 8; s++) {
#if BFLOAD
	Cchildnbody[cellp][s] = 0;
#endif /* BFLOAD */
	Cchild[cellp][s] = NOCELL;
    }
    return (cellp);
}

static void
reset_cellbuf(int n, int npp)
{
    /* guess necessary # of cells from n & npp */
    if (cellpmax < (int)(n*0.6/npp)) {
	fprintf(stderr, "reset_cellbuf ");
	cellpmax = n*0.6/npp+100;
	reallocate_cellbuf(npp);
    }
    cellp = -1;
}

#include "debug.c"
