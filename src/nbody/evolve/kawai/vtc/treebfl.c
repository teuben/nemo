/*
 * create tree structure using breadth-first algorithm
 * in Makefile, set BFLOAD to 1 to activate this part
 *
 * probably this is slower than depth-first algorithm
 * on most platforms
 */

/* breadth-first construction */
static Cell /* root cell */
create_tree(Forceinfo *fi, Nbodyinfo *nb)
{
    int root;
    int i;
    int ncelltmp; /* # of cells in the level in question */
    int ncelltmp_new; /* # of cells to be processed in the next level */
    int nbody; /* # of particles in the level in question */
    int nbody_new; /* # of particles to be processed in the next level */
    double theta2;
    static int bbufsize = 0;
    static Body *bodyp;
#if DYNMEM
    static Body *body = NULL; /* particles in the level in question */
    static Body *body_new = NULL; /* particles to be processed in the next level */
    static Cell *bpid = NULL; /* particle's parent cell */
    static int *bsubid = NULL; /* bpis's child cell which the particle belongs to */
#else /* !DYNMEM */
    static Body body[NBODYMAX];
    static Body body_new[NBODYMAX];
    static Cell bpid[NBODYMAX];
    static int bsubid[NBODYMAX];
#endif /* DYNMEM */

    root = create_root_cell(fi, nb);
    Cisleaf[root] = FALSE;
    if (bbufsize < nb->n) {
#if DYNMEM
	size_t isize = sizeof(int)*(nb->n+1);
	size_t bsize = sizeof(Body)*(nb->n+1);
	size_t csize = sizeof(Cell)*(nb->n+1);
	size_t ksize = sizeof(Mortonkey)*(nb->n+1);

	Bindex = (Body *)realloc(Bindex, bsize);
	Bkey = (Mortonkey *)realloc(Bkey, ksize);
	body = (Body *)realloc(body, bsize);
	body_new = (Body *)realloc(body, bsize);
	bpid = (Cell *)realloc(bpid, csize);
	bsubid = (int *)realloc(bsubid, isize);
	if (NULL == bpid || NULL == body || NULL == bsubid) {
	    perror("create_tree");
	    exit(1);
	}
#else /* !DYNMEM */
    if (nb->n > NBODYMAX) {
	fprintf(stderr, "# of reqired body (%d) exceeds NBODYMAX (%d)\n",
		nb->n, NBODYMAX);
	exit(1);
    }
#endif /* DYNMEM */
	bbufsize = nb->n;
	Cfprintf(stderr, "create_tree buf reallocation done.\n");
    }
    ncelltmp = 1;
    celltmp[0] = root;
    nbody = nb->n;
    for (i = 0; i < nbody; i++) {
	bpid[i] = root;
	body[i] = i;
    }
    bodyp = Bindex;
    create_key(body, Bkey, nb, Csize[root]);

    /*
    for (i = 0; i < nbody; i++) {
	fprintf(stderr, "before i: %d key: %016llx x: %6.5f %6.5f %6.5f\n",
		i, Bkey[i], nb->x[0][i], nb->x[1][i], nb->x[2][i]);
    }
    */

    sort_body(body, Bkey, nbody);

    /*
    for (i = 0; i < nbody; i++) {
	fprintf(stderr, "after  i: %d key: %016llx x: %6.5f %6.5f %6.5f\n",
		i, Bkey[i], nb->x[0][i], nb->x[1][i], nb->x[2][i]);
    }
    */

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

	if (fi->test_id == 0) {
	    viewtree_test(nb);
	}
    }
    return (root);
}

static void
assign_bodies(int n, /* # of bodies */
	      Body *body, /* bodies to be assiged */
	      Cell *bpid, /* cell the body currently assigned to */
	      int *bsubid, /* subid of the child of the cell to which
			      the body should be assigned */
	      double (*bpos)[3]) /* the body position */
{
    int i, sid;

    for (i = 0; i < n; i++) {
	sid = 0;
	if (bpos[body[i]][0] > Cpos[bpid[i]][0]) {
	    sid += 1;
	}
	if (bpos[body[i]][1] > Cpos[bpid[i]][1]) {
	    sid += 2;
	}
	if (bpos[body[i]][2] > Cpos[bpid[i]][2]) {
	    sid += 4;
	}
	bsubid[i] = sid;
    }
}

static void
count_particles(int n, /* # of bodies */
		Cell *bpid, /* cell the body currently assigned to */
		int *bsubid) /* subid of the child of the cell to which
				the body should be assigned */
{
#if 1
    int i, is, ip;

    for (i = 0; i < n;) {
	is = bsubid[i];
	ip = bpid[i];
	while (bsubid[i] == is && bpid[i] == ip) {
	    Cchildnbody[ip][is]++;
	    i++;
	}
    }
#else
    int i;

    for (i = 0; i < n; i++) {
	Cchildnbody[bpid[i]][bsubid[i]]++;
    }
#endif
}

static int /* # of cells newly opened */
open_cells(int ndiv, /* min # of particle to process cell division */
	   int ncell0,
#if DYNMEM
	   Cell **cell0, /* cells in question */
	   Cell **cell1, /* cells newly opened */
#else /* !DYNMEM */
	   Cell (*cell0)[NCELLMAX], /* cells in question */
	   Cell (*cell1)[NCELLMAX], /* cells newly opened */
#endif /* DYNMEM */
	   Body **bodypp, /* points to Bindex[] in use */
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
	    cn1 = Cchildnbody[c0][s];
	    if (0 == cn1) {
		continue; /* empty cell. do not process anymore */
	    }
	    /* open this cell */
	    Cchild[c0][s] = c1 = make_cell();
	    Csize[c1] = newsize;
	    Ccrit[c1] = newsize*newsize/theta2;
	    for (k = 0; k < 3; k++) {
		Cpos[c1][k] = Cpos[c0][k] + offset[k];
	    }
	    if (cn1 > ndiv) {
		(*cell1)[ncell1] = c1;
		ncell1++;
		Cbody[c1] = *bodypp;
		Cisleaf[c1] = FALSE;
	    }
	    else {
		Cbody[c1] = *bodypp;
		(*bodypp) += cn1;
	    }
	}
    }
    return (ncell1);
}

static int /* # of bodies to be processed in the next level */
add_particles(int n0, /* # of bodies to be processed */
	      Body *body0, /* bodies to be processed */
	      Body *body1, /* bodies left to be processed in the next level */
	      Cell *bpid, /* cell the body currently assigned to */
	      int *bsubid) /* subid of the child of the cell to which
			      the body should be added */
{
    int i, j;
    int n1;
    Cell c0, c1;
    static int *processed = NULL;
    static n0save = 0;

    if (n0save < n0) {
	n0save = n0;
	processed = (int *)realloc(processed, sizeof(int)*n0);
    }
    if (NULL == processed) {
	perror("add_particle:");
	exit(1);
    }
    for (i = 0; i < n0; i++) {
	processed[i] = TRUE;
    }
    n1 = 0;
    for (i = 0; i < n0;) {
	c0 = bpid[i];
	c1 = Cchild[c0][bsubid[i]];

	if (Cisleaf[c1]) {
	    int nchild = Cnbody[c1] = Cchildnbody[c0][bsubid[i]];
	    for (j = 0; j < nchild; j++) {
		Cbody[c1][j] = body0[i+j];
	    }
	    i += nchild;
	}
	else {
	    processed[i] = FALSE;
	    i++;
	}
    }
    for (i = 0; i < n0; i++) {
	if (!processed[i]) {
	    c0 = bpid[i];
	    c1 = Cchild[c0][bsubid[i]];
	    body1[n1] = body0[i]; /* prepare bodies for the next level */
	    bpid[n1] = c1;
	    n1++;
	}
    }
    return(n1);
}

static int /* length of the list */
get_descendant_list(Cell c0, Forceinfo *fi, Nbodyinfo *nb, Body *list)
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
	    c1 = Cchild[c0][s];
	    if (NOCELL == c1) {
		continue;
	    }
	    len0 += get_descendant_list(c1, fi, nb, list+len0);
	}
    }
    return (len0);
}

#define LEVELMAX (32)
static void
m2m_1storder(Cell c0, /* dummy */
	     Forceinfo *fi, Nbodyinfo *nb)
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
	Ccmpos[c][0] = 0.0;
	Ccmpos[c][1] = 0.0;
	Ccmpos[c][2] = 0.0;
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

    cmp0 = (double *)&(Ccmpos[0][0]);
    cmp1 = (double *)&(Ccmpos[0][1]);
    cmp2 = (double *)&(Ccmpos[0][2]);
    nbm = nb->m;
    nbx0 = (double *)&(nb->x[0][0]);
    nbx1 = (double *)&(nb->x[0][1]);
    nbx2 = (double *)&(nb->x[0][2]);
    for (lv = level; lv >= 0; lv--) {
	/* leaf cells */
	for (c = first[lv]; c < first[lv+1]; c++) {
	    int i;
	    cb = Cbody[c];
	    if (!Cisleaf[c]) {
		continue;
	    }
	    for (i = 0; i < Cnbody[c]; i++) {
		Body b;
		b = cb[i];
		m1 = nbm[b];
		cmp0[c*3] += nbx0[b*3]*m1;
		cmp1[c*3] += nbx1[b*3]*m1;
		cmp2[c*3] += nbx2[b*3]*m1;
		Ccmmass[c] += m1;
	    }
	    Cndescendant[c] = Cnbody[c];
	}

	/* non-leaf cells */
	for (s = 0; s < 8; s++) {
	    int *ccs;
	    ccs = (int *)&(Cchild[0][s]);
#pragma loop novrec
	    for (c = first[lv]; c < first[lv+1]; c++) {
		Cell c1;
		c1 = ccs[c*8];
		if (Cisleaf[c]) {
		    continue;
		}
		if (NOCELL == c1) {
		    continue;
		}
		m1 = Ccmmass[c1];
		cmp0[c*3] += cmp0[c1*3]*m1;
		cmp1[c*3] += cmp1[c1*3]*m1;
		cmp2[c*3] += cmp2[c1*3]*m1;
		Ccmmass[c] += m1;
		Cndescendant[c] += Cndescendant[c1];
	    }
	}

	for (c = first[lv]; c < first[lv+1]; c++) {
	    cmp0[c*3] /= Ccmmass[c];
	    cmp1[c*3] /= Ccmmass[c];
	    cmp2[c*3] /= Ccmmass[c];
	}
    }
}
