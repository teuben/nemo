/*
 * functions below are for debugging purpose only
 */

static void
viewtree_test(Nbodyinfo *nb)
{
#if USEX11
    vtc_viewtree(nb->n, nb->x, cellp, Cpos, Csize);
#endif /* USEX11 */
}

static void
viewlist_test(Nbodyinfo *nb, int ni, int nj, Cell c0,
	      Body *nodelist, int *iscell)
{
#if USEX11
    int j, k;
    int nbody = nj;
    int nc = nj-ni+1;
    double (*bpos)[3] = NULL;
    double (*cpos)[3] = NULL;
    double *csize = NULL;

    if (!bpos) {
	bpos = (double (*)[3])calloc(nbody, sizeof(double)*3);
	cpos = (double (*)[3])calloc(nc, sizeof(double)*3);
	csize = (double *)calloc(nc, sizeof(double));
    }
    nbody = 0;
    nc = 0;
    for (k = 0; k < 3; k++) {
	cpos[nc][k] = Cpos[c0][k];
	csize[nc] = Csize[c0];
    }
    nc++;
    for (j = 0; j < nj; j++) {
	int node = nodelist[j];
	if (iscell[j]) { /* node is an index to a cell */
	    for (k = 0; k < 3; k++) {
		cpos[nc][k] = Cpos[node][k];
	    }
	    csize[nc] = Csize[node];
	    nc++;
	}
	else {
	    for (k = 0; k < 3; k++) {
		bpos[nbody][k] = nb->x[node][k];
	    }
	    nbody++;
	}
    }
    vtc_viewlist(nb->n, nb->x, nbody, bpos, nc, cpos, csize, ni);
    free(bpos);
    free(cpos);
    free(csize);
#endif /* USEX11 */
}

#if 1
#define Sfprintf
#else
#define Sfprintf fprintf
#endif

static int
sanity_check(Cell c0, int level, Nbodyinfo *nb)
{
    int i, k, s;
    int err = 0;
    static int ntot = 0;

    if (level == 0) {
	ntot = 0;
    }
    ps(stderr, level);
    Sfprintf(stderr, "cell%d size: %f pos: %6.5f %6.5f %6.5f\n",
	    c0, Csize[c0], Cpos[c0][0], Cpos[c0][1], Cpos[c0][2]);
    ps(stderr, level);
    Sfprintf(stderr, "cmmass: %f cmpos: %6.5f %6.5f %6.5f ndescendant: %d\n",
	    Ccmmass[c0], Ccmpos[c0][0], Ccmpos[c0][1], Ccmpos[c0][2], Cndescendant[c0]);
    if (Cisleaf[c0]) {
	/* check if all particles are in the cell */
	Body * bp0 = Cbody[c0];
	int cn0 = Cnbody[c0];

	ps(stderr, level);
	Sfprintf(stderr, "is a leaf (Cnbody: %d)\n", cn0);
	for (i = 0; i < cn0; i++) {
	    Body b = Cbody[c0][i];
	    ps(stderr, level);
	    Sfprintf(stderr, "b: %d key: %016llx x: %6.5f %6.5f %6.5f",
		    b, Bkey[Cbody[c0]+i-Bindex], nb->x[b][0], nb->x[b][1], nb->x[b][2]);
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
	    Cell c1 = Cchild[c0][s];
	    if (NOCELL == c1) {
		continue;
	    }
	    if (Csize[c0]*0.5 != Csize[c1]) {
		fprintf(stderr, " NG (Csize)");
		err++;
	    }
	    for (k = 0; k < 3; k++) {
		if (fabs(Cpos[c0][k]-Cpos[c1][k]) != Csize[c1]*0.5) {
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
	fprintf(stderr, "root cell mass: %f pos: %f %f %f\n",
		Ccmmass[c0], Ccmpos[c0][0], Ccmpos[c0][1], Ccmpos[c0][2]);
    }
    return err;
}

static void
ps(FILE *fp, int len)
{
    int i;
    for (i = 0; i < len; i++) {
	Sfprintf(fp, "  ");
    }
}

static int
inbox(Cell c, Body b, Nbodyinfo *nb)
{
    int k;
    int isin = 1;
    double len = Csize[c]*0.5;
    double (*x)[3] = nb->x;

    for (k = 0; k < 3; k++) {
	Cfprintf(stderr, ">>> %6.5f < %6.5f < %6.5f ?\n",
		 Cpos[c][k]-len, x[b][k], Cpos[c][k]+len);
	if (Cpos[c][k]+len < x[b][k]) {
	    isin = 0;
	}
	if (Cpos[c][k]-len > x[b][k]) {
	    isin = 0;
	}
    }
    return (isin);
}
