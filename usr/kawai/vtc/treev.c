/*
 * vectorized functions
 * in Makefile, set VECRORIZED to 1 to activate this part
 *
 * this may be faster on some platforms
 */

static int /* # of nodes packed */
pack_interaction_list(Nbodyinfo *nb, int jlen, int npp, int usenegativemass,
		      Node *nodelist, int *iscell, double *mlist, double (*xlist)[3])
{
    int i, j;
    int pplen = 0;
    Node node;
    double *nbm, *nbx0, *nbx1, *nbx2;
    double *cpos0, *cpos1, *cpos2;
    double *cpos20, *cpos21, *cpos22;

    nbm = nb->m;
    nbx0 = (double *)&(nb->x[0][0]);
    nbx1 = (double *)&(nb->x[0][1]);
    nbx2 = (double *)&(nb->x[0][2]);

    if (usep2m2) {
	cpos0 = (double *)&(Cpppos[0][0]);
	cpos1 = (double *)&(Cpppos[0][1]);
	cpos2 = (double *)&(Cpppos[0][2]);
	for (j = 0; j < jlen; j++) {
	    node = nodelist[j];
	    if (!iscell[j]) { /* node is an index to a particle */
		mlist[pplen] = nbm[node];
		xlist[pplen][0] = nbx0[node*3];
		xlist[pplen][1] = nbx1[node*3];
		xlist[pplen][2] = nbx2[node*3];
		pplen++;
	    }
	}
#if 1
	for (j = 0; j < jlen; j++) {
	    node = nodelist[j];
	    if (iscell[j]) { /* node is an index to a cell */
		for (i = 0; i < npp; i++) {
		    mlist[pplen] = Cppmass[npp*node+i];
		    xlist[pplen][0] = cpos0[(npp*node+i)*3];
		    xlist[pplen][1] = cpos1[(npp*node+i)*3];
		    xlist[pplen][2] = cpos2[(npp*node+i)*3];
		    pplen++;
		}
	    }
	}
#else /* this may be faster */
	for (i = 0; i < npp; i++) {
	    for (j = 0; j < jlen; j++) {
		if (iscell[j]) { /* node is an index to a cell */
		    node = nodelist[j];
		    mlist[pplen] = Cppmass[npp*node+i];
		    xlist[pplen][0] = cpos0[(npp*node+i)*3];
		    xlist[pplen][1] = cpos1[(npp*node+i)*3];
		    xlist[pplen][2] = cpos2[(npp*node+i)*3];
		    pplen++;
		}
	    }
	}
#endif
    }
    else if (usenegativemass) { /* no P2M2 */
	cpos0 = (double *)&(Ccmpos[0][0]);
	cpos1 = (double *)&(Ccmpos[0][1]);
	cpos2 = (double *)&(Ccmpos[0][2]);
	cpos20 = (double *)&(Ccmpos2[0][0]);
	cpos21 = (double *)&(Ccmpos2[0][1]);
	cpos22 = (double *)&(Ccmpos2[0][2]);
	for (j = 0; j < jlen; j++) {
	    node = nodelist[j];
	    if (!iscell[j]) { /* node is an index to a particle */
		mlist[pplen] = nbm[node];
		xlist[pplen][0] = nbx0[node*3];
		xlist[pplen][1] = nbx1[node*3];
		xlist[pplen][2] = nbx2[node*3];
		pplen++;
	    }
	}
	for (j = 0; j < jlen; j++) {
	    node = nodelist[j];
	    if (iscell[j] && Ccmmass[node] != 0.0) { /* node is an index to a cell */
		mlist[pplen] = Ccmmass[node];
		xlist[pplen][0] = cpos0[node*3];
		xlist[pplen][1] = cpos1[node*3];
		xlist[pplen][2] = cpos2[node*3];
		pplen++;
	    }
	}
	for (j = 0; j < jlen; j++) {
	    node = nodelist[j];
	    if (iscell[j] && Ccmmass2[node] != 0.0) { /* node is an index to a cell */
		mlist[pplen] = Ccmmass2[node];
		xlist[pplen][0] = cpos20[node*3];
		xlist[pplen][1] = cpos21[node*3];
		xlist[pplen][2] = cpos22[node*3];
		pplen++;
	    }
	}
    }
    else { /* no P2M2, no negative mass particle */
	cpos0 = (double *)&(Ccmpos[0][0]);
	cpos1 = (double *)&(Ccmpos[0][1]);
	cpos2 = (double *)&(Ccmpos[0][2]);

	/* note that particles must be packed first since in traverse_tree()
	   the first ilen nodes are assumed to be i-particles */
	for (j = 0; j < jlen; j++) {
	    node = nodelist[j];
	    if (!iscell[j]) { /* node is an index to a particle */
		mlist[pplen] = nbm[node];
		xlist[pplen][0] = nbx0[node*3];
		xlist[pplen][1] = nbx1[node*3];
		xlist[pplen][2] = nbx2[node*3];
		pplen++;
	    }
	}
	for (j = 0; j < jlen; j++) {
	    node = nodelist[j];
	    if (iscell[j]) { /* node is an index to a cell */
		mlist[pplen] = Ccmmass[node];
		xlist[pplen][0] = cpos0[node*3];
		xlist[pplen][1] = cpos1[node*3];
		xlist[pplen][2] = cpos2[node*3];
		pplen++;
	    }
	}
    }
    return (pplen);
}

static int /* length of the list excluding i-particles */
get_interaction_list(Cell icell, Cell root, double theta2, Nbodyinfo *nb,
		     Node *list, int *iscell, int npp)
{
    int i, j, k, s;
    Cell jcell0, jcell1;
    int nlist;
    int ncelltmp; /* # of cells in the level in question */
    int ncelltmp_new; /* # of cells to be processed in the next level */
    double isize, ipos[3];
    double xmin;
#if DYNMEM
#if !BFLOAD    
    static Cell *celltmp = NULL;
    static Cell *celltmp_new = NULL;
#endif /* BFLOAD */
    static int *qualified = NULL;
    static double *dr2 = NULL;
#else /* !DYNMEM */
#if !BFLOAD    
    static Cell celltmp[NCELLMAX];
    static Cell celltmp_new[NCELLMAX];
#endif
    static int qualified[NCELLMAX];
    static double dr2[NCELLMAX];
#endif /* DYNMEM */

#if DYNMEM
    if (qualified == NULL) {
#if !BFLOAD    
	celltmp = (Cell *)realloc(celltmp, sizeof(Cell)*nb->n);
	celltmp_new = (Cell *)realloc(celltmp_new, sizeof(Cell)*nb->n);
#endif /* BFLOAD */
	qualified = (int *)realloc(qualified, sizeof(int)*nb->n);
	dr2 = (double *)realloc(dr2, sizeof(double)*nb->n);
    }
#endif /* DYNMEM */

    isize = Csize[icell];
    xmin = (isize)*0.5;
    for (k = 0; k < 3; k++) {
	ipos[k] = Cpos[icell][k];
    }
    nlist = 0;
    ncelltmp = ncelltmp_new = 1;
    celltmp[0] = root;

    while (ncelltmp > 0) {
	double *cp0, *cp1, *cp2;
	ncelltmp = ncelltmp_new;
	ncelltmp_new = 0;

	cp0 = (double *)&(Cpos[0][0]);
	cp1 = (double *)&(Cpos[0][1]);
	cp2 = (double *)&(Cpos[0][2]);

	for (i = 0; i < ncelltmp; i++) {
	    double dx, r;

	    r = 0.0;
	    dx = ipos[0]-cp0[celltmp[i]*3];
	    dx = dx > 0.0 ? dx:-dx;
	    dx -= xmin;
	    dx = dx > 0.0 ? dx:0.0;
	    r += dx*dx;

	    dx = ipos[1]-cp1[celltmp[i]*3];
	    dx = dx > 0.0 ? dx:-dx;
	    dx -= xmin;
	    dx = dx > 0.0 ? dx:0.0;
	    r += dx*dx;

	    dx = ipos[2]-cp2[celltmp[i]*3];
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
	    if ((qualified[i] == FALSE) &&
		(Cisleaf[celltmp[i]] || Cndescendant[celltmp[i]] <= npp)) {
	      qualified[i] = TRUE;
	    }
	}
	/* leaf cell. handle particle force directly */
	for (i = 0; i < ncelltmp; i++) {
	    if (qualified[i] == TRUE) {
		Body *cb = Cbody[celltmp[i]];
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
		    jcell1 = Cchild[celltmp[i]][s];
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
