
static void
m2m_1storder_md(Cell c0, Forceinfo *fi, Nbodyinfo *nb)
{
    int i, k, s;
    Cell c1;
    Body b;
    double m1;

    for (k = 0; k < 3; k++) {
	Ccmpos[c0][k] = 0.0;
	Ccmpos2[c0][k] = 0.0;
    }
    Ccmmass[c0] = 0.0;
    Ccmmass2[c0] = 0.0;
    Cndescendant[c0] = 0;
    Cnpositive[c0] = 0;
    Cnnegative[c0] = 0;
    if (Cisleaf[c0]) {
	for (i = 0; i < Cnbody[c0]; i++) {
	    b = Cbody[c0][i];
	    m1 = nb->m[b];
	    if (m1 > 0.0) {
		for (k = 0; k < 3; k++) {
		    Ccmpos[c0][k] += (nb->x[b][k])*m1;
		}
		Ccmmass[c0] += m1;
		Cnpositive[c0]++;
	    }
	    else {
		for (k = 0; k < 3; k++) {
		    Ccmpos2[c0][k] += (nb->x[b][k])*m1;
		}
		Ccmmass2[c0] += m1;
		Cnnegative[c0]++;
	    }
	}
	Cndescendant[c0] = Cnbody[c0];
    }
    else {
	for (s = 0; s < 8; s++) {
	    c1 = Cchild[c0][s];
	    if (NOCELL == c1) {
		continue;
	    }
	    m2m_1storder_md(c1, fi, nb);

	    /* for positive ions */
	    m1 = Ccmmass[c1];
	    for (k = 0; k < 3; k++) {
		Ccmpos[c0][k] += Ccmpos[c1][k]*m1;
	    }
	    Ccmmass[c0] += m1;

	    /* for negative ions */
	    m1 = Ccmmass2[c1];
	    for (k = 0; k < 3; k++) {
		Ccmpos2[c0][k] += Ccmpos2[c1][k]*m1;
	    }
	    Ccmmass2[c0] += m1;

	    Cndescendant[c0] += Cndescendant[c1];
	    Cnpositive[c0] += Cnpositive[c1];
	    Cnnegative[c0] += Cnnegative[c1];
	}
    }
    if (Ccmmass[c0] == 0.0) {
	for (k = 0; k < 3; k++) {
	    Ccmpos[c0][k] = 0.0;
	}
    }
    else {
	for (k = 0; k < 3; k++) {
	    Ccmpos[c0][k] /= Ccmmass[c0];
	}
    }
    if (Ccmmass2[c0] == 0.0) {
	for (k = 0; k < 3; k++) {
	  Ccmpos2[c0][k] = 0.0;
	}
    }
    else {
	for (k = 0; k < 3; k++) {
	    Ccmpos2[c0][k] /= Ccmmass2[c0];
	}
    }
}

static void
printmatrix(double a[3][3])
{
    int i, j;

    for (i = 0; i < 3; i++)
    {
	for (j = 0; j < 3; j++)
	{
	    fprintf(stderr, "%16.14f ", a[i][j]);
	}
	fprintf(stderr, "\n");
    }
}

#define SMALLVAL (1e-16)
static int /* 0: success -1: fail */
calc_pppos_2ndorder(double q[3][3], double cmpos[3], double cmmass,
		    double pppos[3][3], int c0)
{
    int i, j, k;
    int nrot;
    double val;
    double x[3][3];
    double eval[3];
    double p[3][3];

    /* q[][] --> diagonalized quad-tensor p[][] */
    if (0 > jacobi(q, eval, p, &nrot)) {
	fprintf(stderr, "c0: %d\n", c0);
	printmatrix(q);
	return(-1);
    }
    eigenvalsort(eval, p);

    /* calc pos in diagonalized coordinate system */
    x[0][0] = 0.0;
    val = (eval[0]+2.0*eval[1])/3.0/cmmass;
    if (val < 0.0) {
	if (-val < SMALLVAL) {
	    x[1][0] = 0.0;
	}
	else {
	    fprintf(stderr, "c0: %d x[1][0]^2 < 0.0: %e\n", c0, val);
	    fprintf(stderr, "eigen val: %e %e %e\n", eval[0], eval[1], eval[2]);
	    fprintf(stderr, "eigen vec:\n");
	    printmatrix(p);
	    return (-1);
	}
    }
    else {
	x[1][0] = 2.0*sqrt(val);
    }
    x[2][0] = 0.0;
    val = (2.0*eval[0]+eval[1])/cmmass;
    if (val < 0.0) {
	if (-val < SMALLVAL) {
	    x[0][1] = 0.0;
	}
	else {
	    fprintf(stderr, "c0: %d x[0][1]^2 < 0.0: %e\n", c0, val);
	    return (-1);
	}
    }
    else {
	x[0][1] = sqrt(val);
    }
    x[1][1] = -x[1][0]/2.0;
    x[2][1] = 0.0;
    x[0][2] = -x[0][1];
    x[1][2] = x[1][1];
    x[2][2] = 0.0;

    for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	    pppos[j][i] = 0.0;
	    for (k = 0; k < 3; k++) {
		pppos[j][i] += p[i][k] * x[k][j];
	    }
	}
    }
    for (k = 0; k < 3; k++) {
	for (i = 0; i < 3; i++) {
	    pppos[i][k] += cmpos[k];
	}
    }
    return (0);
}
#undef SMALLVAL

static void
err_2ndorder(Nbodyinfo *nb, int sign, int c0, double q[3][3],
	     double xx, double yy, double zz, double xy, double yz, double zx)
{
    int i;

    if (sign > 0) {
	fprintf(stderr, "positive cmmass: %e\n", Ccmmass[c0]);
    }
    else {
	fprintf(stderr, "negative cmmass: %e\n", Ccmmass2[c0]);
    }
    fprintf(stderr, "c0: %d isleaf: %d\n", c0, Cisleaf[c0]);
    printmatrix(q);
    fprintf(stderr, "xx: %e yy: %e zz: %e\n", xx, yy, zz);
    fprintf(stderr, "xy: %e yz: %e zx: %e\n", xy, yz, zx);
    if (Cisleaf[c0]) {
	fprintf(stderr, "nbody: %d\n", Cnbody[c0]);
	for (i = 0; i < Cnbody[c0]; i++) {
	    int b = Cbody[c0][i];
	    fprintf(stderr, "%d: pos: %e %e %e m: %e\n",
		    i, nb->x[b][0], nb->x[b][1], nb->x[b][2], nb->m[b]);
	}
    }
    else {
	fprintf(stderr, "ndescendant: %d\n", Cndescendant[c0]);
    }
}

static void
copy_pppos(Nbodyinfo *nb, Cell c0, int off, int nparticle, int sign,
	   int npp, /* # of pseudoparticles per cell */
	   int npe) /* # of pseudoparticles per expansion
		       (npp==2*npe if fi->negativemass, otherwise npp==npe)*/
{
    int i, k, s;
    int index = 0;
    double mchild;

    if (Cisleaf[c0]) {
	for (i = 0; i < Cnbody[c0]; i++) {
	    Body b = Cbody[c0][i];
	    mchild = nb->m[b];
	    if (sign * mchild >= 0.0) {
		Cppmass[npp*c0+index+off] = mchild;
		for (k = 0; k < 3; k++) {
		    Cpppos[npp*c0+index+off][k] = nb->x[b][k];
		}
		index++;
	    }
	}
    }
    else { /* non-leaf */
	for (s = 0; s < 8; s++) {
	    Cell c1 = Cchild[c0][s];
	    if (NOCELL == c1) {
		continue;
	    }
	    for (i = 0; i < npe; i++) {
		mchild = Cppmass[npp*c1+i+off];
		if (mchild == 0.0) {
		    continue;
		}
		Cppmass[npp*c0+index+off] = mchild;
		for (k = 0; k < 3; k++) {
		    Cpppos[npp*c0+index+off][k] = Cpppos[npp*c1+i+off][k];
		}
		index++;
	    }
	}
	if (index != nparticle) {
	    fprintf(stderr, "copy_pppos: index: %d != nparticle: %d", index, nparticle);
	    exit(1);
	}
    }
    for (i = index; i < npe; i++) {
	Cppmass[npp*c0+i+off] = 0.0;
	for (k = 0; k < 3; k++) {
	    Cpppos[npp*c0+i+off][k] = 0.0;
	}
    }
}

static void
m2m_2ndorder(Cell c0, Forceinfo *fi, Nbodyinfo *nb)
{
    int err;
    int i, j, k, s;
    Cell c1;
    int npp = fi->npp;
    double mchild, pchild[3];
    double xx, yy, zz, xy, yz, zx;
    double cmpos[3];
    double q[3][3];
    double pppos[3][3];

    /*
     * set Cppmass
     */
    Cppmass[npp*c0+0] = Cppmass[npp*c0+1] = Cppmass[npp*c0+2] = Ccmmass[c0]/3.0;

    /*
     * set Cpppos
     */
    /* calc quadrupole moment tensor q[][] at Ccmpos */
    xx = yy = zz = xy = yz = zx = 0.0;
    if (Cisleaf[c0]) {
	for (i = 0; i < Cnbody[c0]; i++) {
	    Body b = Cbody[c0][i];
	    mchild = nb->m[b];
	    for (k = 0; k < 3; k++) {
		pchild[k] = nb->x[b][k] - Ccmpos[c0][k];
	    }
	    xx += mchild*pchild[0]*pchild[0];
	    yy += mchild*pchild[1]*pchild[1];
	    zz += mchild*pchild[2]*pchild[2];
	    xy += mchild*pchild[0]*pchild[1];
	    yz += mchild*pchild[1]*pchild[2];
	    zx += mchild*pchild[2]*pchild[0];
	}
    }
    else {
	for (s = 0; s < 8; s++) {
	    c1 = Cchild[c0][s];
	    if (NOCELL == c1) {
		continue;
	    }
	    m2m_2ndorder(c1, fi, nb);
	    for (j = 0; j < npp; j++) {
		mchild = Cppmass[npp*c1+j];
		for (k = 0; k < 3; k++) {
		    pchild[k] = Cpppos[npp*c1+j][k] - Ccmpos[c0][k];
		}
		xx += mchild*pchild[0]*pchild[0];
		yy += mchild*pchild[1]*pchild[1];
		zz += mchild*pchild[2]*pchild[2];
		xy += mchild*pchild[0]*pchild[1];
		yz += mchild*pchild[1]*pchild[2];
		zx += mchild*pchild[2]*pchild[0];
	    }
	}
    }
    if (Ccmmass[c0] == 0.0) {
        for (k = 0; k < 3; k++) {
	    for (j = 0; j < 3; j++) {
	        Cpppos[npp*c0+j][k] = 0.0;
	    }
	}
    } else if (Cndescendant[c0] <= 3) { /* do not apply quadpole expansion.
					   (c0 contains not more than 3 particles) */
	copy_pppos(nb, c0, 0, Cndescendant[c0], +1, npp, npp);
    }
    else { /* apply quadpole expansion */
	q[0][0] = xx-0.5*yy-0.5*zz;
	q[1][1] = yy-0.5*zz-0.5*xx;
	q[2][2] = zz-0.5*xx-0.5*yy;
	q[0][1] = q[1][0] = 1.5 * xy;
	q[1][2] = q[2][1] = 1.5 * yz;
	q[2][0] = q[0][2] = 1.5 * zx;

	for (k = 0; k < 3; k++) {
	    cmpos[k] = Ccmpos[c0][k];
	}
	err = calc_pppos_2ndorder(q, cmpos, Ccmmass[c0], pppos, c0);
	if (err) {
	    err_2ndorder(nb, 1, c0, q, xx, yy, zz, xy, yz, zx);
	    exit(1);
	}
	for (k = 0; k < 3; k++) {
	    for (j = 0; j < 3; j++) {
		Cpppos[npp*c0+j][k] = pppos[j][k];
	    }
	}
    }

#if 0 /* !!! do not activate for normal operation */
    Cppmass[npp*c0+0] = Ccmmass[c0];
    Cppmass[npp*c0+1] = Cppmass[npp*c0+2] = 0.0;

    Cppmass[npp*c0+0] = Cppmass[npp*c0+1] = Cppmass[npp*c0+2] = Ccmmass[c0]/3.0;

    for (k = 0; k < 3; k++) {
	for (j = 0; j < 3; j++) {
	    Cpppos[npp*c0+j][k] = Ccmpos[c0][k];
	}
    }
#endif
 

}

static void
m2m_2ndorder_md(Cell c0, Forceinfo *fi, Nbodyinfo *nb)
{
    int err;
    int i, j, k, s;
    Cell c1;
    int npp = fi->npp;
    double mchild, pchild[3];
    double xx, yy, zz, xy, yz, zx;
    double xx2, yy2, zz2, xy2, yz2, zx2;
    double cmpos[3];
    double q[3][3];
    double pppos[3][3];

    /*
     * set Cppmass
     */
    Cppmass[npp*c0+0] = Cppmass[npp*c0+1] = Cppmass[npp*c0+2] = Ccmmass[c0]/3.0;
    Cppmass[npp*c0+3] = Cppmass[npp*c0+4] = Cppmass[npp*c0+5] = Ccmmass2[c0]/3.0;

    /*
     * set Cpppos
     */
    /* calc quadrupole moment tensor q[][] at Ccmpos */
    xx = yy = zz = xy = yz = zx = 0.0;
    xx2 = yy2 = zz2 = xy2 = yz2 = zx2 = 0.0;
    if (Cisleaf[c0]) {
	for (i = 0; i < Cnbody[c0]; i++) {
	    Body b = Cbody[c0][i];
	    mchild = nb->m[b];
	    if (mchild >= 0.0) {
		for (k = 0; k < 3; k++) {
		    pchild[k] = nb->x[b][k] - Ccmpos[c0][k];
		}
		xx += mchild*pchild[0]*pchild[0];
		yy += mchild*pchild[1]*pchild[1];
		zz += mchild*pchild[2]*pchild[2];
		xy += mchild*pchild[0]*pchild[1];
		yz += mchild*pchild[1]*pchild[2];
		zx += mchild*pchild[2]*pchild[0];
	    }
	    else if (mchild < 0.0) {
		for (k = 0; k < 3; k++) {
		    pchild[k] = nb->x[b][k] - Ccmpos2[c0][k];
		}
		xx2 -= mchild*pchild[0]*pchild[0];
		yy2 -= mchild*pchild[1]*pchild[1];
		zz2 -= mchild*pchild[2]*pchild[2];
		xy2 -= mchild*pchild[0]*pchild[1];
		yz2 -= mchild*pchild[1]*pchild[2];
		zx2 -= mchild*pchild[2]*pchild[0];
	    }
	}
    }
    else { /* non-leaf cell */
	for (s = 0; s < 8; s++) {
	    c1 = Cchild[c0][s];
	    if (NOCELL == c1) {
		continue;
	    }
	    m2m_2ndorder_md(c1, fi, nb);
	    if (Cnpositive[c0] > 3) {
		for (j = 0; j < 3; j++) {
		    mchild = Cppmass[npp*c1+j];
		    for (k = 0; k < 3; k++) {
			pchild[k] = Cpppos[npp*c1+j][k] - Ccmpos[c0][k];
		    }
		    xx += mchild*pchild[0]*pchild[0];
		    yy += mchild*pchild[1]*pchild[1];
		    zz += mchild*pchild[2]*pchild[2];
		    xy += mchild*pchild[0]*pchild[1];
		    yz += mchild*pchild[1]*pchild[2];
		    zx += mchild*pchild[2]*pchild[0];
		}
	    }
	    if (Cnnegative[c0] > 3) {
		for (j = 0; j < 3; j++) {
		    mchild = Cppmass[npp*c1+j+3];
		    for (k = 0; k < 3; k++) {
			pchild[k] = Cpppos[npp*c1+j+3][k] - Ccmpos2[c0][k];
		    }
		    xx2 -= mchild*pchild[0]*pchild[0];
		    yy2 -= mchild*pchild[1]*pchild[1];
		    zz2 -= mchild*pchild[2]*pchild[2];
		    xy2 -= mchild*pchild[0]*pchild[1];
		    yz2 -= mchild*pchild[1]*pchild[2];
		    zx2 -= mchild*pchild[2]*pchild[0];
		}
	    }
	}
    }

    if (Ccmmass[c0] == 0.0) {
	for (k = 0; k < 3; k++) {
	    for (j = 0; j < 3; j++) {
		Cpppos[npp*c0+j][k] = 0.0;
	    }
	}
    }
    else if (Cnpositive[c0] <= 3) { /* do not apply quadpole expansion
				  (c0 contains not more than 3 positive-mass
				  particles) */
	copy_pppos(nb, c0, 0, Cnpositive[c0], +1, npp, npp/2);
    }
    else { /* process positive-mass particles */
	q[0][0] = xx-0.5*yy-0.5*zz;
	q[1][1] = yy-0.5*zz-0.5*xx;
	q[2][2] = zz-0.5*xx-0.5*yy;
	q[0][1] = q[1][0] = 1.5 * xy;
	q[1][2] = q[2][1] = 1.5 * yz;
	q[2][0] = q[0][2] = 1.5 * zx;
	for (k = 0; k < 3; k++) {
	    cmpos[k] = Ccmpos[c0][k];
	}
	err = calc_pppos_2ndorder(q, cmpos, Ccmmass[c0], pppos, c0);
	if (err) {
	    err_2ndorder(nb, 1, c0, q, xx, yy, zz, xy, yz, zx);
	    exit(1);
	}
	for (k = 0; k < 3; k++) {
	    for (j = 0; j < 3; j++) {
		Cpppos[npp*c0+j][k] = pppos[j][k];
	    }
	}
    }

    if (Ccmmass2[c0] == 0.0) {
	for (k = 0; k < 3; k++) {
	    for (j = 0; j < 3; j++) {
		Cpppos[npp*c0+j+3][k] = 0.0;
	    }
	}
    }
    else if (Cnnegative[c0] <= 3) { /* do not apply quadpole expansion
				  (c0 contains not more than 3 negative-mass
				  particles) */
	copy_pppos(nb, c0, 3, Cnnegative[c0], -1, npp, npp/2);
    }
    else { /* process negative-mass particles */
	q[0][0] = (xx2-0.5*yy2-0.5*zz2);
	q[1][1] = (yy2-0.5*zz2-0.5*xx2);
	q[2][2] = (zz2-0.5*xx2-0.5*yy2);
	q[0][1] = q[1][0] = (1.5 * xy2);
	q[1][2] = q[2][1] = (1.5 * yz2);
	q[2][0] = q[0][2] = (1.5 * zx2);
	for (k = 0; k < 3; k++) {
	    cmpos[k] = Ccmpos2[c0][k];
	}
	err = calc_pppos_2ndorder(q, cmpos, -Ccmmass2[c0], pppos, c0);
	if (err) {
	    err_2ndorder(nb, -1, c0, q, xx2, yy2, zz2, xy2, yz2, zx2);
	    exit(1);
	}
	for (k = 0; k < 3; k++) {
	    for (j = 0; j < 3; j++) {
		Cpppos[npp*c0+j+3][k] = pppos[j][k];
	    }
	}
    }
}


/* 
 * form p2m2 expansion at c0 from the one at c1.
 * if c1 == NOCELL, form p2m2 expansion from physical particles inside c0.
 */
static void
expand_by_pp(Forceinfo *fi, Nbodyinfo *nb, Cell c0, Cell c1, double *dm)
{
    int i, j, k, q;
    int npp = fi->npp;
    double rr, cosg;
    double r1, m1;
    double x1[3]; /* physical particle pos */
    static double x0[NPPMAX][3]; /* pseudoparticle pos */
    static double r0[NPPMAX];
    static double pln[TDESIGNMAX/2+1];

    for (j = 0; j < npp; j++) {
	dm[j] = 0.0;
	r0[j] = 0.0;
	for (k = 0; k < 3; k++) {
	    x0[j][k] = fi->pppos[j][k] * Csize[c0];
	    r0[j] += x0[j][k]*x0[j][k];
	}
	r0[j] = sqrt(r0[j]);
    }
    if (NOCELL == c1) { /* form expansion from physical particles inside c0 */
	for (i = 0; i < Cnbody[c0]; i++) {
	    Body b = Cbody[c0][i];
	    m1 = nb->m[b];
	    r1 = 0.0;
	    for (k = 0; k < 3; k++) {
		x1[k] = nb->x[b][k] - Ccmpos[c0][k];
		r1 += x1[k]*x1[k];
	    }
	    r1 = sqrt(r1);
	    for (j = 0; j < npp; j++) {
		rr = r1/r0[j];
		cosg = (x0[j][0]*x1[0]+x0[j][1]*x1[1]+x0[j][2]*x1[2])/r0[j]/r1;
		plgndr0(fi->tdesign/2+1, cosg, pln);
		for (q = 0; q < fi->tdesign/2+1; q++) {
		    dm[j] = dm[j] + pow(rr, (double)q)*m1*(2*q+1)*pln[q];
		}
	    }
	}
    }
    else { /* form expansion from the one at c1 */
	for (i = 0; i < npp; i++) {
	    m1 = Cppmass[npp*c1+i];
	    r1 = 0.0;
	    for (k = 0; k < 3; k++) {
		x1[k] = Cpppos[npp*c1+i][k] - Ccmpos[c0][k];
		r1 += x1[k]*x1[k];
	    }
	    r1 = sqrt(r1);
	    for (j = 0; j < npp; j++) {
		rr = r1/r0[j];
		cosg = (x0[j][0]*x1[0]+x0[j][1]*x1[1]+x0[j][2]*x1[2])/r0[j]/r1;
		plgndr0(fi->tdesign/2+1, cosg, pln);
		for (q = 0; q < fi->tdesign/2+1; q++) {
		    dm[j] = dm[j] + pow(rr, (double)q)*m1*(2*q+1)*pln[q];
		}
	    }
	}
    }
    for (j = 0; j < npp; j++) {
	dm[j] /= npp;
    }
}

static void
m2m_anyorder(Cell c0, Forceinfo *fi, Nbodyinfo *nb)
{
    Cell c1;
    int i, k, s;
    int npp = fi->npp;
    double dm[NPPMAX];

    for (s = 0; s < 8; s++) {
	c1 = Cchild[c0][s];
	if (NOCELL == c1) {
	    continue;
	}
	m2m_anyorder(c1, fi, nb);
    }

    if (Cndescendant[c0] <= npp) { /* do not apply p2m2 expansion.
				    (c0 contains not more than npp particles) */
	copy_pppos(nb, c0, 0, Cndescendant[c0], +1, npp, npp);
    }
    else { /* apply p2m2 expansion */
	for (i = 0; i < npp; i++) {
	    for (k = 0; k < 3; k++) {
		Cpppos[npp*c0+i][k] = Ccmpos[c0][k] + fi->pppos[i][k] * Csize[c0];
	    }
	    Cppmass[npp*c0+i] = 0.0;
	}
	if (Cisleaf[c0]) {
	    expand_by_pp(fi, nb, c0, NOCELL, dm);
	    for (i = 0; i < npp; i++) {
		Cppmass[npp*c0+i] += dm[i];
	    }
	}
	else {
	    for (s = 0; s < 8; s++) {
		c1 = Cchild[c0][s];
		if (NOCELL == c1) {
		    continue;
		}
		expand_by_pp(fi, nb, c0, c1, dm);
		for (i = 0; i < npp; i++) {
		    Cppmass[npp*c0+i] += dm[i];
		}
	    }
	}
    }
}

static void
m2m_anyorder_md(Cell c0, Forceinfo *fi, Nbodyinfo *nb)
{
    fprintf(stderr, "m2m_anyorder_md not implemented yet\n");
    fprintf(stderr, "requested expansion order: %d\n", fi->p);
    exit(1);
}
