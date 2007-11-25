/****************************************************************************/
/* ANALYSIS.C: routines to analize calculated data                          */
/* Copyright (c) 2000 by Jin Koda, Tokyo, JAPAN.                            */
/****************************************************************************/

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "/usr/local/pgplot/cpgplot.h"

/*
 * local routines and variables
 */


void getdata(char *);
void transform(int, real *, real *, real *, real *, real *, real, real);
void calcmesh(int , real *, real *, real *, int, int, real *, real, real);
void calcmesh2(int, real *, real *, real *, real *, int, int, real *, real, real);
void plot(real *, int, int, real, real, int, real *, int);
void palett(int, real, real);
void setvp(void);
void fiddle(real, real, char);
void explain(void);

#define NMESH1 128
#define NMESH2 NMESH1

int nbody, ndim;
real tnow, *m, *x, *y, *vx, *vy;
real *xt, *yt, *zt, *vxt, *vyt, *vzt;
real mscale, rscale, tscale, pscale;

/*
 * MAIN: main function.
 */

int main(int argc, string argv[])
{

#define NALEV 21
    
    int pal, ii, jj, cnt;
    real theta, psi, omgb;
    real meshd[NMESH2][NMESH1], xs, yl, r, l, alev[NALEV];
    char fname[32], ch;


    if (argc < 2) {
        printf("No output file name\n"); 
        printf("Usage:                                 \n"); 
        printf("   # analysis (file) (pattarn speed)   \n"); 
        printf("  (ex.) # analysis data 30.0           \n"); 
    }

    /*                                       */
    /* Set input file name and get body data */
    /*                                       */

    strcpy(fname, argv[1]);
    printf("\n\tInput file = %s\n\n", fname);
    getdata(fname);

    /*                   */
    /* Set pattarn speed */
    /*                   */ 

    if (argc < 3) {
        omgb = 0.0;
        printf("No pattarn speed is set.\n");
    } else {
        sscanf(argv[2], "%f", &omgb);
        printf("Pattarn Speed  = %.2f [km/s/kpc]\n", omgb);
    }
    omgb = 1.02271202e-9 * omgb;

    /*                      */
    /* Set bar horizontally */
    /*                      */

    theta = omgb * tnow;
    transform(nbody, x, y, x, y, zt, -theta, 0.0);
    transform(nbody, vx, vy, vx, vy, vzt, -theta,  0.0);

    /*               */
    /* Open pgdevice */
    /*               */

    cpgbeg(0, "/xwindow", 1, 1);

    /*                                    */
    /* Initialize code and explain usages */
    /*                                    */

    pscale = 3.0 * rscale;
    ch = 'S';
    xs = 0.5;
    yl = 1.0;
    theta = 0.0;
    psi = 0.0;
    explain();

    /*               */
    /* Show graphics */
    /*               */

    while (1) {

	if (ch == 'X') {          /* break interaction */
	    break;
	} else if (ch == 'Z') {   /* zoom box          */
	    pscale /= 1.5;
	} else if (ch == 'E') {   /* expand box        */
	    pscale *= 1.5;
	} else if (ch == 'H') {   /* show help message */
	    explain();
	} else if (ch == 'S') {   /* draw snap shot    */
	    calcmesh(nbody,x,y,m,NMESH1,NMESH2,(real *) meshd,pscale,pscale);
	    plot((real *) meshd,NMESH1,NMESH2,pscale,pscale,0,alev,2);
	} else if (ch == 'P') {   /* draw PV diagram   */
	    transform(nbody, x, y, xt, yt, zt, -theta, 0.0);
	    transform(nbody, vx, vy, vxt, vyt, vzt, -theta, 0.0);
	    calcmesh(nbody,yt,vxt,m,NMESH1,NMESH2,(real *) meshd,pscale,800.0);
	    plot((real *) meshd, NMESH1, NMESH2, pscale, 800.0, 0, alev, 2);
	} else if (ch == 'V') {   /* draw vel. field   */
	    transform(nbody, x, y, xt, yt, zt, -theta, psi);
	    transform(nbody, vx, vy, vxt, vyt, vzt, -theta, psi);
	    calcmesh2(nbody,xt,yt,vzt,m,NMESH1,NMESH2,(real *) meshd,pscale,pscale);
	    for (ii=0; ii<NALEV; ii++)
		alev[ii] = 30.0 * (real) (ii-10.0);
	    plot((real *) meshd, NMESH1, NMESH2, pscale, pscale, NALEV, alev, 3);
	} else {                  /* fiddle            */
	    fiddle(xs, yl, ch);
	}
	
	if ( cpgcurs(&xs, &yl, &ch) != 1) {
	    printf("ERROR: cpgcurs\n");
	    break;
	}

	ch = toupper(ch);
	r = sqrt(xs*xs + yl*yl);
	theta = acos( xs / r );
	if (yl < 0.0) theta = 2.0 * PI - theta;
	l = pscale / 4.0;
	psi = acos( l / sqrt(r*r + l*l) );
	printf("theta = %.2f;  incli = %.2f\n",  (theta/PI*180.0), (psi/PI*180.0) );
	
    }
    
    /*            */
    /* End pgplot */
    /*            */

    cpgend();
}


/*
 * TRANSFORM: rotate coordinate. (xx, yy) -> (xxt, yyt, zzt)
 */

void transform(int nn, real *xx, real *yy,
	       real *xxt, real *yyt, real *zzt, real theta, real psi)
{
    int i;
    real xx0, yy0, zz0;

    for(i=0; i<nn; i++) {
	xx0 = xx[i];
	yy0 = yy[i];
	xxt[i] = cos(theta) * xx0 - sin(theta) * yy0;
	yyt[i] = sin(theta) * xx0 + cos(theta) * yy0;
    }
    if (psi > 0.0) {
	for (i=0; i<nn; i++) {
	    zz0 = 0.0;
	    xx0 = xxt[i];
	    zzt[i] = cos(-psi) * zz0 - sin(-psi) * xx0;
	    xxt[i] = sin(-psi) * zz0 + cos(-psi) * xx0;
	}
	for(i=0; i<nn; i++) {
	    xx0 = xxt[i];
	    yy0 = yyt[i];
	    xxt[i] = cos(-theta) * xx0 - sin(-theta) * yy0;
	    yyt[i] = sin(-theta) * xx0 + cos(-theta) * yy0;
	}
    }
}

/*
 * CALCMESH: Assign particle to mesh with Cloud-In-Cell (CIC) method.
 */

void calcmesh(int nn, real *xx, real *yy, real *weight,
	      int ni, int nj, real *f, real xsize, real ysize)
{

    int i, j, p;
    real xx0, yy0, ddx, ddy, ttx, tty, dx, dy, rsq;

    dx = xsize / ((real) NMESH1);
    dy = ysize / ((real) NMESH2);

    for (i=0; i<NMESH1; i++)
	for (j=0; j<NMESH2; j++)
	    f[ni*i+j] = 0.0;

    for (p=0; p<nn; p++) {
	xx0 = (xx[p] + 0.5 * xsize - 0.5 * dx)/dx;
	yy0 = (yy[p] + 0.5 * ysize - 0.5 * dy)/dy;
	i = (int) xx0;
	j = (int) yy0;
	if (i >= 0 && i < NMESH1-1 && j >= 0 && j < NMESH2-1) {
	    ddx = xx0 - (real) i;
	    ddy = yy0 - (real) j;
	    ttx = 1.0 - ddx;
	    tty = 1.0 - ddy;
	    f[ni*(j  ) + (i  )] += ttx*tty*weight[p];
	    f[ni*(j  ) + (i+1)] += ddx*tty*weight[p];
	    f[ni*(j+1) + (i  )] += ttx*ddy*weight[p];
	    f[ni*(j+1) + (i+1)] += ddy*ddy*weight[p];

	}
    }
}


/*
 * CALCMESH2: Assign particle to mesh, and normalize it
 * with secondary weighting.
 */

void calcmesh2(int nn, real *xx, real *yy, real *weight1, real *weight2,
	       int ni, int nj, real *f, real xsize, real ysize)
{

    int i, j;
    real *w;

    w = malloc(ni * nj * sizeof(real));

    calcmesh(nn, xx, yy, weight1, ni, nj, (real *) f, pscale, pscale);
    calcmesh(nn, xx, yy, weight2, ni, nj, (real *) w, pscale, pscale);

    for (i=0; i<ni; i++) {
	for (j=0; j<nj; j++) {
	    if (w[ni*j+i] > 0.0 * weight2[0]) {
		f[ni*j+i] /= (w[ni*j+i] / weight2[0]);
	    } else {
		f[ni*j+i] = 1.0e9;
	    }
	}
    }
    free(w);
}


/*
 * PLOT: plot data
 */

void plot(real *f, int ni, int nj, real xsize, real ysize, int nalev, real *alev, int pal)
{
 
    int   ii, jj, ll, c1, c2, nc, nij;
    real fmin, fmax, tr[6], contra, bright;
    char  val[32], s[32];

    /*                  */
    /*  Set color index */
    /*                  */

    cpgqcir(&c1, &c2);                /* inquire color index range  */
    nc = c2-c1+1;                     /* define the number of color */


    /*                                  */
    /*  Maximum and minimum in array F. */ 
    /*                                  */

    nij  = ni  * nj;
    fmin =  1.0e9;
    fmax = -1.0e9;
    for (ii=0; ii<nij; ii++) {
	if (*(f+ii) > -1.0e9) fmin = (*(f+ii) < fmin) ? *(f+ii) : fmin;
	if (*(f+ii) <  1.0e9) fmax = (*(f+ii) > fmax) ? *(f+ii) : fmax;
    }

    /**********************************/
    /*  Simple transformation matrix  */
    /**********************************/

    /*                                           */
    /* Set the coordinate transformation matrix: */
    /*                                           */

    tr[0] = 0.0;
    tr[1] = 1.0;
    tr[2] = 0.0;
    tr[3] = 0.0;
    tr[4] = 0.0;
    tr[5] = 1.0;

    /*                                               */
    /* Clear the screen. Set up window and viewport. */
    /*                                               */

    cpgpap(6.0, 1.3);
    cpgpage();
    setvp();
    cpgwnad(0.0, 1.0+ni, 0.0, 1.0+nj);

    /*                       */
    /* Set up the color map. */
    /*                       */

    bright = 0.5;
    contra = 1.0;
    palett(pal, contra, bright);
    
    /*                           */
    /* Draw the map with PGIMAG. */
    /*                           */
    
    cpgimag(f,ni,nj,1,ni,1,nj,fmin,fmax,tr);

    /*               */
    /* Draw contours */
    /*               */

    for (ii=0; ii<nalev; ii++)
	cpgcont(f, ni, nj, 1, ni, 1, nj, &alev[ii], -1, tr);

    /*                    */
    /* Annotate the plot. */
    /*                    */

    cpgsch(0.6);
    cpgswin(-xsize/2.0, xsize/2.0, -ysize/2.0, ysize/2.0);
    cpgsch(1.0);
    cpgbox("bcntsi",0.0,0,"bcntsi",0.0,0);

    /*               */
    /* Draw a wedge. */
    /*               */

    cpgwedg("bi", 3.0, 3.0,fmin,fmax,"pixel value");
    cpgsch(1.0); 

    cpgask(0);
}


/*
 * PALETT: set up a color table.
 */

void palett(int type, real contra, real bright)
{
    /*                      */
    /* Define color palett. */
    /*                      */
    const real gl[] = {0.0, 1.0};
    const real gr[] = {0.0, 1.0};
    const real gg[] = {0.0, 1.0};
    const real gb[] = {0.0, 1.0};

    const real rl[] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
    const real rr[] = { 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
    const real rg[] = { 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
    const real rb[] = { 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};

    const real ril[] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
    const real rir[] = { 1.0, 1.0,  1.0,  1.0,  0.6,  0.0,  0.0, 0.0, 0.0};
    const real rig[] = { 1.0, 0.0,  0.6,  1.0,  1.0,  1.0,  0.0, 0.0, 0.0};
    const real rib[] = { 1.0, 0.0,  0.0,  0.0,  0.3,  1.0,  0.8, 0.3, 0.0};

    const real hl[] = {0.0, 0.2, 0.4, 0.6, 1.0};
    const real hr[] = {0.0, 0.5, 1.0, 1.0, 1.0};
    const real hg[] = {0.0, 0.0, 0.5, 1.0, 1.0};
    const real hb[] = {0.0, 0.0, 0.0, 0.3, 1.0};
    
    const real wl[] = {0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0};
    const real wr[] = {0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0};
    const real wg[] = {0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0};
    const real wb[] = {0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0};
    
    const real al[] ={0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0};
    const real ar[] = {0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    const real ag[] = {0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0};
    const real ab[] = {0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    /*                                                */
    /* Install the color table to be used by cpgimag. */
    /*                                                */
    switch (type) {
    case 1: /* gray scale */
	cpgctab(gl, gr, gg, gb,  2, contra, bright);
	break;
    case 2: /* rainbow */
	cpgctab(rl, rr, rg, rb,  9, contra, bright);
	break;
    case 3: /* rainbow inverse */
	cpgctab(ril, rir, rig, rib,  9, contra, bright);
	break;
    case 4: /* heat */
	cpgctab(hl, hr, hg, hb,  5, contra, bright);
	break;
    case 5: /* weird IRAF */
	cpgctab(wl, wr, wg, wb, 10, contra, bright);
	break;
    case 6: /* AIPS */
	cpgctab(al, ar, ag, ab, 10, contra, bright);
	break;
    }
}

/*
 * SETVP: set up viewport.
 */

void setvp(void)
{
    real d, dx, dy, vpx1, vpx2, vpy1, vpy2;

    cpgsvp(0.0,1.0,0.0,1.0);
    cpgqvp(1, &vpx1, &vpx2, &vpy1, &vpy2);

    dx = vpx2-vpx1;
    dy = vpy2-vpy1;
    d  = (dx < dy) ? dx/40.0 : dy/40.0;    /* min(dx,dy)/40.0 */

    vpx1 = vpx1 + 2.0 * d;
    vpx2 = vpx2;
    vpy1 = vpy1;
    vpy2 = vpy2;

    cpgvsiz(vpx1, vpx2, vpy1, vpy2);

}



/***************************/
/* Fiddle the color table. */
/***************************/

void fiddle(real xs, real yl, char ch)
{
    static int pal=2;
    float contra, bright, sign;
    float x1, x2, y1, y2, b1, b2, c1, c2;
    
    contra = 1.0;
    bright = 0.5;
    sign = +1.0;

    cpgqwin(&x1, &x2, &y1, &y2);
    b1 = 0.0;
    b2 = 1.0;
    c1 = 0.0;
    c2 = 10.0;
    cpgswin(b1,b2,c1,c2);

    /* Change contrast and brightness */
    if ( ch=='A' ) {
	bright = (b2 < xs) ? b2 : xs;
	bright = (b1 > bright) ? b1 : bright;
	contra = (c2 < yl) ? c2 : yl;
	contra = (c1 > contra) ? c1 : contra;
    }

    /* Initialize contrast and brightness */
    if (ch=='I' ) {
	bright = 0.5;
	contra = 1.0;
	xs = 0.5;
	yl = 1.0;
    }

    /* Change palett */
    switch (ch) {
    case '-':
	sign = -sign;
	break;
    case '1':
	pal = 1;
	break;
    case '2':
	pal = 2;
	break;
    case '3':
	pal = 3;
	break;
    case '4':
	pal = 4;
	break;
    case '5':
	pal = 5;
	break;
    case 'P':
	pal = 1 + pal % 5;
	break;
    }
    palett(pal, sign*contra, bright);
}

/**********************/
/* print help message */
/**********************/

void explain(void)
{

    printf(" Cursor usages:                                  \n");
    printf("                                                 \n");
    printf("    Key A: adjusts contrast and brightness       \n");
    printf("    Key I: initialize contrast and brightness    \n");
    printf("    Key X: exit                                  \n");
    printf("                                                 \n");
    printf("    Key -: reverse contrast                      \n");
    printf("    Key 1: gray scale                            \n");
    printf("    Key 2: rainbow                               \n");
    printf("    Key 3: rainbow inverse                       \n");
    printf("    Key 4: heat                                  \n");
    printf("    Key 5: weird IRAF                            \n");
    printf("    Key 6: AIPS                                  \n");
    printf("                                                 \n");
    printf("    Key E: expand box                            \n");
    printf("    Key Z: zoom box                              \n");
    printf("                                                 \n");
    printf("    Key S: draw face-on shot                     \n");
    printf("    Key P: draw PV-diagram from cursor point     \n");
    printf("    Key V: draw Velocity field from cursor point \n\n");

}

/*****************/
/* get body data */
/*****************/

void getdata(char *fname)
{

    int i;
    FILE *fp;

    /*                 */
    /* Open input file */
    /*                 */

    if ((fp = fopen(fname,"r")) == NULL) {
	fclose(fp);
	printf("error: file \"%s\" cannot be opened/n", fname);
	exit(1);
    }

    /*                     */
    /* Read number of body */
    /*                     */

    fread(&nbody, sizeof(int), 1, fp);

    /*                                    */
    /* Allocate memory to store body data */
    /*                                    */

    m = malloc(nbody * sizeof(real));
    x = malloc(nbody * sizeof(real));
    y = malloc(nbody * sizeof(real));
    vx = malloc(nbody * sizeof(real));
    vy = malloc(nbody * sizeof(real));
    xt = malloc(nbody * sizeof(real));
    yt = malloc(nbody * sizeof(real));
    zt = malloc(nbody * sizeof(real));
    vxt = malloc(nbody * sizeof(real));
    vyt = malloc(nbody * sizeof(real));
    vzt = malloc(nbody * sizeof(real));

    /*                */
    /* Read body data */
    /*                */

    fread(&ndim, sizeof(int), 1, fp);
    fread(&mscale, sizeof(real), 1, fp);
    fread(&rscale, sizeof(real), 1, fp);
    fread(&tscale, sizeof(real), 1, fp);

    fread(&tnow, sizeof(real), 1, fp);
    tnow *= tscale;

    for (i=0; i<nbody; i++) {
	fread(&m[i], sizeof(real), 1, fp);
	m[i] *= mscale;
    }

    for (i=0; i<nbody; i++) {
	fread(&x[i], sizeof(real), 1, fp);
	fread(&y[i], sizeof(real), 1, fp);
	x[i] *= rscale;
	y[i] *= rscale;
    }
    
    for (i=0; i<nbody; i++) {
	fread(&vx[i], sizeof(real), 1, fp);
	fread(&vy[i], sizeof(real), 1, fp);
	vx[i] *= rscale / tscale * 9.77792e8;
	vy[i] *= rscale / tscale * 9.77792e8;
    }

    /*                  */
    /* Close input file */
    /*                  */

    fclose(fp);
}
