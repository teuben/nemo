/****************************************************************************/
/* MOVIE.C: routines to make gif images                                     */
/* Copyright (c) 2000 by Jin Koda, Tokyo, JAPAN.                            */
/****************************************************************************/

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "/usr/local/pgplot/cpgplot.h"

void getdata(FILE *);
void transform(int, real *, real *, real, real);
void calcmesh(int, real *, real *, real *, real *, real, real);
void plot(real *, int , int , int, real, real);
void palett(int, real, real);
void setvp(void);

#define NMESH1 128
#define NMESH2 NMESH1

int nbody, ndim;
real tnow, *m, *x, *y, *vx, *vy;
real mscale, rscale, tscale;
int flag = 0;

int main(int argc, string argv[])
{

    int istep;
    real omgb, theta, pscale;
    real meshd[NMESH2][NMESH1];
    char fname[32], fnum[32], s[10];
    FILE *fp1;

    /*                   */
    /* Set pattarn speed */
    /*                   */ 
    if (argc < 2)
        printf("No output file name\n"); 

    if (argc < 3) {
	omgb = 0.0;
	printf("No pattarn speed is set.\n");
    } else {
	sscanf(argv[2], "%f", &omgb);
	printf("Pattarn Speed  = %.2f [km/s/kpc]\n", omgb);
    }
    omgb = 1.02271202e-9 * omgb;

    /*                 */
    /* Loop over steps */
    /*                 */
    nbody = -1;
    for (istep=0; istep<1000; istep++) {

	/*                     */
	/* Set input file name */
	/*                     */

        printf("istep = %.4d: ", istep);
	sprintf(fnum, "%.3d", istep);
 	strcpy(fname, argv[1]);
        strcat(fname, "_");
        strcat(fname, fnum);
        printf("input file = %s: ", fname);

	/*                         */
	/* Open file and read data */
	/*                         */
	if ((fp1 = fopen(fname,"r")) == NULL) {
	    printf("\n\n    !!! Finished !!!\n\n");
	    exit(1);
	}
	getdata(fp1);
	printf("tnow = %4.1f Myr\n", tnow/1.0e6);
	fclose(fp1);

	/*                                     */
	/* Open pgdevice; 0: /XWINDOW, 1: /GIF */
	/*                                     */
	if (1) {
	    strcat(fname, ".gif/GIF");
	    cpgbeg(0, fname, 1, 1);
	} else {
	    cpgbeg(0, "/xwindow", 1, 1);
	}

	/*                      */
	/* Set bar horizontally */
	/*                      */
	theta = omgb * tnow;
	transform(nbody, x, y, -theta,    1.0);
	transform(nbody, vx, vy, -theta,  1.0);

	/*                                   */
	/* Plot; 0: Snapshot, 1: P-V diagram */
	/*                                   */
	pscale = 2.0 * rscale;
	if (0) {
	    calcmesh(nbody, x, vy, m, (real *) meshd, pscale, 600.0);
	    plot((real *) meshd, NMESH1, NMESH2, 1, pscale, 600.0);
	} else {
	    calcmesh(nbody, x, y, m, (real *) meshd, pscale, pscale);
	    plot((real *) meshd, NMESH1, NMESH2, 1, pscale, pscale);
	}
    }

}

void transform(int nn, real *xx, real *yy, real tt, real sign)
{

    int i;
    real xx0, yy0;

    for(i=0; i<nn; i++) {
	xx0 = xx[i];
	yy0 = yy[i];
	xx[i] = sign * (cos(tt) * xx0 - sin(tt) * yy0);
	yy[i] = sign * (sin(tt) * xx0 + cos(tt) * yy0);
    }
}

void calcmesh(int nbody, real *xx, real *yy, real *weight,
	      real *mm, real xsize, real ysize)
{

#define F(j,i) mm[(j) * NMESH1 + (i)]

    int i, j, p;
    real xx0, yy0, ddx, ddy, ttx, tty, dx, dy, rsq;

    dx = xsize / ((real) NMESH1);
    dy = ysize / ((real) NMESH2);

    for (i=0; i<NMESH1; i++)
	for (j=0; j<NMESH2; j++)
	    F(i,j) = 0.0;	

    for (p=0; p<nbody; p++) {
	xx0 = (xx[p] + 0.5 * xsize - 0.5 * dx)/dx;
	yy0 = (yy[p] + 0.5 * ysize - 0.5 * dy)/dy;
	i = (int) xx0;
	j = (int) yy0;
	if (i >= 0 && i < NMESH1-1 && j >= 0 && j < NMESH2-1) {
	    ddx = xx0 - (real) i;
	    ddy = yy0 - (real) j;
	    ttx = 1.0 - ddx;
	    tty = 1.0 - ddy;
	    F(j  ,i  ) = F(j  ,i  ) + ttx*tty*weight[p];
	    F(j  ,i+1) = F(j  ,i+1) + ddx*tty*weight[p];
	    F(j+1,i  ) = F(j+1,i  ) + ttx*ddy*weight[p];
	    F(j+1,i+1) = F(j+1,i+1) + ddy*ddy*weight[p];
	}
    }
}

void plot(real *data, int ni, int nj, int nk, real bsize1, real bsize2)
{
 
    int   ii, jj, kk, ll, c1, c2, nc, nij;
    real *f;
    real fmin, fmax, tr[6], contra, bright;
    char  val[32], s[32];

    /*                            */
    /*  Open device for graphics. */
    /*                            */
    cpgqcir(&c1, &c2);                /* inquire color index range  */
    nc = c2-c1+1;                     /* define the number of color */

    /*                           */
    /*  Summing up along z-axis. */
    /*                           */
    nij  = ni  * nj;
    f    = malloc( sizeof(real) * nij );     /* allocate memory    */
    for (ii=0; ii<nij; ii++) *(f+ii) = 0.0;   /* initialize array f */
    for (kk=0; kk<nk ; kk++)
	for (ii=0; ii<nij; ii++)
	    *(f+ii) += *(data+nij*kk+ii);         /* summation along z  */

    /*                                  */
    /*  Maximum and minimum in array F. */ 
    /*                                  */
    fmin =  1.0e9;
    fmax = -1.0e9;
    for (ii=0; ii<nij; ii++) {
	fmin = (*(f+ii) < fmin) ? *(f+ii) : fmin;
	fmax = (*(f+ii) > fmax) ? *(f+ii) : fmax;
    }
    fmin = 0.0;

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
    cpgpap(6.0, 1.2);
    cpgpage();
    setvp();
    cpgwnad(0.0, 1.0+ni, 0.0, 1.0+nj);
    

    /*                       */
    /* Set up the color map. */
    /*                       */
    bright = 0.5;
    contra = 1.0;
    palett(3, contra, bright);
    
    /*                           */
    /* Draw the map with PGIMAG. */
    /*                           */
    
    cpgimag(f,ni,nj,1,ni,1,nj,fmin,fmax,tr);
    /*cpgline(NMESH1, pospv, velpv);*/

    /*cpgsci(4);
    cpgline(NMESH1/2, &pospv[NMESH1/4], &velenv[NMESH1/4]);
    cpgsci(1);*/
    
    /*                    */
    /* Annotate the plot. */
    /*                    */

    cpgsch(0.6);
    cpgswin(-bsize1/2.0, bsize1/2.0, -bsize2/2.0, bsize2/2.0);
    cpgsch(1.0);
    cpgbox("bcntsi",0.0,0,"bcntsi",0.0,0);
    cpgsch(1.5);
    sprintf(s, "%4.1f Myr\0", tnow/1.0e6);
    cpgmtxt("t",1.0,0.5,0.5,s);
    cpgsch(1.0);
    /*cpgmtxt("b",3.0,0.5,0.,"kpc");*/

    /*               */
    /* Draw a wedge. */
    /*               */

    /* cpgwedg("bi",4.0,5.0,fmin,fmax,"pixel value"); */
    cpgsch(1.0); 

    /*                 */
    /* Release memory. */
    /*                 */
    free(f);
    cpgend();
}


/*************************/
/* Set up a color table. */
/*************************/

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
    case 3: /* heat */
	cpgctab(hl, hr, hg, hb,  5, contra, bright);
	break;
    case 4: /* weird IRAF */
	cpgctab(wl, wr, wg, wb, 10, contra, bright);
	break;
    case 5: /* AIPS */
	cpgctab(al, ar, ag, ab, 10, contra, bright);
	break;
    }
}

/********************/
/* Set up viewport. */
/********************/

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

void getdata(FILE *fpin)
{
    int i, k, nb=0;

    fread(&nb, sizeof(int), 1, fpin);

    if (flag == 0) {
        nbody = nb;
        m = malloc(nbody * sizeof(real));
        x = malloc(nbody * sizeof(real));
        y = malloc(nbody * sizeof(real));
        vx = malloc(nbody * sizeof(real));
        vy = malloc(nbody * sizeof(real));
        flag++;
    } else if (nb != nbody) {
        free(m);
        free(x);
        free(y);
        free(vx);
        free(vy);
        exit(1);
    }

    fread(&ndim, sizeof(int), 1, fpin);
    fread(&mscale, sizeof(real), 1, fpin);
    fread(&rscale, sizeof(real), 1, fpin);
    fread(&tscale, sizeof(real), 1, fpin);

    fread(&tnow, sizeof(real), 1, fpin);
    tnow *= tscale;

    for (i=0; i<nbody; i++) {
	fread(&m[i], sizeof(real), 1, fpin);
	m[i] *= mscale;
    }

    for (i=0; i<nbody; i++) {
	fread(&x[i], sizeof(real), 1, fpin);
	fread(&y[i], sizeof(real), 1, fpin);
	x[i] *= rscale;
	y[i] *= rscale;
    }
    
    for (i=0; i<nbody; i++) {
	fread(&vx[i], sizeof(real), 1, fpin);
	fread(&vy[i], sizeof(real), 1, fpin);
	vx[i] *= rscale / tscale * 9.77792e8;
	vy[i] *= rscale / tscale * 9.77792e8;
    }
}
