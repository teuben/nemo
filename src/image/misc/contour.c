/*
 * CONTOUR: contouring drawing or extraction routine 
 *	
 *	    based on a PDP-11 version by M.D. Jones
 *              Ref. [conrec.bas - Byte  Jun 87, pp. 143]
 *          documented by PJT
 *
 *        4)                      3)
 *          x--------------------x      Points 1..4 form the 'cell', and 0
 *          |  .             .   |      the center of a 'cell', all seen from
 *          |     .       .      |      pixel 1. Hence looping over the image
 *          |        . . 0)      |      must be from pixel ix=0..nx-1 and 
 *          |       .   .        |      iy=0..ny-1
 *          |    .         .     |
 *          | .               .  |
 *          x--------------------x
 *        1)                      2)
 *
 *
 *      17-Jun-87   P.J.Teuben - created with mongo-like calling convention
 *                  Array 'a' needs to be 1-D. a[][] doesn't work
 *      29-jun-87   removed ydim as second parameter  ...
 *       8-jul-87   bug: one column (x-pos) too many
 *       1-feb-89   FORDEF/CDEF definition of matrices implemented
 *      16-mar-90   made GCC happier    PJT
 *      10-sep-90   added the undef feature PJT
 *       7-mar-92   gcc2.0 happy        PJT
 *      22-aug-92   using idx, m1, r1,r2 for faster computation    PJT
 *                  also corrected bug: ish[0] had never been set
 *                  documented the code more, plus a diagram.
 *		    Can't figure the speed decrease out yet, program
 *                  is 10-20% slower now...
 *	30-mar-97   finally fixed the -DSINGLEPREC bug with lineto_proc
 *      10-jan-03   SGN -> SIGN
 *	29-jul-03   SIGN -> SGN 
 */

/*CREF
 *
 * Bourke, P.D., 1987, A contouring subroutine: Byte, 12, p. 143-150.
 *
 * Bregoli, L.J., 1982, A BASIC plotting subroutine sophisticated 
 * plotting with your MX-80: Byte, 7(3), p. 142-156.
 *
 * Daulton, T., 1987, Three-dimensional perspective plotting: 
 * Byte, 12(14), p. 307-314.
 * 
 * Enns, S., 1986, Free-form curves on your micro: 
 * Byte, 11(13), p. 225-230.
 * 
 * Simons, S.L.Jr., 1983, Make fast and simple plots on a 
 * microcomputer: Byte, 8, p. 487-492.
 * 
 * Watson, D. 1992,
 * CONTOURING: A guide to the analysis and display of spatial data:
 * Pergamon Press, [ISBN 0 08 040286 0]
 */


#include <stdinc.h>
#include <matdef.h>   /* this one is purely to get the FORDEF/CDEF define */

local   int  im[4]= {0, 1, 1, 0};       /* vertex -> x-offset lookup table  */
local   int  jm[4]= {0, 0, 1, 1};       /* vertex -> y-offset lookup table  */
  /* the following lookup table handles all possible levels in a triangle   */
  /* where the three values at the vertices could be less then, equal or    */
  /* more than the requested contour level.                                 */
local   int  castab[3][3][3] = {0,0,8,0,2,5,7,6,9,
                                0,3,4,1,3,1,4,3,0,
                                9,6,7,5,2,0,8,0,0};

local real   undef_value=0.0;
local int    undef_mode=0;     /* default: there are not undef's */

/* 
 * A helper proc: otherwise (if real==float) these are passed as 
 * double's by e.g. gcc even if cast'd ! 
 */
typedef void (*lineto_proc)(real,real,real,real);

int contour (real *a, int nx, int ny, real *z, int nc, 
             real xmin, real ymin, real xmax, real ymax, lineto_proc lineto)
{
    real h[5];                  /* height of the cell above contour */
    int  ish[5];                /* sign of h[] */
    real xh[5], yh[5];          /* cell coordinates */
    real deltax, deltay;        /* size of plot area, assuming equidistant */
    real x1,x2,y1,y2;
    real v1,v2,v3,v4;
    int  i,j,k,m,m1,m2,m3, icase, idx, ncnt=0;
    real dmin, dmax, d1, d2;    /* keeps datamin/datamax, scaling factors */

    if (nc<=0 || nx <=0 || ny <= 0)     /* if no contours, or map */
        return 0;                       /* return immediately */
    if (a==NULL) return 0;              /* also return for other */
    if (z==NULL) return 0;              /* obvious errors */

    deltax = (xmax-xmin)/(nx-1);        /* phyical pixel size (assumed */
    deltay = (ymax-ymin)/(ny-1);        /* linear pixels               */
    if (deltax==0 || deltay==0) {
        warning("contour: illegal range(s) xmin=xmax or ymin=ymax");
        return 0;
    }
    for (k=1; k<nc; k++)
        if (z[k-1] > z[k]) {
            warning("contour: (#%d) contours need to be sorted",k);
            return 0;
        }

    for (j=0; j<ny-1; j++) {            /* loop over all cells */
        for (i=0; i<nx-1; i++) {        /* in the image */

                        /* get value of 'a' at the 4 corners corners */
#if defined(CDEF)
            idx = i*ny+j;               /* offset for lower left cell */
            v1 = *(a + idx);            /* CDEF of matrix */
            v2 = *(a + idx + ny);    
            v3 = *(a + idx + 1);       
            v4 = *(a + idx + ny + 1);
#else
            idx = j*nx+i;               /* offset for lower left cell */
            v1 = *(a + idx);            /* FORDEF of matrix */
            v2 = *(a + idx + nx);    
            v3 = *(a + idx + 1);       
            v4 = *(a + idx + nx + 1);
#endif
            if (undef_mode &&           /* if undef to be used */
                    (v1==undef_value || v2==undef_value || 
                     v3==undef_value || v4==undef_value))
                continue;               /* don't contour here */

            if (v1<v2)   dmin=v1;          /* find lowest vertex */
            else         dmin=v2;
            if (v3<dmin) dmin=v3;
            if (v4<dmin) dmin=v4;

            if (v1>v2)   dmax=v1;          /* find highest vertex */
            else         dmax=v2;
            if (v3>dmax) dmax=v3;
            if (v4>dmax) dmax=v4;

            if ((dmax<z[0]) || (dmin>z[nc-1]))  /* using sorted contours */
                continue;                     /* none in cell -> next cell */

            for (k=0; k<nc; k++) {      /* loop for all contour values */
                if ((z[k]<dmin) || (z[k]>dmax)) /* check if this contour */
                    continue;                   /* no: next one */
                for (m=4; m>=0; m--) {
                    if (m==0) {          /* m=0: special case, handled last */
                        h[0]=0.25*(h[1]+h[2]+h[3]+h[4]);    /* center point */
                        xh[0]=xmin + (i+0.5)*deltax;
                        yh[0]=ymin + (j+0.5)*deltay;
                        ish[0] = SGN(h[0]) + 1;
                    } else {
                        m1 = m-1;
#if defined(CDEF)
                        h[m] = *(a+(i+im[m1])*ny+j+jm[m1]) - z[k];
#endif
#if defined(FORDEF)
                        h[m] = *(a+(i+im[m1])+(j+jm[m1])*nx) - z[k];
#endif
                        xh[m] = xmin + (i+im[m1])*deltax;
                        yh[m] = ymin + (j+jm[m1])*deltay;
                    }
                    ish[m] = SGN(h[m]) + 1;   /* make sure ish = 0,1,2 */
                }  /* end m-loop */

                for (m=1; m<=4; m++) {  /* handle each of 4 triangles in cell */

                    /* m2 is always the cell centerpoint */
                    /* m1 loops around, and m3 is the next one around cell */

                    m1 = m; m2 = 0; m3 = (m==4) ? 1 : m+1;

                    icase = castab[ish[m1]][ish[m2]][ish[m3]];
                    switch (icase) {
                      case 0:   /* nothing to plot */
                        break;
                      case 1:   /* line between vertices m1 and m2 */ 
                        x1=xh[m1]; y1=yh[m1]; x2=xh[m2]; y2=yh[m2];
                        break;
                      case 2:  /* line between vertices m2 and m3 */
                        x1=xh[m2]; y1=yh[m2]; x2=xh[m3]; y2=yh[m3];
                        break;
                      case 3:  /* line between vertices m3 and m1 */
                        x1=xh[m3]; y1=yh[m3]; x2=xh[m1]; y2=yh[m1];
                        break;
                      case 4:   /* from m1 to opposing vertex */
                        d2=1/(h[m3]-h[m2]);
                        x1=xh[m1]; y1=yh[m1];
                        x2=(h[m3]*xh[m2]-h[m2]*xh[m3])*d2;
                        y2=(h[m3]*yh[m2]-h[m2]*yh[m3])*d2;
                        break;
                      case 5:   /* from m2 to opposing vertex */
                        d2=1/(h[m1]-h[m3]);
                        x1=xh[m2]; y1=yh[m2];
                        x2=(h[m1]*xh[m3]-h[m3]*xh[m1])*d2;
                        y2=(h[m1]*yh[m3]-h[m3]*yh[m1])*d2;
                        break;
                      case 6:   /* from m3 to opposing vertex */
                        d2=1/(h[m2]-h[m1]);
                        x1=xh[m3]; y1=yh[m3];
                        x2=(h[m2]*xh[m1]-h[m1]*xh[m2])*d2;
                        y2=(h[m2]*yh[m1]-h[m1]*yh[m2])*d2;
                        break;
                      case 7:   /* cross-over 'parallel' to m1-m3 */
                        d1=1/(h[m2]-h[m1]);
                        d2=1/(h[m3]-h[m2]);
                        x1=(h[m2]*xh[m1]-h[m1]*xh[m2])*d1;
                        y1=(h[m2]*yh[m1]-h[m1]*yh[m2])*d1;
                        x2=(h[m3]*xh[m2]-h[m2]*xh[m3])*d2;
                        y2=(h[m3]*yh[m2]-h[m2]*yh[m3])*d2;
                        break;
                      case 8:   /* cross-over 'parallel' to m1-m2 */
                        d1=1/(h[m3]-h[m2]);
                        d2=1/(h[m1]-h[m3]);
                        x1=(h[m3]*xh[m2]-h[m2]*xh[m3])*d1;
                        y1=(h[m3]*yh[m2]-h[m2]*yh[m3])*d1;
                        x2=(h[m1]*xh[m3]-h[m3]*xh[m1])*d2;
                        y2=(h[m1]*yh[m3]-h[m3]*yh[m1])*d2;
                        break;
                      case 9:   /* cross-over 'parallel' to m2-m3 */
                        d1=1/(h[m1]-h[m3]);
                        d2=1/(h[m2]-h[m1]);
                        x1=(h[m1]*xh[m3]-h[m3]*xh[m1])*d1;
                        y1=(h[m1]*yh[m3]-h[m3]*yh[m1])*d1;
                        x2=(h[m2]*xh[m1]-h[m1]*xh[m2])*d2;
                        y2=(h[m2]*yh[m1]-h[m1]*yh[m2])*d2;
                        break;
                      default:
                        error("contour: impossible triangle case %d",icase);
                        break;
                    }
                    if (icase>0) {
			lineto(x1,y1,x2,y2); /* return/draw segment */
			ncnt++;
		    }
                }  /* end m-loop over triangles */
            } /* end k-loop of contours */
        } /* end x-loop over all cells in a row in the image */
    } /* and y-loop over all rows in the image */
    return 1;
}

void contour_setdef(int mode,real value)
{
    if (mode) {
        undef_mode = 1;
        undef_value = value;
    } else
        undef_mode = 0;
}



/* believe it or not: there's no TESTBED in contour ! */
