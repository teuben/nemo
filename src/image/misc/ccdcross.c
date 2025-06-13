/* 
 * CCDSTACK: cross-correlate images, first one in the reference image
 *
 *   11-apr-2022:    derived from ccdstack
 *
 * @todo   gaussian fit to peak in corr image?
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <extstring.h>
#include <ctype.h>

string defv[] = {
  "in=???\n       Input image files, first image sets the WCS",
  "out=.\n        Output cross correlated image",
  "center=\n      X-Y Reference center (0-based pixels)",
  "box=\n         Half size of correlation box",
  "n=\n          Half size of box inside correlation box to find center",
  "clip=\n        Only use values above this clip level from input",
  "clop=\n        Clip from correllation image to finder center",
  "bad=0\n        bad value to ignore",
  "VERSION=0.6\n  6-jun-2024 PJT",
  NULL,
};

string usage = "cross correlate images to find the offset between images";


#ifndef HUGE
# define HUGE 1.0e35
#endif

#define MAXIMAGE 100

imageptr iptr[MAXIMAGE];	/* pointers to (input) images */
imageptr optr = NULL;           /* could do array */

real     badval;
int      nimage;                /* actual number of input images */
int      center[2];              /* center of reference point image */
int      box;
bool     Qclip, Qclop;
real     clip[MAXIMAGE];       // @todo
int      nclip;
real     clop;


local void do_cross(int l0, int l, int n);

extern int minmax_image (imageptr iptr);

void nemo_main()
{
    string *fnames;
    stream  instr;                      /* input files */
    stream  outstr;                     /* output file */
    int     nx, ny, nc;
    int     ll, n;
    
    badval = getrparam("bad");
    Qclip = hasvalue("clip");
    if (Qclip) clip[0] = getrparam("clip");
    Qclop = hasvalue("clop");
    if (Qclop) clop = getrparam("clop");
 
    fnames = burststring(getparam("in"), ", ");  /* input file names */
    nimage = xstrlen(fnames, sizeof(string)) - 1;
    if (nimage > MAXIMAGE) error("Too many images %d > %d", nimage, MAXIMAGE);
    dprintf(0,"Using %d images\n",nimage);

    instr   = stropen(fnames[0],"r");    /* open file */
    iptr[0] = NULL;        /* make sure to init it right */
    read_image (instr, &iptr[0]);
    dprintf(0,"Image %d: pixel size %g %g   minmax %g %g [%s]\n",
	    0, Dx(iptr[0]),Dy(iptr[0]),MapMin(iptr[0]),MapMax(iptr[0]),fnames[0]);
    strclose(instr);

    if (hasvalue("center")) {
      //  @todo  :   allow "max" and "ref" ???
      nc = nemoinpi(getparam("center"),center,2);
      if (nc != 2) error("center= needs 2 values");
    } else {
      center[0] = (Nx(iptr[0])-0)/2;
      center[1] = (Ny(iptr[0])-0)/2;
      dprintf(0,"center=%d,%d\n",center[0],center[1]);
    }

    if (hasvalue("box"))
      box = getiparam("box");
    else {
      box = (Nx(iptr[0])-1)/2;
      if (Ny(iptr[0]) < Nx(iptr[0]))
	  box = (Ny(iptr[0])-1)/2;
      // @todo this can fail if off-center is used
    }
    if (hasvalue("n"))
      n = getiparam("n");
    else
      n = box/2;

    nx = ny = 2*box + 1;
    dprintf(0,"Full box size %d (square), n=%d\n",nx,n);
    
    optr = NULL;
    create_cube(&optr, nx, ny, 1);
    CRPIX1(optr) = (nx-1)/2;  // 0 based!
    CRPIX2(optr) = (ny-1)/2;
    //CDELT1(optr) = CDELT1(iptr[0]);
    //CDELT2(optr) = CDELT2(iptr[0]);
    outstr = stropen(getparam("out"),"w");  /* open output file first ... */

    for (ll=1; ll<nimage; ll++) {
        instr   = stropen(fnames[ll],"r");    /* open file */
        iptr[ll] = NULL;        /* make sure to init it right */
        read_image (instr, &iptr[ll]);
	dprintf(0,"Image %d: pixel size %g %g   minmax %g %g [%s]\n",
		ll, Dx(iptr[ll]),Dy(iptr[ll]),MapMin(iptr[ll]),MapMax(iptr[ll]),fnames[ll]);
        strclose(instr);        /* close input file */
	do_cross(0,ll,n);
    }

    write_image(outstr,optr); 
    strclose(outstr);

}


/* 
 *  combine input maps into an output map  --
 *
 */
local void do_cross(int l0, int l, int n)
{
    real   sum, imval;
    int    i, j, ix, iy, nx, ny;
    int    ix1,iy1,ix2,iy2;
    int    cx,cy;
    int    badvalues;
    
    badvalues = 0;		/* count number of bad operations */

    nx = Nx(iptr[l0]);
    ny = Ny(iptr[l0]);

    dprintf(0,"%d %d  %d  %d %d\n", nx,ny,box,l0,l);
    cx = center[0];
    cy = center[1];

    // make the cross corr image
    for (j=-box; j<=box; j++) {
      for (i=-box; i<=box; i++) {
	sum = 0.0;
	for (iy=-box; iy<=box; iy++) {
	  iy1 = cy + iy;
	  iy2 = iy1 + j;
	  if (iy1 < 0 || iy2 < 0 || iy1 >= ny || iy2 >= ny) continue;	  
	  for (ix=-box; ix<=box; ix++) {
	    ix1 = cx + ix;;
	    ix2 = ix1 + i;
	    if (ix1 < 0 || ix2 < 0 || ix1 >= nx || ix2 >= nx) continue;
	    // @todo   handle bad values
	    if (Qclip) {
	      if (CubeValue(iptr[l0],ix1,iy1,0) < clip[0]) continue;
	      if (CubeValue(iptr[l], ix2,iy2,0) < clip[0]) continue;
	    }
	    sum += CubeValue(iptr[l0],ix1,iy1,0) * CubeValue(iptr[l],ix2,iy2,0);
	  }
	}
	CubeValue(optr,i+box,j+box,0) = sum;
      }
    }

    // find where the peak is
    int ix0=0, iy0=0;
    minmax_image(optr);
    dprintf(0,"New min and max in correlation image are: %f %f\n",MapMin(optr) ,MapMax(optr) );
    for (iy=0; iy<Ny(optr); iy++)
      for (ix=0; ix<Nx(optr); ix++)
	if (CubeValue(optr,ix,iy,0) == MapMax(optr)) {
	  ix0 = ix;
	  iy0 = iy;
	  dprintf(0,"Max cross %g @ %d %d\n",MapMax(optr),ix,iy);
	}

    // find the weighted mean center
    real sum0=0, sumx=0, sumy=0, sumxx=0, sumxy=0, sumyy=0;
    real smin=0, smax=0;
    int dx, dy;
    sum = 0;
    for (dy=-n; dy<=n; dy++) {        // look around (ix0,iy0) by just a few pixels 
      iy = iy0+dy;  
      for (dx=-n; dx<=n; dx++) {      // to find the 
	ix = ix0+dx;
	imval = CubeValue(optr,ix,iy,0);
	if (Qclop && imval < clop) continue;
	dprintf(1,"clop %g %g\n",sum0,imval);
	if (sum0 == 0)
	  smin = smax = imval;
	smin = MIN(smin, imval);
	smax = MAX(smax, imval);
	sum0  = sum0  + 1;
	sum   = sum   + imval;
	sumx  = sumx  + imval * dx;
	sumy  = sumy  + imval * dy;
	sumxx = sumxx + imval * dx*dx;
	sumxy = sumxy + imval * dx*dy;
	sumyy = sumyy + imval * dy*dy;
      }
    }
    sumx /= sum;   sumy /= sum;
    sumxx = sqrt(sumxx/sum - sumx*sumx);
    sumyy = sqrt(sumyy/sum - sumy*sumy);
    dprintf(0,"MinMax in %g pts: %g %g, use clop= to clip the corr out image\n", sum0, smin, smax);
    dprintf(0,"X,Y mean: %g %g\n", sumx,  sumy);
    dprintf(0,"X,Y sig:  %g %g\n", sumxx, sumyy);    
    real xcen = ix0-box+sumx;
    real ycen = iy0-box+sumy;
    dprintf(0,"Center at: %g %g from center\n",xcen,ycen);
    // @todo (ix0-xcen-box, iy0-ycen-box) is a value expected to be around (0,0) for no offsets
    // 
    printf("%g %g  %g %g  %d %d  %g\n",xcen,ycen,sumx,sumy, ix0,iy0, MapMax(optr));	   
    
    if (badvalues)
    	warning("There were %d bad operations in dofie",badvalues);
}



