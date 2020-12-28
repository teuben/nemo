/* 
 * CCDSUB: get a subimage from an image
 *
 *	quick and dirty: 10-aug-93
 *	something...	  6-sep-93
 *	n[xy]aver=	  5-nov-93 	puzzling.....
 *		something is wrong here                     
 * 1.3  added some z stuff for Martin Bureau       1-may-2002      PJT
 * 1.4  added reorder=  (args, largely unfinished)
 * 1.5  added a more proper WCS when simple subsetting is done     PJT
 * 2.0  rearranged a lot, removed useless options, enabled others  PJT
 * 2.1  added moving=t averaging for nxaver only (for now)         PJT
 * 2.2  fixed WCS on output 

    TODO:  wcs is wrong on output
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
  "in=???\n       Input image file",
  "out=???\n      Output image file",
  "x=\n		  X coordinates to select (1..nx)",
  "y=\n		  Y coordinates to select (1..ny)",
  "z=\n		  Z coordinates to select (1..nz)",
  "nxaver=1\n	  Number X to aver (size remains same)",
  "nyaver=1\n	  Number Y to aver (size remains same)",
  "nzaver=1\n	  Number Z to aver (size remains same)",
  "centerbox=\n   xf,yf,zf - if given, use this centered fraction of this axis",
  "dummy=t\n      Retain dummy axis?",
  "reorder=\n     New coordinate ordering",
  "moving=f\n     Moving average in n{x,y,z}aver= ?",
  "VERSION=2.4\n  24-dec-2020 PJT",
  NULL,
};

string usage = "sub/average of an image, with reorder option";

string cvsid="$Id$";



#define MAXDIM 16384

int  ix[MAXDIM], iy[MAXDIM], iz[MAXDIM];

local int  ax_index(string , int , int , int *);
local void ax_shift(imageptr iptr);
local void ax_copy(imageptr i0, imageptr i1);
local void ax_swap_xy(imageptr iptr);
local void ax_swap_yz(imageptr iptr);
local void ax_swap_xz(imageptr iptr);

#define SWAP(a,b,tmp) {tmp=a; a=b; b=tmp;}
#define LOOP(i,n)     for(i=0;i<n;i++)
#define CV(p,i,j,k)   CubeValue(p,i,j,k)

void nemo_main(void)
{
    stream  instr, outstr;
    int     nx, ny, nz;        /* size of scratch map */
    int     nx1,ny1,nz1;
    int     nxaver, nyaver,nzaver;
    int     i,j,k, i0,j0,k0, i1,j1,k1, l;
    int     ncb, n1,n2;
    real    centerbox[3];
    imageptr iptr=NULL, iptr1=NULL;      /* pointer to images */
    real    sum, tmp, zzz;
    real    *row;
    bool    Qreorder = FALSE;
    bool    Qdummy, Qsample, Qmoving;
    string  reorder;

    instr = stropen(getparam("in"), "r");
    nxaver=getiparam("nxaver");
    nyaver=getiparam("nyaver");
    nzaver=getiparam("nzaver");
    Qdummy = getbparam("dummy");
    Qmoving = getbparam("moving");

    nx1 = nemoinpi(getparam("x"),ix,MAXDIM);
    ny1 = nemoinpi(getparam("y"),iy,MAXDIM);
    nz1 = nemoinpi(getparam("z"),iz,MAXDIM);
    if (nx1<0 || ny1<0 || nz1<0) error("Error parsing x,y,z=");
    Qsample = nx1>0 || ny1>0 || nz1>0;

    ncb = nemoinpr(getparam("centerbox"),centerbox,3);

    read_image( instr, &iptr);

    nx = Nx(iptr);	                   /* old cube size */
    ny = Ny(iptr);      
    nz = Nz(iptr);
    
    if (ncb > 0) {               // See if centerbox was used
      Qsample = TRUE;
      //warning("Using centerbox %d",ncb);
      if (ncb>0) {
	n1 = (int)rint(0.5*nx*(1 - centerbox[0]));
	n2 = (int)rint(0.5*nx*(1 + centerbox[0]));
	if (n1<1) n1=1;
	dprintf(1,"x: %d:%d\n",n1,n2);
	nx1 = n2-n1+1;
	for (i=0; i<nx1; i++)
	  ix[i] = i + n1;
      }
      if (ncb>1) {
	n1 = (int)rint(0.5*ny*(1 - centerbox[1]));
	n2 = (int)rint(0.5*ny*(1 + centerbox[1]));
	if (n1<1) n1=1;	
	dprintf(1,"y: %d:%d\n",n1,n2);	
	ny1 = n2-n1+1;
	for (i=0; i<ny1; i++)
	  iy[i] = i + n1;
      }
      if (ncb>2) {
	n1 = (int)rint(0.5*nz*(1 - centerbox[2]));
	n2 = (int)rint(0.5*nz*(1 + centerbox[2]));
	if (n1<1) n1=1;	
	dprintf(1,"z: %d:%d\n",n1,n2);	
	nz1 = n2-n1+1;
	for (i=0; i<nz1; i++)
	  iz[i] = i + n1;
      }
    }
    if (Qsample) warning("Sampling will be done");
    
    nx1 = ax_index("x",nx,nx1,ix);         /* initialize new cube axes */
    ny1 = ax_index("y",ny,ny1,iy);         /* before any reordering */
    nz1 = ax_index("z",nz,nz1,iz);
    dprintf(1,"old nx,y,z: %d %d %d\n",nx ,ny ,nz );    
    dprintf(1,"new nx,y,z: %d %d %d\n",nx1,ny1,nz1);

    Qreorder = hasvalue("reorder");
    if (Qreorder) {
      reorder = getparam("reorder");
      if (strlen(reorder) != 3) error("Reorder must have 3 letters (e.g. xzy)");
    } 

    outstr = stropen(getparam("out"), "w");

    if (nxaver>1 && Qmoving) {
      warning("work in progress; nxaver=%d moving=t", nxaver);   /* this can be done a lot more efficient */
      row = (real *) allocate(nx*sizeof(real));
      LOOP(k,nz) {
	LOOP(j,ny) {
	  LOOP(i,nx) { 
	    row[i] = CV(iptr, i, j, k);
	    CV(iptr, i, j, k);
	  }
	  LOOP(i,nx) { 
	    CV(iptr, i, j, k) = 0.0;
	    for(i1=-(nxaver-1)/2; i1<=(nxaver-1)/2; i1++) {
	      l = i+i1;
	      if (l<0 || l>=nx) continue;
	      CV(iptr, i, j, k) += row[l];
	    }
	    CV(iptr, i, j, k) /= nxaver;	    
	  }
	}
      }
      free(row);
      write_image(outstr, iptr);
    } else if (nxaver>1 || nyaver>1 || nzaver>1) {  /*  averaging, but retaining size */
        if (Qmoving) error("moving=t not implemented in this mode");
        dprintf(0,"Averaging map %d * %d * %d pixels; mapsize %d * %d * %d\n",
                   nxaver,nyaver,nzaver,nx,ny,nz);
        nx1 = nx/nxaver;  if (nx % nxaver) warning("X binning not even");
        ny1 = ny/nyaver;  if (ny % nyaver) warning("Y binning not even");
        nz1 = nz/nzaver;  if (nz % nzaver) warning("X binning not even");
	LOOP(k1,nz1) {
	  k = k1*nzaver;
	  LOOP(j1,ny1) {
            j = j1*nyaver;
	    LOOP(i1,nx1) {
	      i = i1*nxaver;
	      sum = 0.0;
	      LOOP(k0,nzaver) LOOP(j0,nyaver) LOOP(i0,nxaver) sum += CV(iptr, i+i0, j+j0, k+k0);
	      sum /= (real) (nxaver*nyaver*nzaver);
	      LOOP(k0,nzaver) LOOP(j0,nyaver) LOOP(i0,nxaver) CV(iptr, i+i0, j+j0, k+k0) = sum;
            }
	  }
	}
	if (!Qdummy) ax_shift(iptr);
        write_image(outstr, iptr);
    } else if (Qreorder) {            	/* reordering */
      warning("new reordering mode %s",reorder);
      if (streq(reorder,"xyz")) {
	create_cube(&iptr1,nx1,ny1,nz1);
	ax_copy(iptr,iptr1);
	LOOP(k,nz1) LOOP(j,ny1) LOOP(i,nx1) CV(iptr1,i,j,k) = CV(iptr,i,j,k);
      } else if (streq(reorder,"xzy")) {
	create_cube(&iptr1,nx1,nz1,ny1);
	ax_copy(iptr,iptr1);
	ax_swap_yz(iptr1);
	LOOP(k,nz1) LOOP(j,ny1) LOOP(i,nx1) CV(iptr1,i,k,j) = CV(iptr,i,j,k);
      } else if (streq(reorder,"yxz")) {
	create_cube(&iptr1,ny1,nx1,nz1);
	ax_copy(iptr,iptr1);
	ax_swap_xy(iptr1);
	LOOP(k,nz1) LOOP(j,ny1) LOOP(i,nx1) CV(iptr1,i,j,k) = CV(iptr,j,i,k);
      } else if (streq(reorder,"yzx")) {
	create_cube(&iptr1,ny1,nz1,nx1);
	ax_copy(iptr,iptr1);
	ax_swap_xy(iptr1);
	ax_swap_yz(iptr1);
	LOOP(k,nz1) LOOP(j,ny1) LOOP(i,nx1) CV(iptr1,i,j,k) = CV(iptr,j,k,i);
      } else if (streq(reorder,"zxy")) {
	create_cube(&iptr1,nz1,nx1,ny1);
	ax_copy(iptr,iptr1);
	ax_swap_xy(iptr1);
	ax_swap_xz(iptr1);
	LOOP(k,nz1) LOOP(j,ny1) LOOP(i,nx1) CV(iptr1,i,j,k) = CV(iptr,k,i,j);
      } else if (streq(reorder,"zyx")) {
	create_cube(&iptr1,nz1,ny1,nx1);
	ax_copy(iptr,iptr1);
	ax_swap_xz(iptr1);
	LOOP(k,nz1) LOOP(j,ny1) LOOP(i,nx1) CV(iptr1,i,j,k) = CV(iptr,k,j,i);
      }
      if (!Qdummy) ax_shift(iptr1);
      write_image(outstr, iptr1);
    } else if (Qsample) {            	/* straight sub-sampling */
      create_cube(&iptr1,nx1,ny1,nz1);
      LOOP(k,nz1)
	LOOP(j,ny1)
	  LOOP(i,nx1)
	    CV(iptr1,i,j,k) = CV(iptr,ix[i],iy[j],iz[k]);
      warning("Attempting to fix the WCS");

      Xmin(iptr1) = Xmin(iptr) + ix[0]*Dx(iptr);
      Dx(iptr1)   = (ix[1]-ix[0]) * Dx(iptr);

      Ymin(iptr1) = Ymin(iptr) + iy[0]*Dy(iptr);
      Dy(iptr1)   = (iy[1]-iy[0]) * Dy(iptr);

      Zmin(iptr1) = Zmin(iptr) + iz[0]*Dz(iptr);
      Dz(iptr1)   = (iz[1]-iz[0]) * Dz(iptr);
      dprintf(0,"WCS Corner: %g %g %g\n",Xmin(iptr1),Ymin(iptr1),Zmin(iptr1));

      if (!Qdummy) ax_shift(iptr1);
      write_image(outstr, iptr1);
    } else {                            /* nothing really done, still a great benchmark */
      warning("No x=,y=,z= selection applied");
      if (!Qdummy) ax_shift(iptr);
      write_image(outstr, iptr);
    }
}

/*
 * either initialize idx array, if not done, or normalize to 0..n-1
 */
 
int ax_index(string name, int n, int n1, int *idx)
{
    int i;
    
    if (n1==0) {		/* copy array */
        n1=n;
        for (i=0; i<n; i++) idx[i] = i;
    } else {
        for (i=0; i<n1; i++) {
            if (idx[i] < 1 || idx[i] > n)
                error("Index %d illegal in %s axis; max %d",idx[i],name,n);
            idx[i] -= 1;
        }
    }
    return n1;
}


void ax_copy(imageptr i0, imageptr i1)
{
  Dx(i1) = Dx(i0);
  Dy(i1) = Dy(i0);
  Dz(i1) = Dz(i0);
  Xmin(i1) = Xmin(i0);
  Ymin(i1) = Ymin(i0);
  Zmin(i1) = Zmin(i0);
  Xref(i1) = Xref(i0);
  Yref(i1) = Yref(i0);
  Zref(i1) = Zref(i0);
  Beamx(i1) = Beamx(i0);
  Beamy(i1) = Beamy(i0);
  Beamz(i1) = Beamz(i0);
  if (Namex(i0))  Namex(i1) = strdup(Namex(i0));
  if (Namey(i0))  Namey(i1) = strdup(Namey(i0));
  if (Namez(i0))  Namez(i1) = strdup(Namez(i0));
}

/* ax_shift:  
 *
 */

void ax_shift(imageptr iptr)
{
  if (Nx(iptr)==1 && Ny(iptr)>1 && Nz(iptr)==1) {            /* permute {x,y} */
    dprintf(0,"Permuting x,y axis for dummy=f\n");
    ax_swap_xy(iptr);
  } else if (Nx(iptr)>1 && Ny(iptr)==1 && Nz(iptr)>1) {      /* permute {y,z} */
    dprintf(0,"Permuting y,z axis for dummy=f\n");
    ax_swap_yz(iptr);
  } else if (Nx(iptr)==1) {                                  /* permute (x,z) */
    dprintf(0,"Permuting x,z axis for dummy=f\n");
    ax_swap_xz(iptr);
  }
}



void ax_swap_xy(imageptr iptr) 
{
  int    tmpi;
  real   tmpr;
  string tmps;

  //SWAP(Nx(iptr),Ny(iptr),tmpi);
  SWAP(Dx(iptr),Dy(iptr),tmpr);
  SWAP(Xmin(iptr),Ymin(iptr),tmpr);
  SWAP(Xref(iptr),Yref(iptr),tmpr);
  SWAP(Beamx(iptr),Beamy(iptr),tmpr);
  SWAP(Namex(iptr),Namey(iptr),tmps);
}

void ax_swap_xz(imageptr iptr) 
{
  int    tmpi;
  real   tmpr;
  string tmps;

  //SWAP(Nx(iptr),Nz(iptr),tmpi);
  SWAP(Dx(iptr),Dz(iptr),tmpr);
  SWAP(Xmin(iptr),Zmin(iptr),tmpr);
  SWAP(Xref(iptr),Zref(iptr),tmpr);
  SWAP(Beamx(iptr),Beamz(iptr),tmpr);
  SWAP(Namex(iptr),Namez(iptr),tmps);
}

void ax_swap_yz(imageptr iptr) 
{
  int    tmpi;
  real   tmpr;
  string tmps;

  //SWAP(Nz(iptr),Ny(iptr),tmpi);
  SWAP(Dz(iptr),Dy(iptr),tmpr);
  SWAP(Zmin(iptr),Ymin(iptr),tmpr);
  SWAP(Zref(iptr),Yref(iptr),tmpr);
  SWAP(Beamz(iptr),Beamy(iptr),tmpr);
  SWAP(Namez(iptr),Namey(iptr),tmps);
}
