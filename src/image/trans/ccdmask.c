/* 
 * CCDMASK: mask one or more maps into a bitmask
 *
 *     21-may-2013    Initial version, cloned off ccdmath      PJT
 *     29             report final mask count
 *                      
 */

#include <stdinc.h>
#include <ctype.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <strlib.h>
#include <image.h>

string defv[] = {
  "in=???\n        Input file(s), separated by comma's (optional)",
  "out=???\n       Output file",
  "clip=0.0\n      Clip level(s) for each input file",
  "VERSION=0.2\n   29-may-2013 PJT",
  NULL,
};

string usage = "image masking using clip levels";

string cvsid="$Id$";

#ifndef HUGE
# define HUGE 1.0e20
#endif

#define MAXIMAGE 64

extern string *burststring(string,string);

imageptr iptr, optr;              /* images */


void nemo_main ()
{
  string  *fnames;
  stream   instr, outstr;                /* files */
  int      i,j,k,l,nx,ny,nz;    
  int      nclip, nimage;
  real     cv, clip[MAXIMAGE];
  bool     Qclip = hasvalue("clip");
  int      imask,nmask;

  fnames = burststring(getparam("in"), ", ");  /* input file names */
  nimage = xstrlen(fnames, sizeof(string))-1;  /* number of files */
  if (nimage > MAXIMAGE) error("cannot handle %d images, MAXIMAGE=%d",nimage,MAXIMAGE);

  /* read clip levels, one per file. if not enough, last one repeats */
  nclip = nemoinpr(getparam("clip"),clip,MAXIMAGE);
  if (nclip<1) error("error parsing clip=%s",getparam("clip"));
  for (l=nclip; l<nimage; l++) clip[l] = clip[l-1];


  outstr = stropen (getparam("out"),"w");  /* open output file first ... */

  dprintf(0,"%d input file(s)\n",nimage);

  instr = stropen(fnames[0],"r");    /* open first file */
  read_image (instr, &optr);  
  nx = Nx(optr);	
  ny = Ny(optr);
  nz = Nz(optr);
  
  nmask = 0;
  imask = 1;
  for (k=0; k<nz; k++) {
    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) { 
	cv = CubeValue(optr,i,j,k);
	if (cv > clip[0]) {
	  CubeValue(optr,i,j,k) = imask;
	  nmask++;
	} else
	  CubeValue(optr,i,j,k) = 0;
      }
    }
  }
  dprintf(0,"%d/%d masked in map 1 @ clip %g\n",nmask,nx*ny*nz,clip[0]);
  strclose(instr);

  for (l=1; l<nimage; l++) {
    instr = stropen(fnames[l],"r");    /* open first file */
    read_image (instr, &iptr);  
    imask *= 2;
    nmask = 0;
    for (k=0; k<nz; k++) {
      for (j=0; j<ny; j++) {
	for (i=0; i<nx; i++) { 
	  cv = CubeValue(iptr,i,j,k);
	  if (cv > clip[0]) {
	    CubeValue(optr,i,j,k) += imask;
	    nmask++;
	  }
	}
      }
    }
    dprintf(0,"%d/%d masked in map %d @ clip %g\n",nmask,nx*ny*nz,l+1,clip[l]);
    strclose(instr);    
  }

  nmask = 0;
  for (k=0; k<nz; k++)
      for (j=0; j<ny; j++)
	for (i=0; i<nx; i++) 
	  if (CubeValue(iptr,i,j,k) > 0) nmask++;
  dprintf(0,"%d/%d masked in final map\n",nmask,nx*ny*nz);
  
  write_image (outstr,optr);         /* write image to file */
  strclose(outstr);
}
