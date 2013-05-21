/* 
 * CCDMASK: mask one or more maps into a bitmask
 *
 *     21-may-2013    Initial version, cloned off ccdmath      PJT
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
  "VERSION=0.1\n   21-may-2013 PJT",
  NULL,
};

string usage = "image masking";

string cvsid="$Id$";

#ifndef HUGE
# define HUGE 1.0e20
#endif

#define MAXIMAGE 20

imageptr iptr;	/* pointers to (input) images */
imageptr optr;	/* pointers to (input) images */

int      nimage;                /* actual number of input images */

#define MAXNAX 3

double crval[MAXNAX], crpix[MAXNAX], cdelt[MAXNAX];
int nwcs = 0;



extern  int debug_level;		/* see initparam() */
extern    void    dmpfien();
extern    int     inifien();
extern    void    dofien();

extern string *burststring(string,string);


void nemo_main ()
{
  string *fnames, fie;
  stream  instr;                /* input files */
  stream  outstr;                     /* output file */
  int     i,j,k,l,nx,ny,nz;     /* size of scratch map */
  int     noper;                      /* number of images needed from oper */
  int     nclip;
  real    cv, clip[MAXIMAGE];
  bool    Qclip = hasvalue("clip");
  int     imask,nmask;

  fnames = burststring(getparam("in"), ", ");  /* input file names */
  nimage = xstrlen(fnames, sizeof(string))-1;  /* number of files */

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
  dprintf(0,"%d/%d masked\n",nmask,nx*ny*nz);
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
    dprintf(0,"%d/%d masked\n",nmask,nx*ny*nz);
    strclose(instr);    
  }
  
  write_image (outstr,optr);         /* write image to file */
  strclose(outstr);
}
