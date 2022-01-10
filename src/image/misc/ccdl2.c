/* 
 * CCDL2: Find the Lagrangian points in a 2D potential
 *
 *	quick and dirty:  9-jan-2022	pjt
 *
 * @todo   it's off by 1/2 pixel ... hmmm
 */



#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <image.h>
#include <moment.h>

string defv[] = {
  "in=???\n       Input image file",
  "VERSION=0.2\n  9-jan-2022 PJT",
  NULL,
};

string usage = "find Lagrangian points in a 2D potential";

string cvsid="$Id$";

#define CV(i,j) MapValue(iptr,i,j)

void nemo_main()
{
  stream   instr;
  imageptr iptr=NULL;
  int      i,j,nx,ny,nz;
  int      sx1,sx2,sy1,sy2;
  real     x,y;

  warning("Draft program");

  instr = stropen(getparam("in"), "r");
  
  read_image( instr, &iptr);
  nx = Nx(iptr);  ny = Ny(iptr);  nz = Nz(iptr);
  if (nz > 1) error("Cannot handle cubes");
  dprintf(0,"Dx,Dy= %g %g\n", Dx(iptr), Dy(iptr));
  printf("###: ix  iy     x   y   Potential\n");

  for (j=1; j<ny-1; j++) {
    y = Ymin(iptr) + (j-Yref(iptr)) * Dy(iptr);
    for (i=1; i<nx-1; i++) {
      x = Xmin(iptr) + (i-Xref(iptr)) * Dx(iptr);
      sx1 = SGN(CV(i+1,j)-CV(i,j));
      sx2 = SGN(CV(i-1,j)-CV(i,j));
      sy1 = SGN(CV(i,j+1)-CV(i,j));
      sy2 = SGN(CV(i,j-1)-CV(i,j));
      if (sx1==1 && sx2==1 && sy1==1 && sy2==1)
	printf("L11: %d %d   %g %g  %g  %g %g %g %g\n", i, j, x, y, CV(i,j),
	       (CV(i+1,j)-CV(i,j)),
	       (CV(i-1,j)-CV(i,j)),
	       (CV(i,j+1)-CV(i,j)),
	       (CV(i,j-1)-CV(i,j)));
      else if (sx1==-1 && sx2==-1 && sy1==-1 && sy2==-1)
	printf("L00: %d %d   %g %g  %g  %g %g %g %g\n", i, j, x, y, CV(i,j),
	       (CV(i+1,j)-CV(i,j)),
	       (CV(i-1,j)-CV(i,j)),
	       (CV(i,j+1)-CV(i,j)),
	       (CV(i,j-1)-CV(i,j)));
      else if (sx1==1 && sx2==1 && sy1==-1 && sy2==-1)
	printf("L10: %d %d   %g %g  %g  %g %g %g %g\n", i, j, x, y, CV(i,j),
	       (CV(i+1,j)-CV(i,j)),
	       (CV(i-1,j)-CV(i,j)),
	       (CV(i,j+1)-CV(i,j)),
	       (CV(i,j-1)-CV(i,j)));
      else if (sx1==-1 && sx2==-1 && sy1==1 && sy2==1)
	printf("L01: %d %d   %g %g  %g  %g %g %g %g\n", i, j, x, y, CV(i,j),
	       (CV(i+1,j)-CV(i,j)),
	       (CV(i-1,j)-CV(i,j)),
	       (CV(i,j+1)-CV(i,j)),
	       (CV(i,j-1)-CV(i,j)));      
    }
  }
}

