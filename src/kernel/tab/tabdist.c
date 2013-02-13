/*
 *  Compare distance between points in two tables
 *
 *      13-feb-2013     0.1    Created, Q&D
 */

#include <stdinc.h>
#include <getparam.h>

string defv[] = {
    "in1=???\n		First dataset",
    "in2=???\n		Second dataset",
    "col1=1,2\n         X and Y Column from 1st dataset",
    "col2=1,2\n         X and Y Column from 2nd dataset",
    "VERSION=0.1\n	13-feb-2013 PJT",
    NULL,
};

string usage="Compare distances between points in two tables";
string cvsid="$Id$";

void compare_1d(int npt1, int npt2, real *x1, real *x2);

inline real dist1(real x1,real x2) {
  real d = x1-x2;
  if (d > 0.0) return d;
  return -d;
}

inline real dist2(real x1,real y1,real x2,real y2) {
  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
}

inline real dist3(real x1,real y1,real z1,real x2,real y2,real z2) {
  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
}


nemo_main()
{
  int i, npt, npt1, npt2;
  int col1 = getiparam("col1");
  int col2 = getiparam("col2");
  real *coldat1[2], *coldat2[2];
  int colnr1[2], colnr2[2];
  real *x1 = NULL, *x2 = NULL;
  real t, f, tu, tp, tprob, fprob, tuprob, tpprob, ks2, ks2prob;
  string input1 = getparam("in1");
  string input2 = getparam("in2");
  stream instr1, instr2;
  
  
  instr1 = stropen(input1,"r");
  instr2 = stropen(input2,"r");

  npt1 = nemo_file_lines(input1,0);
  npt2 = nemo_file_lines(input2,0);

  x1 = (real *) allocate(npt1 * sizeof(real));
  x2 = (real *) allocate(npt2 * sizeof(real));

  
  coldat1[0] = x1;    colnr1[0] = col1;
  coldat2[0] = x2;    colnr2[0] = col2;

  npt1 = get_atable(instr1,1,colnr1,coldat1,npt1);
  npt2 = get_atable(instr2,1,colnr2,coldat2,npt2);

  dprintf(0,"Npt1=%d Npt2=%d\n",npt1,npt2);
  
  strclose(instr1);
  strclose(instr2);

  compare_1d(npt1, npt2, x1, x2);

}




void compare_1d(int npt1, int npt2, real *x1, real *x2)
{
  int i, j, jmin;
  real d, dmin;
  
  for (i=0; i<npt1; i++) {
    dmin = dist1(x1[i],x2[0]);
    jmin = 0;
    for (j=1; j<npt2; j++) {
      d = dist1(x1[i],x2[j]);
      if (d<dmin) {
	dmin = d;
	jmin = j;
      }
    }
    printf("%d %g   %d %g  %g\n",i+1,x1[i],jmin+1,x2[jmin], dmin);
  }
  
}
