/*
 *  CCDREORDER - re-order the axes of a cube.
 *
 *       Currently just a benchmark program
 */

#include <nemo.h>
#include <mdarray.h>

string defv[] = {
  "in=\n               Input file (ignored for now)",
  "out=\n              Output file  (ignored for now)", 
  "order=2,1,0\n       New order of output cube (0=first dimension)",
  "dims=100,100,100\n  Dimensions of the cube",
  "seed=123\n          Random seed for [0,1] in the created cube",
  "iter=1\n            Times to repeat the reorder",
  "VERSION=0.1\n       26-Dec-2019 PJT",
  NULL,
};

string usage="benchmark reordering (with optional openmp)";

string cvsid="$Id:$";


void reorder3(int ndims, int *dims,  int *doms, int *order, int iter);

void nemo_main()
{
  string fin = getparam("in");
  string fout = getparam("out");
  int seed = init_xrandom(getparam("seed"));
  int iter = getiparam("iter");
  int order[MDMAXDIM], dims[MDMAXDIM], doms[MDMAXDIM];
  int norder = nemoinpi(getparam("order"), order, MDMAXDIM);
  int ndims  = nemoinpi(getparam("dims"),  dims,  MDMAXDIM);
  int i;

  if (hasvalue("in"))  warning("benchmark mode, no in=");
  if (hasvalue("out")) warning("benchmark mode, no out=");

  if (ndims != norder) error("dims= and order= need the same number of arguments");
  for (i=0; i<ndims; i++) {
    doms[i] = dims[order[i]];
    dprintf(0,"dims[%d] = %d   doms[%d] = %d\n",i,dims[i],i,doms[i]);
  }
  
  reorder3(ndims, dims, doms, order, iter);
}

/*
 * First a dumb approach. It's worth looking at MIRIAD's trnio.for routine
 * which uses a pattern and scratch array
 * These assume things won't fit in memory. We do want to try out memory
 * based models here
 */

//            order   iter   
//  10^2 cube 0,1,2  *1000  =  1.9"
//  10^2 cube 2,1,0  *1000  =  4.5"
//  10^3 cube 0,1,2  *1     = 14"
//  10^3 cube 2,1,0  *1     = 60"


void reorder3(int ndims, int *dims,  int *doms, int *order, int iter)
{
  int i0, i1, i2, idx[3];
  if (ndims != 3) error("ndims=%d illegal, must be 3",ndims);
#if 0
  // something wrong here
  real x[dims[2]][dims[1]][dims[0]];
  real y[doms[2]][doms[1]][doms[0]];
#else
  mdarray3 x = allocate_mdarray3(dims[2],dims[1],dims[0]);
  mdarray3 y = allocate_mdarray3(doms[2],doms[1],doms[0]);
#endif  

  for (i2=0; i2 < dims[2]; ++i2) 
    for (i1=0; i1 < dims[1]; ++i1)
      for (i0=0; i0 < dims[0]; ++i0) 
	x[i2][i1][i0] = xrandom(0.0,1.0);

  while (iter-- > 0) {
#if 1
    for (i2=0; i2 < dims[2]; ++i2) {
      idx[2] = i2;
      for (i1=0; i1 < dims[1]; ++i1){
	idx[1] = i1;
	for (i0=0; i0 < dims[0]; ++i0) {
	  idx[0] = i0;
	  y[idx[order[2]]][idx[order[1]]][idx[order[0]]] = x[i2][i1][i0];
	  // y[i2][i1][i0] = x[i2][i1][i0];	  	    
	}
      }
    }
#else    
    for (i2=0; i2 < dims[2]; ++i2) 
      for (i1=0; i1 < dims[1]; ++i1)
	for (i0=0; i0 < dims[0]; ++i0) 
	  y[i2][i1][i0] = x[i2][i1][i0];
    
#endif    
  } // while
}
