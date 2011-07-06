/*
 * AHF2GRP:  convert AHF to halo's group list
 *
 *       1-jul-11  V0.1 created, just for fun	PJT
 *
 * 
 */

#include <stdinc.h>
#include <getparam.h>

string defv[] = {
  "in=???\n           Dataset",
  "out=???\n          Output",
  "nbody=30\n        Maximum index particles in input (e.g. 128**3)",
  "nmax=10000\n       Max lines in data, if pipe",
  "VERSION=0.2\n      6-jul-11 PJT",
  NULL,
};

string usage="AHF to halo group";

void convert(int nbody, int n,int *idx);

static int *ihalo = NULL;

nemo_main()
{
  int i, npt, npt1, nmax = getiparam("nmax");
  int col1 = 1;
  int *coldat1[2];
  int colnr1[2];
  real *x1 = NULL;
  int *idx1 = NULL;
  string input1 = getparam("in");
  string input2 = getparam("out");
  stream instr1, instr2;
  int nbody = getiparam("nbody");

  
  instr1 = stropen(input1,"r");
  instr2 = stropen(input2,"w");
  npt1 = nemo_file_lines(input1,nmax);
  idx1 = (int *) allocate(npt1 * sizeof(int));

  dprintf(0,"nbody=%d\n",nbody);
  ihalo =  (int *) allocate(nbody*sizeof(int));

  coldat1[0] = idx1;    
  colnr1[0] = col1;
  npt1 = get_itable(instr1,1,colnr1,coldat1,npt1);
  if (npt1 < 0)
    error("%d: Could not read all data from %s",npt1, input1);
  strclose(instr1);

  convert(nbody, npt1,idx1);

  fprintf(instr2,"%d\n",nbody);
  for (i=0; i<nbody; i++)
    fprintf(instr2,"%d\n",ihalo[i]);
  strclose(instr2);

}


/*
 *  Nh, {Np_h1, i_1, i_2, .... i_Np_h1}, {Np_h2, ....}, .... {}   }
 * 
 *  to
 *  Nbody, i_H, .....
 */

void convert(int nbody, int ndata,int *kdata)
{
  int i,j,k,n,nh,np,p, nuse = 0;
  
  for (i=0; i<nbody; i++)      /* for all stars */
    ihalo[i] = 0;              /* initialize halo membership */

  k = 0;                       /* counts lines in kdata[] */
  nh = kdata[k++];             /* number of halos */
  dprintf(1,"Number of halos: %d\n",nh);
  for (j=0; j<nh; j++)   {     /* loop over the halos */
    np = kdata[k++];           /* number of particles in this halo */
    dprintf(2,"Number of stars in halo %d: %d\n",j+1,np);
    for (i=0; i<np; i++) {     /* loop over particles in this halo */
      p = kdata[k++];          /* particle number */
      if (p<0 || p>nbody) error("bad particle %d on line %d\n",p,k);
      if (k>ndata) error("file incomplete for k=%d ndata=%d\n",k,ndata);
      ihalo[p-1] = j+1;        /* assign halo group # to particle */
      nuse++;                  /* count stars in halos */
    }
  }
  dprintf(0,"%d/%d stars in halos\n",  nuse, nbody);
  
}
