
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <falcON_C.h>

#define twothird 0.666666666666666666666666
#define Pi       3.14159265358979323846264338328
#define TPi      (2*Pi)
#define PiH      (0.5*Pi)

#if defined(SINGLE_DOUBLE) || defined(DOUBLE_DOUBLE)
#  define INPUT_TYPE double
#else
#  define INPUT_TYPE float
#endif

#ifndef NDIM
#  define NDIM 3
#endif

void dmintmin(const INPUT_TYPE R[NDIM], const INPUT_TYPE V[NDIM],
	      const INPUT_TYPE T, INPUT_TYPE *dmin, INPUT_TYPE *tmin)
{
  int  i;
  INPUT_TYPE RV,Vq,tmp;
  RV = 0.;
  Vq = 0.;
  for(i=0; i<NDIM; i++) {
    RV += R[i]*V[i];
    Vq += V[i]*V[i];
  }
  *tmin = -RV/Vq;
  if(*tmin < 0) *tmin = 0;
  if(*tmin > T) *tmin = T;
  RV = 0.;
  for(i=0; i<NDIM; i++) {
    tmp = R[i]+ *tmin * V[i];
    RV += tmp*tmp;
  }
  *dmin = sqrt(RV);
}

int main(int argc, char* argv[])
{
  clock_t     cpu0, cpu1;
  INPUT_TYPE  time_prep, time_tree, time_stic;
  int         N=0,Ncrit=6,NS,p,*F,Z[NDIM]={0},*I1,*I2,NA=0;
  long        Seed;
  INPUT_TYPE *M,*S,*X[NDIM],*V[NDIM],TAU,SIZE,dX[NDIM],dV[NDIM];
  int         i,j;
  INPUT_TYPE  mass=1./(INPUT_TYPE)(N),tmin,dmin;
  INPUT_TYPE  Mr,r,rq,P,Pot,veq,f0,v,vq,fs,R,cth,phi;
  FILE       *fp_out;

  if(argc < 6 || argc > 7) {
    printf("\n Testing tree::sticky_interactions with a Plummer sphere\n\n");
    printf(" \"TestPair N S C Ns s tau [Nc]\" with \n");
    printf(" N              : number of particles\n");
    printf(" S = long int   : seed for RNG\n");
    printf(" Ns             : first Nsub bodies to make sticky\n");
    printf(" s              : size of sticky particles\n");
    printf(" tau            : time step; <0 -> SPH\n");
    printf(" Nc (default 6) : don't split cells with N <= Nc bodies\n");
    exit(1);
  }

  cpu0 = clock();
  p    = 1;
  N    = atoi(argv[p++]);
  Seed = atoi(argv[p++]);
  NS   = atoi(argv[p++]);
  SIZE = atof(argv[p++]);
  TAU  = atof(argv[p++]);
  if(argc > p) Ncrit= atoi(argv[p++]);

  if(NS>N) NS=N;
  srand48(Seed);

  F     = (int*)  malloc(sizeof(int)*N);
  I1    = (int*)  malloc(sizeof(int)*N);
  I2    = (int*)  malloc(sizeof(int)*N);
  for(j=0; j<NDIM; j++) {
    X[j]  = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
    V[j]  = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  }
  S     = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  M     = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);

  for(i=0; i<N; i++) {
    Mr  = pow(drand48(),twothird);
    rq  = Mr/(1-Mr);
    Pot =-1./sqrt(rq+1.);
    r   = sqrt(rq);
    veq =-2*Pot;
    do {
      P   = veq*pow(drand48(),1./4.5);
      vq  = veq-P;
      v   = sqrt(vq);
      f0  = v*pow(P,3.5);
      fs  = veq*f0;
      f0 *= v;
    } while (fs* drand48()>f0);
#if NDIM == 3
    cth      = 2*drand48()-1.;
    phi      = TPi*drand48();
    R        = r * sqrt(1.-cth*cth);
    X[0][i]  = R * sin(phi);
    X[1][i]  = R * cos(phi);
    X[2][i]  = r * cth;
    cth      = 2*drand48()-1.;
    phi      = TPi*drand48();
    R        = v * sqrt(1.-cth*cth);
    V[0][i]  = R * sin(phi);
    V[1][i]  = R * cos(phi);
    V[2][i]  = v * cth;
#else
    phi      = TPi*drand48();
    X[0][i]  = r * sin(phi);
    X[1][i]  = r * cos(phi);
    phi      = TPi*drand48();
    V[0][i]  = v * sin(phi);
    V[1][i]  = v * cos(phi);
#endif
    M[i]     = mass;
    S[i]     = SIZE;
    F[i]     = (i<NS)? 13 : 1;
  }
  cpu1 = clock();
  time_prep = (cpu1 - cpu0) / (INPUT_TYPE)(CLOCKS_PER_SEC);
  cpu0 = cpu1;
  printf(" time needed for set up of (X,V)_i: %f\n",time_prep);

  falcON_initialize(F,M,X[0],X[1],X[2],0,0,0,0,0,0,N,0,0.6,1);
  falcON_grow(Ncrit);

  cpu1 = clock();
  time_tree = (cpu1 - cpu0) / (INPUT_TYPE)(CLOCKS_PER_SEC);
  cpu0 = cpu1;
  printf(" time needed for growth of tree: %f\n",time_tree);

  falcON_iactionlist(I1,I2,N,&NA,S,TAU,V[0],V[1],V[2]);

  cpu1 = clock();
  time_stic = (cpu1 - cpu0) / (INPUT_TYPE)(CLOCKS_PER_SEC);
  cpu0 = cpu1;
  printf(" time needed for pair finding: %f\n",time_stic);

  if(TAU<0) {
    fp_out = fopen("pairs.C","w");
    for(i=0; i<N && i<NA; i++) {
      rq = 0;
      for(j=0; j<NDIM; j++) {
	r   = X[j][I1[i]] - X[j][I2[i]];
	rq += r*r;
      }
      r = sqrt(rq);
      fprintf(fp_out," %d  %d  %f\n",I1[i],I2[i],r);
    }
  } else {
    fp_out = fopen("pairs.C","w");
    for(i=0; i<N && i<NA; i++) {
      for(j=0; j<NDIM; j++) {
	dX[j] = X[j][I1[i]] - X[j][I2[i]];
	dV[j] = V[j][I1[i]] - V[j][I2[i]];
      }
      dmintmin(dX,dV,TAU,&dmin,&tmin);
      fprintf(fp_out," %d  %d  %f  %f \n",I1[i],I2[i],tmin,dmin);
    }
  }
  if(NA < N) printf(" %d pairs found; %d written to file \"pairs.C\"\n",NA,N);
  else       printf(" %d pairs found and written to file \"pairs.C\"\n",NA);
  
  return 0;
}
