
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <falcON_C.h>

#define Pi   3.14159265358979323846264338328
#define TPi  (2*Pi)
#define PiH  (0.5*Pi)

double ms_Force_gamma(double gamma)
{
  double x=5-3*gamma,f;
  if(x < 0.) return 0.;
  f =1/x;
  f-=4/(x+=1);
  f+=6/(x+=1);
  f-=4/(x+=1);
  f+=1/(x+=1);
  return (3-gamma)*f;
}

double ms_density_gamma(double gamma)
{
  double x=3-3*gamma,f;
  if(x < 0.) return 0.;
  f = 1/x;
  f-= 8/(x+=1);
  f+=28/(x+=1);
  f-=56/(x+=1);
  f+=70/(x+=1);
  f-=56/(x+=1);
  f+=28/(x+=1);
  f-= 8/(x+=1);
  f+= 1/(x+=1);
  return pow(3-gamma,3)*f/pow(4*Pi,2);
}

#if defined(SINGLE_DOUBLE) || defined(DOUBLE_DOUBLE)
#  define INPUT_TYPE double
#else
#  define INPUT_TYPE float
#endif

int main(int argc, char* argv[])
{
  clock_t     cpu0, cpu1;
  INPUT_TYPE  time_setup, time_grow, time_tree;
  INPUT_TYPE  rmax, GAM;
  int         MOD, C, K, N, Ncrit=6, p, *F, ig,Ng;
  long        S;
  INPUT_TYPE *M,*X[3],*A[3],*PH,*RH,*RHO,*AT[3], EPS,TH ;

  int         i,j;
  INPUT_TYPE  g1,g2,ig2,g3,ig3, tm,r,R,phi,fr,P,cth,mi,rh;
  INPUT_TYPE  ase,dacc,damx,tmp,sf[3];

  /* 
   * read parameters
   */
  if(argc < 6 || argc > 10) {
    printf("\"TestGravC MOD gamma N S EPS [Ng THE K Nc] \" with \n");
    printf(" MOD = 0/1/2/3/4: hom. sphere / Plummer / gamma-model / Kuzmin disk / hom. disk\n");
    printf(" gamma               : parameter of model (if any)\n");
    printf(" N                   : number of positions\n");
    printf(" S = long            : seed for RNG\n");
    printf(" EPS                 : max / fixed softening length\n");
    printf(" Ng (default   0)    : # additional grow()s\n");
    printf(" THE(default 0.6)    : opening angle, if THE<0, use theta=const\n");
    printf(" K  (default   1)    : P_n softening kernel (P0=Plummer)\n");
    printf(" Nc (default   6)    : don't split cells with N <= Nc bodies\n");
    exit(1);
  }

  cpu0 = clock();
  rmax = 1.e3;
  p    = 1;
  MOD  = atoi(argv[p++]);
  GAM  = atof(argv[p++]);
  N    = atoi(argv[p++]);
  S    = atoi(argv[p++]);
  EPS  = atof(argv[p++]);
  if(argc>p) Ng   = atoi(argv[p++]); else Ng=0;
  if(argc>p) TH   = atof(argv[p++]); else TH=0.6;
  if(argc>p) K    = atoi(argv[p++]); else K =1;
  if(argc>p) Ncrit= atoi(argv[p++]); else Ncrit=6;

  /*
   * set up positions
   */
  F     = (int*)  malloc(sizeof(int)*N);
  X[0]  = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  X[1]  = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  X[2]  = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  M     = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  A[0]  = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  A[1]  = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  A[2]  = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  AT[0] = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  AT[1] = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  AT[2] = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  PH    = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  RH    = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  RHO   = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  srand48(S);

  g1  = 1.-GAM;
  g2  = 2.-GAM;
  ig2 = 1./g2;
  g3  = 3-GAM;
  ig3 = 1./g3;
  mi  = 1./(INPUT_TYPE)(N);

  for(i=0; i<N; i++) {
    do {
      switch(MOD) {
      case 0:  /* homogeneous */
	r  = pow(drand48(),1./3.);
	P  = 0.5*(3*r*r-1);
	fr =-1.;
	rh = 0.75/Pi;
	break;
      case 1:  /* Plummer */
	tm = pow(drand48(),2./3.);
	tm/= 1.-tm;                    // r^2
	r  = sqrt(tm);
	P  =-1./sqrt(tm+1);
	fr = P*P*P;
	rh = 0.75/Pi*pow(-P,5);
	break;
      case 2:  /* gamma-model */
	tm = pow(drand48(),ig3);
	r  = tm / (1.-tm);
	P  = (GAM==2.)?  log(tm) : -ig2*(1-pow(tm,g2));
	fr =-pow(tm,g1) * (1-tm) * (1-tm) / r;
	rh = 0.25*g3/Pi * pow(1.-tm,4)/pow(tm,GAM);
	break;
      case 3:  /* Kuzmin disk */
	tm = drand48();
	r  = sqrt((2-tm)*tm)/(1.-tm);
	P  =-1./sqrt(r*r+1);
	fr = P*P*P;
	rh =-0.5/Pi*fr;
	break;
      case 4: default: /* homogeneous disk */
	tm = drand48();
	r  = sqrt(tm);
	P  = 2*(r*r/3.-2);
	fr =-4*r/3.;
	rh = 1.0/Pi;
	break;
      }
    } while(r>rmax);
    cth      = (MOD>=3) ? 0. : 2*drand48()-1;
    R        = r*sqrt(1.-cth*cth);
    phi      = TPi*drand48();
    F[i]     = 1;
    X[0][i]  = R * sin(phi);
    X[1][i]  = R * cos(phi);
    X[2][i]  = r * cth;
    M[i]     = mi;
    for(j=0; j<3; j++) AT[j][i] = fr*X[j][i]; /* "true" force   */
    RH[i] = rh;                               /* "true" density */
  }

  cpu1       = clock();
  time_setup = (cpu1 - cpu0) / (INPUT_TYPE)(CLOCKS_PER_SEC);
  cpu0       = cpu1;
  printf("\n time needed for set up of X_i: %f\n",time_setup);

  /*
   * computing forces
   */

  falcON_initialize(F,M,X[0],X[1],X[2],
#ifdef falcON_INDI
		    0,
#endif
		    A[0],A[1],A[2],PH,RHO,N,EPS,TH,K);
  for(ig=0; ig<=Ng; ++ig) {
    falcON_grow(Ncrit);
    cpu1      = clock();
    time_grow = (cpu1 - cpu0) / (INPUT_TYPE)(CLOCKS_PER_SEC);
    cpu0      = cpu1;
    printf(" time needed by tree growth:    %f\n",time_grow);
  }

  falcON_approx_grav();
  cpu1      = clock();
  time_tree = (cpu1 - cpu0) / (INPUT_TYPE)(CLOCKS_PER_SEC);
  cpu0      = cpu1;
  printf(" time needed by tree forces:    %f\n",time_tree);
  falcON_stats();
  printf("\n");


  /*
   * outputs
   */
  
  ase  = 0.;
  damx = 0.;
  for(j=0; j<3; j++) sf[j]=0.;
  for(i=0; i<N; i++) {
    dacc = 0.;
    for(j=0; j<3; j++) {
      tmp = (A[j][i]-AT[j][i]);
      dacc += tmp*tmp;
      sf[j]+= M[i]*A[j][i];
    }
    ase += M[i] * dacc;
    if(dacc > damx) damx = dacc;
/*       printf("%f %f %f %f %f %f\n", X[0][i],X[1][i],X[2][i], */
/* 	                            A[0][i],A[1][i],A[2][i]);   */
  }
  switch(MOD){
  case 1: ase *= 13.125; break;
  case 3: ase *= 17.5;   break;
  case 4: ase *= 1.125;  break;
  case 0: ase /= 0.6;    break;
  case 2: if(GAM<1.66) ase /= ms_Force_gamma(GAM); break;
  default:  break;
  }
  if(MOD==2 && GAM>=1.66) printf(" ASE(F)           = %13.10f\n",ase);
  else                    printf(" ASE(F)/<F^2>     = %13.10f\n",ase);
  printf(" max (dF)^2       = %13.10f\n",damx);
  printf(" Sum m_i acc_i    = %e, %e, %e\n",sf[0],sf[1],sf[2]);
  return 0;
}
