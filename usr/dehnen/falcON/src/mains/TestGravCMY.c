// -*- C -*-                                                                   |
//-----------------------------------------------------------------------------+
//                                                                             |
// TestGravC.cc                                                                |
//                                                                             |
// Copyright (C) 2000-2005 Walter Dehnen                                       |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <FAlCON_C.h>

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

static double iRMAX;

inline double RNG()
{
  return rand() * iRMAX;
}


#ifdef falcON_DOUBLE
#  define INPUT_TYPE double
#else
#  define INPUT_TYPE float
#endif

typedef INPUT_TYPE INPUT_VECT[3];

int main(int argc, char* argv[])
{
  clock_t     cpu0, cpu1;
  INPUT_TYPE  time_setup, time_grow, time_tree;
  INPUT_TYPE  rmax, GAM;
  int         MOD, C, kern, N, Ncrit=6, p, *F, ig,Ngrow;
  long        S;
  INPUT_TYPE *M,*PH,*RH,*RHO,EPS,theta;
  INPUT_VECT *X,*A,*AT;

  int         i,j;
  INPUT_TYPE  g1,g2,ig2,g3,ig3, tm,r,R,phi,fr,P,cth,mi,rh;
  INPUT_TYPE  ase,dacc,damx,tmp;
  double      sf[3];

  /* 
   * read parameters
   */
  if(argc < 6 || argc > 11) {
    printf("\"TestGravC MOD gamma N S EPS [Ng THE K Nc Rmax] \" with \n");
    printf(" MOD = 0/1/2/3/4: "
	   "hom. sphere / Plummer / gamma-model / Kuzmin disk / hom. disk\n");
    printf(" gamma               : parameter of model (if any)\n");
    printf(" N                   : number of positions\n");
    printf(" S = long            : seed for RNG\n");
    printf(" EPS                 : max / fixed softening length\n");
    printf(" Ng (default    0)   : # additional grow()s\n");
    printf(" THE(default %4.2f)   : opening angle, if THE<0, use theta=const\n",
	   falcON_default_theta());
    printf(" K  (default    %d)   : P_n softening kernel (P0=Plummer)\n",
	   falcON_default_kernel());
    printf(" Nc (default   %2d)   : don't split cells with N <= Nc bodies\n",
	   falcON_default_Ncrit());
    printf(" Rmax (default 1e3)   : max radius\n");
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
  if(argc>p) Ngrow= atoi(argv[p++]); else Ngrow=0;
  if(argc>p) theta= atof(argv[p++]); else theta=falcON_default_theta();
  if(argc>p) kern = atoi(argv[p++]); else kern =falcON_default_kernel();
  if(argc>p) Ncrit= atoi(argv[p++]); else Ncrit=falcON_default_Ncrit();
  if(argc>p) rmax = atoi(argv[p++]);

  /*
   * set up positions
   */
  F   = (int       *)malloc(sizeof(int)       *N);
  X   = (INPUT_VECT*)malloc(sizeof(INPUT_VECT)*N);
  M   = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  A   = (INPUT_VECT*)malloc(sizeof(INPUT_VECT)*N);
  AT  = (INPUT_VECT*)malloc(sizeof(INPUT_VECT)*N);
  PH  = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  RH  = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);
  RHO = (INPUT_TYPE*)malloc(sizeof(INPUT_TYPE)*N);

  srand(S);
  iRMAX = 1.0 / (1.0 + RAND_MAX);

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
	r  = pow(RNG(),1./3.);
	P  = 0.5*(3*r*r-1);
	fr =-1.;
	rh = 0.75/Pi;
	break;
      case 1:  /* Plummer */
	tm = pow(RNG(),2./3.);
	tm/= 1.-tm;                    // r^2
	r  = sqrt(tm);
	P  =-1./sqrt(tm+1);
	fr = P*P*P;
	rh = 0.75/Pi*pow(-P,5);
	break;
      case 2:  /* gamma-model */
	tm = pow(RNG(),ig3);
	r  = tm / (1.-tm);
	P  = (GAM==2.)?  log(tm) : -ig2*(1-pow(tm,g2));
	fr =-pow(tm,g1) * (1-tm) * (1-tm) / r;
	rh = 0.25*g3/Pi * pow(1.-tm,4)/pow(tm,GAM);
	break;
      case 3:  /* Kuzmin disk */
	tm = RNG();
	r  = sqrt((2-tm)*tm)/(1.-tm);
	P  =-1./sqrt(r*r+1);
	fr = P*P*P;
	rh =-0.5/Pi*fr;
	break;
      case 4: default: /* homogeneous disk */
	tm = RNG();
	r  = sqrt(tm);
	P  = 2*(r*r/3.-2);
	fr =-4*r/3.;
	rh = 1.0/Pi;
	break;
      }
    } while(rmax>0 && r>rmax);
    cth      = (MOD>=3) ? 0. : 2*RNG()-1;
    R        = r*sqrt(1.-cth*cth);
    phi      = TPi*RNG();
    F[i]     = 1;
    X[i][0]  = R * sin(phi);
    X[i][1]  = R * cos(phi);
    X[i][2]  = r * cth;
    M[i]     = mi;
    for(j=0; j<3; j++) 
      AT[i][j] = fr*X[i][j];                  /* "true" force   */
    RH[i] = rh;                               /* "true" density */
  }

  cpu1       = clock();
  time_setup = (cpu1 - cpu0) / (INPUT_TYPE)(CLOCKS_PER_SEC);
  cpu0       = cpu1;
  printf("\n time needed for set up of X_i: %f\n",time_setup);

  /*
   * computing forces
   */

/*   falcON_initialize(F,M,(INPUT_TYPE*)X, */
/* #ifdef falcON_INDI */
/* 		    0, */
/* #endif */
/* 		    (INPUT_TYPE*)A,PH,RHO,N,EPS,theta,kern,1); */
  falcON_initialize(F,M,(INPUT_TYPE*)X,
#ifdef falcON_INDI
		    0,
#endif
		    (INPUT_TYPE*)A,PH,RHO,N,0,EPS,theta,kern,1);
  for(ig=0; ig<=Ngrow; ++ig) {
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
      tmp = (A[i][j]-AT[i][j]);
      dacc += tmp*tmp;
      sf[j]+= M[i]*A[i][j];
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
