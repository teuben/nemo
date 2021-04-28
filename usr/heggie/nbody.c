#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#define NMAX 1024
#define NDIM 3
#define loop(idx,last) for (idx = 0; idx < last ; idx++)
#define min(x, y) (x<y?x:y)
typedef double real;
static int i, j, k;

void main(int argc, char **argv) {
    real x[NMAX][NDIM], xdot[NMAX][NDIM], f[NMAX][NDIM],
      /**/fdot[NMAX][NDIM];
    real step[NMAX], tlast[NMAX], m[NMAX], t, dt;
    long int j, i, n, nsteps, noutp, cpu, timenow;
    extern void initialise(), runtime_output(), final_output();
/*     cpu = clock(); */
    time(&cpu);
    srand(cpu);
    assert(argc==4);
    n = atoi(argv[1]); 
    noutp = atoi(argv[2]); 
    dt = atof(argv[3])/(real)noutp; 
    t  = 0;
    /*    printf("\nInitial conditions.\n");*/
    initialise(n, x, xdot, f, fdot, step, tlast, m);
    nsteps=0;
    loop(j, noutp) {
      nsteps += hermite(x, xdot, f, fdot, step, tlast, m, &t, 
			/**/t+dt, n);
      runtime_output(x,t,n);
    }
    final_output(t,cpu,nsteps,noutp);
}

void initialise(int n, real x[NMAX][NDIM], 
	     /**/real xdot[NMAX][NDIM], real f[NMAX][NDIM], 
	     /**/real fdot[NMAX][NDIM], real step[NMAX], 
	     /**/real tlast[NMAX], real m[]) {
      real fi[NDIM],fidot[NDIM], radius, fmod, fdotmod;
      extern void uniform(), ffdot();
      radius = 6./5.;
      uniform(n,radius,x);
      radius = sqrt(5./6.);
      uniform(n,radius,xdot);
      loop(i, n) { m[i] = 1./(real)n;}
      loop(i, n) {
         tlast[i] = 0;
         ffdot(x,xdot,m,i,n,fi,fidot);
         loop(k, NDIM) {
            f[i][k] = fi[k];
            fdot[i][k] = fidot[k];
	 }
         fmod = fdotmod = 0;
         loop(k, NDIM) {
            fmod += pow(fi[k],2);
            fdotmod += pow(fidot[k],2);
	 }
         step[i] = 0.01*sqrt(fmod/fdotmod);
      }
}

void uniform(int n, real a, real x[NMAX][NDIM]) {
  real r, cos_theta, sin_theta, phi; 
  loop(i, n) {
    r = a*pow((real)rand()/((real)RAND_MAX+1), 1./3.); 
/*     r = a*pow(drand48(), 1./3.); */
    cos_theta = 2.*rand()/((real)RAND_MAX+1)-1;
    sin_theta = sqrt(1-pow(cos_theta, 2));
    phi = 2*M_PI*rand()/((real)RAND_MAX+1);
    x[i][0] = r*sin_theta*cos(phi);
    x[i][1] = r*sin_theta*sin(phi);
    x[i][2] = r*cos_theta;
  }
}



int hermite(real x[NMAX][NDIM], real xdot[NMAX][NDIM], 
	    /**/real f[NMAX][NDIM], real fdot[NMAX][NDIM], 
	    /**/real step[], real tlast[], real m[], real *t, 
	    /**/real tend, int n) 
{
  real xtemp[NMAX][NDIM], xdottemp[NMAX][NDIM], fi[NDIM], 
    /**/fidot[NDIM];
  real a2[NDIM], a3[NDIM], modfi, modfidot, moda2, moda3;
  real tmin, dt, dts, dtc, dnmtr, temp; 
  extern void ffdot();
  int imin, nsteps;
  assert(n<=NMAX);
  nsteps=0;
  do {
      tmin = step[0] + tlast[0];
      imin = 0;
      loop(i, n) {
	if (step[i] + tlast[i] < tmin) {
            tmin = step[i] + tlast[i];
            imin = i;
	}
      }
      assert(imin>=0);
      *t = tmin;
      i = imin;
      loop(j, n) {
         dt = *t - tlast[j];
         dts = pow(dt, 2);
         dtc = dts*dt;
         loop(k, NDIM) {
            xtemp[j][k] = x[j][k] + xdot[j][k]*dt + f[j][k]*dts/2  
                        + fdot[j][k]*dtc/6;
            xdottemp[j][k] = xdot[j][k] + f[j][k]*dt + 
	      /**/fdot[j][k]*dts/2;
	 }
      }
      ffdot(xtemp,xdottemp,m,i,n,fi,fidot);
      modfi = modfidot = moda2 = moda3 = 0;
      assert(step[i]!=0);

      loop(k, NDIM) {
         a2[k] = (-6*(f[i][k] - fi[k]) - step[i]*(4*fdot[i][k]  
               + 2*fidot[k]))/pow(step[i],2);
         a3[k] = (12*(f[i][k] - fi[k]) + 6*step[i]*(fdot[i][k] 
	       + fidot[k]))/pow(step[i],3);
         x[i][k] = xtemp[i][k] + pow(step[i],4)*a2[k]/24 
	         + pow(step[i],5)*a3[k]/120;
         xdot[i][k] = xdottemp[i][k] + pow(step[i],3)*a2[k]/6 
                    + pow(step[i],4)*a3[k]/24;
         f[i][k] = fi[k];
         fdot[i][k] = fidot[k];
         modfi += pow(fi[k],2);
         modfidot += pow(fidot[k],2);
         moda2 += pow(a2[k],2);
         moda3 += pow(a3[k],2);
      }
      dnmtr =  moda2 + sqrt(modfidot*moda3);
      assert(dnmtr!=0);
      temp = sqrt(0.02*(sqrt(modfi*moda2) + modfidot)/dnmtr);
      step[i] = min(1.2*step[i], temp);
      tlast[i] = *t;
      nsteps++;
  }
  while(*t<tend);
  return nsteps;
}

void ffdot(real x[NMAX][NDIM], real xdot[NMAX][NDIM], 
	   /**/real m[], int i, int n, real fi[], real fidot[]) {
      real rij[NDIM], rijdot[NDIM], r2, r3, r5, rrdot;
      loop(k, NDIM) { fi[k] = fidot[k] = 0;}
      loop(j, n) {
	if (j!=i) {
            r2 = rrdot = 0;
	    loop(k, NDIM) {
	      rij[k] = x[i][k] - x[j][k];
	      rijdot[k] = xdot[i][k] - xdot[j][k];
               r2 += pow(rij[k],2);
               rrdot += rij[k]*rijdot[k];
	    }
            r3 = pow(r2,1.5);
            r5 = r3*r2;
	    assert(r3!=0&&r5!=0);
            loop(k, NDIM) {
               fi[k] = fi[k] - m[j]*rij[k]/r3;
               fidot[k] = fidot[k] - m[j]*(rijdot[k]/r3 - 
					   /**/3*rrdot*rij[k]/r5);
	    }
	}
      }
}


void runtime_output(real x[NMAX][NDIM], real t, 
		    /**/long int n){
  real ss=0;
  loop (i,n) ss += x[i][0]*x[i][0] + x[i][1]*x[i][1] + 
    /**/x[i][2]*x[i][2];
  printf("%lf %lf\n",t,sqrt(ss)/n);
}

void final_output(real t, long int cpu, long int nsteps, 
		  /**/long int noutp){
  long int timenow;
  time(&timenow);
      printf("\n#%lf, %ld, %ld\n", 
/* 	     t, (clock()-cpu)/CLOCKS_PER_SEC, nsteps); */
	     t, (timenow-cpu), nsteps);
}  

