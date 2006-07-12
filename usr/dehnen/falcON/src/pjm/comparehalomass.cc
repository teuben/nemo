//
// comparehalomass.cc
//------------------------------------------------------------------------
//
// Program which gives mass of halo 
// rho \propto ((r/rscale)^-inner)*((1+(r/rscale))^-(outer-inner))*sech(r/rtrunc) 
//
// with Minc at r_0 same as that of the second halo defined by the user 
// if rtrunc is given as 0, it is assumed to be infinite

// Inputs are:
// argv[1] : inner density exponent of desired halo
// argv[2] : outer density exponent of desired halo
// argv[3] : halo scale radius of desired halo
// argv[4] : halo trunction radius of desired halo (if 0, no truncation)
// argv[5] : mass of comparison halo
// argv[6] : inner density exponent of comparison halo
// argv[7] : outer density exponent of comparison halo
// argv[8] : halo scale radius of comparison halo
// argv[9] : halo trunction radius of comparison halo
// argv[10]: r_0, the radius within which the masses of the 
//           desired and comparison halos are equal

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>

using namespace std;

// Halo Model -------------------------------------------------------------
//
class Halo {
const double b,g,rs,irs,rt;
  double rsont;
public:
  Halo(const double _g,
       const double _b,
       const double _rs,
       const double _rt):
    b   (_b),                 // Outer density exponent
    g   (_g),                 // Inner density exponent
    rs  (_rs),                // Scale radius
    irs (1./rs),
    rt  (_rt)                 // Truncation radius (if required)
  {
    rsont=(rt)? rs/rt : 0;
  }

double operator() (const double r) const
    {
      register double x=r*irs,
	x1=1.+ x;
      if(rt){return pow(x,-g)*pow(x1,(g-b))/(exp(x*rsont)+exp(-x*rsont));}
      else  return pow(x,-g)*pow(x1,(g-b)); // put density here
    }

  ~Halo();

};
//--------------------------------------------------------------------------------
  
Halo *RHO,*RHO2;

  inline double dM(double lr, double M)            // dM/dln r                  
{
  register double r = exp(lr);                   // ln r as independent var   
  return r*r*r*(*RHO)(r);
}

  inline double dM2(double lr, double M)            // dM/dln r                  
{
  register double r = exp(lr);                   // ln r as independent var   
  return r*r*r*(*RHO2)(r);
}
 // Integrator --------------------------------------------------------------------
 //
template<typename scalar_type, typename vector_type>
   inline
   vector_type rk4(const vector_type y,
		  const vector_type dy0,
		  const scalar_type x,
		  const scalar_type h,
		  vector_type(*derivs)(const scalar_type, const vector_type))
  {
    const    scalar_type hh=0.5*h, h6=h/6., xh=x+hh;
    register vector_type yt,dym,dyt;
    dyt = derivs(xh,y+hh*dy0);
    dym = derivs(xh,y+hh*dyt);
    yt  = y+h*dym;
    dym+= dyt;
    dyt = derivs(x+h,yt);
    return y+h6*(dy0+dyt+dym+dym);
  }
  //----------------------------------------------------------------------------
  template<typename scalar_type, typename vector_type>
  inline
  vector_type rk4(const vector_type y,
		  const scalar_type x,
		  const scalar_type h,
		  vector_type(*derivs)(const scalar_type, const vector_type))
  {
    return rk4(y,derivs(x,y),x,h,derivs);
  }

  //----------------------------------------------------------------------------

// Program which gives mass of halo with truncation radius rtrunc and same 
// central density distribution as a halo with rtrunc=60


int main(int argc, char* argv[])
{
  if(argc!=11){
    cerr << "wrong inputs\n";
    exit(1);
  }

  double inner  = atof(argv[1]),
         outer  = atof(argv[2]),
         rscale = atof(argv[3]),
         rtrunc = atof(argv[4]),
         Mhalo2 = atof(argv[5]),
         inner2 = atof(argv[6]),
         outer2 = atof(argv[7]),
         rscale2= atof(argv[8]),
         rtrunc2= atof(argv[9]),
         rcomp  = atof(argv[10]);
  double Mrcomp,Mrcomp2;
  bool   found=false;
  RHO  = new Halo(inner,outer,rscale,rtrunc);
  RHO2 = new Halo(inner2,outer2,rscale2,rtrunc2); 
  const double rmin=0.00001*rcomp, rmax=100000*rcomp;
  if(rcomp=0.) rcomp=rmax;
  int n1 = int(200*log10(rmax/rmin));
  int  n  = 1+n1;
  const double dlr= log(rmax/rmin)/double(n1);
  double M=rmin*rmin*rmin*pow(rmin/rscale,-inner)/(3-inner)
	+ rmin*rmin*rmin*(inner-outer)*pow(rmin/rscale,1.-inner)/(4-inner);
  double M2 = rmin*rmin*rmin*pow(rmin/rscale2,-inner2)/(3-inner2)
	+ rmin*rmin*rmin*(inner2-outer2)*pow(rmin/rscale2,1.-inner2)/(4-inner2);
  double lr = log(rmin);
  for(int i=1; i!=n; ++i) {
    M     = rk4(M,lr,dlr,dM); //Integrate mass outwards
    M2     = rk4(M2,lr,dlr,dM2);
    if(exp(lr)>0.999999*rcomp && !found) {
      Mrcomp  = M;
      Mrcomp2 = M2;
      found=true;
    }
    lr+=dlr;
  }

  cout <<  Mhalo2*M*Mrcomp2/(M2*Mrcomp) << "\n";
  return 0;
}
