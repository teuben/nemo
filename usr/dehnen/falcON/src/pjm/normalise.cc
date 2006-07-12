
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

  double mfull=atof(argv[1]),
    inner=atof(argv[2]),
    outer=atof(argv[3]),
    rscale=atof(argv[4]),
    rtrunc=atof(argv[5]);

  RHO  = new Halo(inner,outer,rscale,rtrunc);
  RHO2 = new Halo(inner,outer,rscale,60); 
  const double rmin=0.001, rmax=10000;
  int n1 = int(200*log10(rmax/rmin));
  int  n  = 1+n1;
  const double dlr= log(rmax/rmin)/double(n1);
  double M=rmin*rmin*rmin*pow(rmin/rscale,-inner)/(3-inner)
	+ rmin*rmin*rmin*(inner-outer)*pow(rmin/rscale,1.-inner)/(4-inner);
  double M2 = M;
  double lr = log(rmin);
  for(int i=1; i!=n; ++i) {
    M     = rk4(M,lr,dlr,dM); //Integrate mass outwards
    M2     = rk4(M2,lr,dlr,dM2);
    lr+=dlr;
  }

  cout << mfull*M/M2 << "\n";
  return 0;
}
