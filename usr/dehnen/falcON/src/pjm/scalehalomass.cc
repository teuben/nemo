//
// scalehalomass.cc
//------------------------------------------------------------------------
//
// Program which gives mass of halo 
// rho \propto ((r/rscale)^-inner)*((1+(r/rscale))^-(outer-inner))*sech(r/rtrunc) 
//
// with Minc at r_0 as defined by user 
// if rtrunc is given as 0, it is assumed to be infinite
// Inputs are:
// argv[1]: Mass inside r = r_0
// argv[2]: the aforementioned r_0
// argv[3]: inner density exponent
// argv[4]: outer density exponent
// argv[5]: halo scale radius
// argv[6]: halo trunction radius (if 0, no truncation)

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
  
Halo *RHO;

  inline double dM(double lr, double M)            // dM/dln r                  
{
  register double r = exp(lr);                   // ln r as independent var   
  return r*r*r*(*RHO)(r);
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


int main(int argc, char* argv[])
{
  if(argc!=7){
    cerr << "wrong inputs\n";
    exit(1);
  }
  double mbyr_0=atof(argv[1]),
    r_0=atof(argv[2]),
    inner=atof(argv[3]),
    outer=atof(argv[4]),
    rscale=atof(argv[5]),
    rtrunc=atof(argv[6]);

  RHO  = new Halo(inner,outer,rscale,rtrunc);
  const double rmin=0.000001*r_0, rmax=1000000*r_0;
  double M2;
  bool found=false;
  int n1 = int(200*log10(rmax/rmin));
  int  n  = 1+n1;
  const double dlr= log(rmax/rmin)/double(n1);
  double M=rmin*rmin*rmin*pow(rmin/rscale,-inner)/(3-inner)
	+ rmin*rmin*rmin*(inner-outer)*pow(rmin/rscale,1.-inner)/(4-inner);
  double lr = log(rmin);
  for(int i=1; i!=n; ++i) {
    M     = rk4(M,lr,dlr,dM); //Integrate mass outwards
    if( exp(lr)>0.999999*r_0 && !found) {
      M2=M;
      found=true;
    }
    lr+=dlr;
  }

  cout << M*mbyr_0/M2 << "\n";
  return 0;
}
