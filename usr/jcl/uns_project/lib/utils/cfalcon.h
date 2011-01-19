// ============================================================================
// Copyright Jean-Charles LAMBERT - 2010                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

/*
  @author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
*/
#ifndef CFALCON_H
#define CFALCON_H
#include <forces.h>                                // falcON      
#include <body.h>                                  // the bodies    
#include <public/neighbours.h>                     // finding neighbours       
namespace jclut {
  
class cfalcon
{
public:
    cfalcon();
    static bool addGravity(const int nbody,
                           const float * pos, const float * mass, 
                           float * acc, float * phi,
                           const float eps, 
                           const float G=1.0,
                           const float theta=0.6,
                           const int kernel_type=1,
                           const int ncrit=6);
    
    static bool addGravity2(const int nbody,
                             const float * pos, const float * mass, 
                             const int nbody_tp,
                             const float * pos_tp,
                             float * acc, float * phi,
                             const bool selfp,
                             const float eps,
                             const float G=1.0,
                             const float theta=0.6,
                             const int kernel_type=1,
                             const int ncrit=6);
    
};
using namespace falcON;
using falcON::real;
// -------------------------------------------------------------
// Class CDensity
// estimate density using falcON engine
class CDensity {
  
public:
  CDensity(const int nbody, float * pos, float * mass);
  ~CDensity();
  void setData(const int nbody, float * pos, float * mass);
  void compute(const int method, const int K,const int N, const int ncrit);
  static real F; ///< normalisation factor for kernel
  static int  N; ///< order of Ferrers kernel
  real * getRho()  { return rho; }
  real * getHsml() { return hsml; }
private:
  real * rho, * hsml;
  int nbody;
  //real RHO,AUX;
  snapshot *  SHOT ;
  void prepare(int n) {
    N = n;
    F = 0.75/Pi;
    for(n=1; n<=N; ++n)
      F *= double(n+n+3)/double(n+n);
  } 
  static void SetDensity(const bodies*B, const OctTree::Leaf*L,
		  const Neighbour*NB, int K)
  {
    real iHq = one/NB[K-1].Q;
    real rho = zero;
    for(int k=0; k!=K-1; ++k) {
      //std::cerr << "rneib["<<k<<"]="<<NB[k].Q<<"\n";
      rho += scalar(NB[k].L) * std::pow(one-iHq*NB[k].Q,N);
    }
    rho *= F * std::pow(sqrt(iHq),3);
    B->rho(mybody(L)) = rho;
    B->aux(mybody(L)) = sqrt(NB[K-1].Q);
    //std::cerr << "neib="<<B->aux(mybody(L))<<"\n";

  }
#define   FRTHRD_PI  4.18879020478639098462
  static void SetDensity2(const bodies*B, const OctTree::Leaf*L,
		  const Neighbour*NB, int K)
  {
    float rn=NB[K-1].Q;
    float sqrn = sqrt(rn);
    real rho=(K-1)/(rn*sqrn*FRTHRD_PI);
    B->rho(mybody(L)) = rho;
    B->aux(mybody(L)) = sqrn;
    //std::cerr << "neib="<<B->aux(mybody(L))<<"\n";
    //fprintf(stdout,"r=%f nb=%d den=%f\n", sqrt(rn),K,rho);
  }
};
}
#endif // CFALCON_H
