// ============================================================================
// Copyright Jean-Charles LAMBERT - 2011                                       
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
#include "cfalcon.h"
#include <vector>

#ifndef CBAR_H
#define CBAR_H
using namespace jclut;
namespace uns_proj {
  class CBar;
  //------------------------------------------------------------------------------
  // CPartVec
  // class to store rho, index and id of particles
  class CVecRho {
  public:
    CVecRho(CBar * _bar, int _index) {
      index = _index;
      bar = _bar;
    }
    static bool sortRho(const CVecRho& a, const CVecRho& b);
    static bool sortId (const CVecRho& a, const CVecRho& b);
    int index;
    CBar * bar;      
  };
  
  class CBar {
  public:
    CBar(const int _nbody, float * _pos, float * _vel, float * mass, 
         float * _rho, float *_hsml);
    ~CBar();
    float * getRho() { return rho;}
    float computeAngle(const float dmin, const float dmax, const bool mvcod=false);
    float computeAngle(const bool mvcod=false);
    void rotateOnX(const float);
    void rotateOnY(const float);
    void saveAllRho(std::string out);
    void save(std::string out, const float timu,const bool mvcod);
  private:
    int nbody;
    float * pos, * vel, *mass, *rho, *hsml;
    CDensity * density;
    int data_histo[100]; // store #particles per percentage
    std::vector <CVecRho> vec_rho;
    
    void sortRho();
    void rotate(const float angle);
   
  }; // end of class
} // namespace

#endif // CBAR_H
