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
#ifndef CROTGAL_H
#define CROTGAL_H
#include <vector>
#include <string>
#include "uns.h"
#include "cfalcon.h"

namespace uns_proj {
  class CRotgal;
  
  //------------------------------------------------------------------------------
  // CPartVec
  // class to store rho, index and id of particles
  class CPartVec {
  public:
    CPartVec(CRotgal * _rotgal, int _index) {
      index = _index;
      rotgal = _rotgal;
    }
    static bool sortRho(const CPartVec& a, const CPartVec& b);
    static bool sortId (const CPartVec& a, const CPartVec& b);
    int index;
    CRotgal * rotgal;  
    float computeR2();
  private:
    
  };
  //------------------------------------------------------------------------------
  // CPartRT
  // class to store radius and theta of selected particles
  class CPartRT {
  public:
    CPartRT(float _diff, float _theta) {
      diff_radius = _diff;
      theta       = _theta;
    } 
    static bool sortRadius(const CPartRT& a, const CPartRT& b);
    static bool sortTheta(const CPartRT& a, const CPartRT& b);
    float diff_radius,theta;
  };
    
  //------------------------------------------------------------------------------
  // CRotgal
  // to manage all the process
  class CRotgal
  {    
  public:
    CRotgal(uns::CunsIn * uns);
    ~CRotgal();
    bool loadData();
    void process();
    void saveSelectPart(std::string,std::vector <CPartVec> *);
    void saveSelectPart(std::vector <CPartVec> *);
    void selectPart();
    void computeRotation();
    void sortRho() {
      // Descending sort particles according to their densities
      std::sort(pvec.begin(),pvec.end(),CPartVec::sortRho);
    }
    void computeRadiusTheta(CPartVec *,CPartVec *);
    jclut::CDensity * getDensity() { return density;}
    int nbody;
    std::vector<float> pos,vel,mass,hsml,rho;
    std::vector<int> id;
    std::vector <CPartVec> pvecselect;
    float time;
    
  private:
    uns::CunsIn * unsin;
    jclut::CDensity * density;
    void clearVectors() {
      pos.clear(); vel.clear(); mass.clear(); hsml.clear(); rho.clear(); id.clear();
    }
    std::vector <CPartVec> pvec;
    std::vector <CPartRT > prtvec;
    
  };
}
#endif // CROTGAL_H
