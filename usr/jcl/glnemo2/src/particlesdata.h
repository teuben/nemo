// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2011                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
/**
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */
#ifndef GLNEMOPARTICLESDATA_H
#define GLNEMOPARTICLESDATA_H

#include <vector>
namespace glnemo {

class PhysicalData {
  public:
  enum ALLOC {New,Malloc};
  enum PHYS {neib,rho,temperature,pressure,temperaturesd};
  PhysicalData(const PHYS, const int _nbody=0, const ALLOC model=New);
  ~PhysicalData();
  const PhysicalData& operator=(const PhysicalData& m);
  void setNbody(const int n) {
    nbody = n;
  }
  int computeMinMax();
  double getMin() const { return min;}
  double getMax() const { return max;}
  bool isValid() const { return valid;}
  int getType()  const { return type;}
  void setType(const PHYS i)  { // set the index of the physical value selected
    if (i!=-1) {
      type = i; 
    }
  }
  float * data;
  int data_histo[100]; // store #particles per percentage
  private:
  int nbody;
  double min,max;
  bool valid;
  ALLOC cmodel;
  PHYS type;
};
class ParticlesData{
public:
    enum ALLOC {New,Malloc};

    ParticlesData(const ALLOC model=New );

    ~ParticlesData();
    const ParticlesData& operator=(const ParticlesData& m);


    int   * nbody, i_max[3], * nemobits;;
    float * pos, 
    * vel, * vel_norm,
    * timu, coo_max[3], coo_min[3];//* rneib, * rho, * temp, * pressure;
    int tree_size_max;
    std::vector <int> id;
    PhysicalData * rneib, * rho, * temp, * pressure;
    static char * mallocate(char *, int, bool force=false);
    void computeVelNorm();
    void computeMaxSize();
    float  getMaxSize() { computeMaxSize(); 
                          return max_size; }
    float getMaxVelNorm() const { return max_vel_norm; }
    void setIpvs(const int i=1)  { // set the index of the physical value selected
      if (i!=-1) {
        ipvs = i; 
      }
    }
    int getIpvs() const {
      return ipvs;
    }
    PhysicalData * getPhysData(int=-1) const;
    private:
     
      float max_vel_norm;
      float max_size;
      ALLOC cmodel;
      int allocVar();
      int ipvs; // Index of the physical value selected
};

}

#endif
