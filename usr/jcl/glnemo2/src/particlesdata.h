// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2009                                  
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

namespace glnemo {

class ParticlesData{
public:
    enum ALLOC {New,Malloc};

    ParticlesData(const ALLOC model=New );

    ~ParticlesData();
    const ParticlesData& operator=(const ParticlesData& m);


    int   * nbody, i_max[3], * nemobits;;
    float * pos, 
    * vel, * vel_norm,
    * timu, coo_max[3],* rneib, * rho, * temp;
    int tree_size_max;
    char * mallocate(char *, int, bool force=false);
    void computeVelNorm();
    float getMaxVelNorm() const { return max_vel_norm; };
    float getMaxRho() const { return max_rho; };
    float getMinRho() const { return min_rho; };
    float getMaxTemp() const { return max_temp; };
    float getMinTemp() const { return min_temp; };
    bool isRhoValid() const { return valid_rho;}
    bool isTempValid() const { return valid_temp;}
    int computeMinMaxRho();
    int computeMinMaxTemp();
    int density_histo[100];
    int temp_histo[100];
    private:
      float max_vel_norm;
      float max_rho, min_rho;
      float max_temp, min_temp;
      bool valid_rho, valid_temp;
      const ALLOC cmodel;
      int allocVar();
};

}

#endif
