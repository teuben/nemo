// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
/**
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */
#ifndef GLNEMOPARTICLESSELECTRANGE_H
#define GLNEMOPARTICLESSELECTRANGE_H

#include "particlesobject.h"

namespace glnemo {

// ----------------------------------------------------------------------------
// class ParticlesSelectRange parse a string (nemoinpi format) and create a vector
// of objects
class ParticlesSelectRange{
  public:
    ParticlesSelectRange(const int _nbody,              // nbody in the snapshot  
                         const std::string _select,     // particles to select
                         const bool _keep_all,          // keep all particles ?   
                         ParticlesObjectVector * pov);  // vector to store objects

    ~ParticlesSelectRange();

  private:
    int nbody;
    std::string select;
    bool keep_all;
    ParticlesObjectVector * pov;
    
    void parseSelect();
    std::string parseString(std::string&);

    void storeObject(const std::string);
};

}

#endif
