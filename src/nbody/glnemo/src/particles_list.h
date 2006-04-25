// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2005                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#ifndef PARTICLES_LIST_H
#define PARTICLES_LIST_H

#include <fstream>
#include <vector>

#include "virtual_particles_select.h"

class ParticlesList: public VirtualParticlesSelect {
public:
    ParticlesList();
    
    ~ParticlesList();
    int getIndex(int);
    int loadFile(const char *,ParticlesSelectVector * );  
    int defaultIndexTab();
 private:
  // variables
  static const std::string header;
  
  // method
  int  parseSelectedString(char * select_sting, const int nbody,
                           ParticlesSelectVector * psv);
};

#endif
// ============================================================================
