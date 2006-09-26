// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2006                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#ifndef PARTICLES_SELECT_H
#define PARTICLES_SELECT_H
#include <vector>

//#include "virtual_particles_select.h"

class VirtualParticlesSelect;
class ParticlesSelect;

typedef std::vector <ParticlesSelect> ParticlesSelectVector;

class ParticlesSelect{
public:
    ParticlesSelect();
    // copy constructors for vectors
    ParticlesSelect(const ParticlesSelect&);
    const ParticlesSelect& operator=(const ParticlesSelect&);
    
    ~ParticlesSelect();
    VirtualParticlesSelect * vps; // virtual object of type:
                                  // ParticlesRange         
                                  // ParticlesList          
                                  
};

#endif
// ============================================================================
