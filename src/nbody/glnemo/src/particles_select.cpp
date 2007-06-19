// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2007                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include "particles_select.h"

#include "particles_range.h" // 
#include "particles_list.h"

// ============================================================================
// constructor                                                                 
ParticlesSelect::ParticlesSelect()
{
  // virtual object
  vps = NULL;
}
// ============================================================================
// destructor                                                                  
ParticlesSelect::~ParticlesSelect()
{
  if (vps) delete vps;
}
// ============================================================================
// copy constructor for vectors                                                
ParticlesSelect::ParticlesSelect(const ParticlesSelect& m)
{
  m.vps->nb_select--; // decrease #object selected because      
                      // next new will increase it for free     
                      // and will change the color of the object
  switch (m.vps->v_type) {
  case 1:     
    vps = new ParticlesRange(); break;
  case 2: 
    vps = new ParticlesList(); break;
  default: 
    std::cerr << "Unexpected type for [ParticlesSelect], aborting....\n";
    std::exit(1);
  }
  *vps = *m.vps; // copy object
}
// ============================================================================
// copy constructor for vectors                                                
const ParticlesSelect& ParticlesSelect::operator=(const ParticlesSelect& m)
{
  switch (m.vps->v_type) {
  case 1: 
    vps = new ParticlesRange(); break;
  case 2: 
    //vps = new B(); break;
  default: 
    std::cerr << "Unexpected type for [ParticlesSelect], aborting....\n";
    std::exit(1);
  }
  *vps = *m.vps; // copy object
  return *this;
}
// ============================================================================
