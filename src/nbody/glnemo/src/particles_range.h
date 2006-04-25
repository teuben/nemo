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
//                                                                             
// ParticlesRange class definition                                             
//                                                                             
// Parse selected string                                                       
// ============================================================================
#ifndef PARITCLES_RANGE_H
#define PARITCLES_RANGE_H
#include <qnamespace.h>
#include <iostream>
#include <vector>

#include "virtual_particles_select.h"

using namespace std;

class ParticlesRange;

class ParticlesRange: public VirtualParticlesSelect {
 public:

  ParticlesRange();

  ~ParticlesRange();

  // method
  
  int getIndex(int);
  int defaultIndexTab();
 private:
  // method
  int  parseSelectedString(char * select_sting, const int nbody,
                           ParticlesSelectVector * psv);
};
#endif
// ============================================================================
