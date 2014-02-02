// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2014                                  
// e-mail:   Jean-Charles.Lambert@lam.fr                                      
// address:  Centre de donneeS Astrophysique de Marseille (CeSAM)              
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include "orbits.h"
#include <assert.h>
#include <QGLWidget>
namespace glnemo {
// ============================================================================
// Constructor                                                                 
Orbits::Orbits(const ParticlesData * p_data, const int index):GLObject()
{
  assert(p_data);
  assert(index<*p_data->nbody);
  
  time  = *(p_data->timu);
  pos[0] = p_data->pos[index*3];
  pos[1] = p_data->pos[index*3+1];
  pos[2] = p_data->pos[index*3+2];

}
// ============================================================================
// Copy Constructor                                                            
Orbits::Orbits(const Orbits & m):GLObject()
{
  time  = m.time;

  pos[0] = m.pos[0];
  pos[1] = m.pos[1];
  pos[2] = m.pos[2];

}
// ============================================================================
// Constructor                                                                 
Orbits::Orbits()
{

}
// ============================================================================
// Copy Constructor                                                            
const Orbits & Orbits::operator=(const Orbits & m)
{
  time  = m.time;

  pos[0] = m.pos[0];
  pos[1] = m.pos[1];
  pos[2] = m.pos[2];
  return *this;
}
// ============================================================================
// Destructor                                                                  
Orbits::~Orbits()
{

}
}
