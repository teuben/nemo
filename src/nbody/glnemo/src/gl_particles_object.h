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
// GLParticlesObject class definition                                          
//                                                                             
// Manage OpenGL Particles Object                                              
// ============================================================================
#ifndef GL_PARTICLES_OBJECT_H
#define GL_PARTICLES_OBJECT_H

#include "gl_object.h"
#include <vector>
#include "virtual_particles_select.h"
#include "frustumculling.h"
class GLParticlesObject;

class GLParticlesObject : public GLObject {
  Q_OBJECT
 public:
  GLParticlesObject() { is_activated=FALSE; }
  GLParticlesObject(const int * _nbody, const float * _pos, 
                    VirtualParticlesSelect *);

  int updateObject(const int * _nbody, const float * _pos, 
                   VirtualParticlesSelect *);

  ~GLParticlesObject();
  // method
  void displayPolygons(const double * mModel,GLuint texture,float u_max, float v_max);
  void displaySprites(GLuint texture);
  void rebuildDisplayList();                              
  void setTextureSize(const float ts) {  texture_size=ts;};
  void setTextureAlphaColor(const int alpha) {  texture_alpha_color=alpha; };
  float coo_max[3];
  int i_max[3];
 private:

  const int * nbody;
  const float * pos;
  VirtualParticlesSelect * vps;
  float texture_size;
  int texture_alpha_color;
  // method
  void  buildDisplayList(const int * nbody, const float * pos, 
                         VirtualParticlesSelect *);

  void computeCooMax();
  FrustumCulling frustum;
};
#endif
// ============================================================================
