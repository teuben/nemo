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
#include "particles_data.h"

class GLParticlesObject;

class GLParticlesObject : public GLObject {
  Q_OBJECT
 public:
  GLParticlesObject() { is_activated=FALSE; }
  GLParticlesObject(const ParticlesData *,
                    VirtualParticlesSelect *,
                    float vel_resize_factor=1.0);

  int updateObject(const ParticlesData *,
                   VirtualParticlesSelect *,
                   float vel_resize_factor=1.0);

  void displayVelVector();
    
  ~GLParticlesObject();
  // method
  void displayPolygons(const double * mModel,GLuint texture,float u_max, float v_max);
  void displaySprites(GLuint texture);
  void rebuildDisplayList(float _vel_resize_factor=-1.);
  void rebuildVelDisplayList(float );                              
  void setTextureSize(const float ts) {  texture_size=ts;};
  void setTextureAlphaColor(const int alpha) {  texture_alpha_color=alpha; };
  float coo_max[3];
  int i_max[3];
 private:

  // Data
  const ParticlesData * p_data;
  VirtualParticlesSelect * vps;
  // texture
  float texture_size;
  int texture_alpha_color;
  // velocity
  float vel_resize_factor;
  GLuint vel_dp_list;
  // method
  void buildDisplayList(const ParticlesData * , 
                         VirtualParticlesSelect *);
  void buildVelDisplayList(const ParticlesData * , 
                         VirtualParticlesSelect *);

  void computeCooMax();
  FrustumCulling frustum;
};
#endif
// ============================================================================
