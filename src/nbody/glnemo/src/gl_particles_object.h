// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004                                       
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

//#include <qgl.h>
#include "gl_object.h"
#include <vector>
#include "particles_range.h"


//using namespace std;

class GLParticlesObject;

//typedef vector <GLParticlesObject> GLParticlesObjectVector;

class GLParticlesObject : public GLObject {
  Q_OBJECT
 public:
  GLParticlesObject() { is_activated=FALSE; }
  GLParticlesObject(const int * _nbody, const float * _pos, 
                    const ParticlesRange * _prv);

  int updateObject(const int * _nbody, const float * _pos, 
                   const ParticlesRange * _prv);

  ~GLParticlesObject();
  //int getNpart() { return npart; };
  // method
  void displayPolygons(const double * mModel,GLuint texture,float u_max, float v_max);                              
  
 private:

  const int * nbody;
  const float * pos;
  const ParticlesRange * prv;
  
  // method
  void  buildDisplayList(const int * nbody, const float * pos, 
                         const ParticlesRange * prv);



};
#endif
// ============================================================================
