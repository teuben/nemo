// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2011                                  
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
#ifndef GLAXESOBJECT_H
#define GLAXESOBJECT_H
#include "globject.h"
#include <GL/glu.h>

namespace glnemo {

class GLAxesObject: public GLObject {
public:
  GLAxesObject();
  ~GLAxesObject();
  void display(const double * mScreen,const double * mScene, const int width, const int height);
private:
  void buildDisplayList();
  void buildDisplayList2();
  void buildArrow(const float length, const float radius, const int nbSubdivisions);
  GLUquadric *quadric;
};

} // namespace
#endif // GLAXESOBJECT_H
