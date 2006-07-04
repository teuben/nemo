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
#ifndef GLCUBE_H
#define GLCUBE_H

#include "gl_object.h"
/**
@author Jean-Charles Lambert
*/
class GLCube: public GLObject {
public:
    GLCube(float,const QColor &c=yellow, bool activated=TRUE);

    ~GLCube();
    void setSquareSize(float);
    void rebuild() { buildDisplayList(); };
private:
  float square_size;  // size of square
  // method
  void  buildDisplayList( ); 
};

#endif
