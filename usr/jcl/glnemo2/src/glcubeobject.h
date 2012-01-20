// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
/**
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
*/
#ifndef GLCUBEOBJECT_H
#define GLCUBEOBJECT_H
#include <QGLWidget>

#include "globject.h"
namespace glnemo {
  
class GLCubeObject : public GLObject {
public:
  GLCubeObject(float,const QColor &c=Qt::yellow, bool activated=TRUE);
    ~GLCubeObject();
    void setSquareSize(float);
    void rebuild() { 
         buildDisplayList(); 
    }
private:
    float square_size;  // size of square
    // method
    void  buildDisplayList(); 
};
}
#endif // GLCUBEOBJECT_H
