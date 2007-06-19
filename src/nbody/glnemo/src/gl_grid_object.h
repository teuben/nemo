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
// GLGridObject class definition                                               
//                                                                             
// Manage OpenGL Grid Object                                                   
// ============================================================================
#ifndef GL_GRID_OBJECT_H
#define GL_GRID_OBJECT_H

#include <qgl.h>
#include "gl_object.h"
#include "global_options.h"

class GLGridObject : public GLObject {
 public:
  GLGridObject(int axe_parm=0,const QColor &c=yellow, bool activated=TRUE);
  
  ~GLGridObject();

  // method
  void setNbSquare(int _nsquare) { nsquare = _nsquare;};
  void setSquareSize(float _square_size) { square_size = _square_size;};
  void rebuild() { buildDisplayList(); };
 private:
  static int nsquare;        // #square       
  static float square_size;  // size of square
  int axe;
  // method
  void  buildDisplayList( );
};
#endif
// ============================================================================
