// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2015                                  
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
/**
	@author Jean-Charles Lambert <jean-charles.lambert@lam.fr>
 */

#include "glcubeobject.h"

namespace glnemo {
// ============================================================================
// Constructor    
GLCubeObject::GLCubeObject(float  _square_size,const QColor &c, bool activated)
{
  square_size = _square_size;
  dplist_index = glGenLists( 1 );
  buildDisplayList();
  setColor(c);
  is_activated=activated;
}
// ============================================================================
// Destructor    
GLCubeObject::~GLCubeObject()
{
  glDeleteLists( dplist_index, 1 );
}
// ============================================================================
// GLCubeObject::buildDisplayList()                                            
// Build Display List                                                          
void GLCubeObject::buildDisplayList()
{

  float hs = square_size/2.; // half square
  // display list
  glNewList( dplist_index, GL_COMPILE );
  
  // bottom square
  glBegin(GL_LINE_STRIP);  
  glVertex3f(-hs , -hs  , -hs );
  glVertex3f( hs , -hs  , -hs );
  glVertex3f( hs ,  hs  , -hs );
  glVertex3f(-hs ,  hs  , -hs );
  glVertex3f(-hs , -hs  , -hs );  
  glEnd();
  
  // top square
  glBegin(GL_LINE_STRIP);  
  glVertex3f(-hs , -hs  , hs );
  glVertex3f( hs , -hs  , hs );
  glVertex3f( hs ,  hs  , hs );
  glVertex3f(-hs ,  hs  , hs );
  glVertex3f(-hs , -hs  , hs );  
  glEnd();
  
  // segment between two squares
  glBegin(GL_LINES);
  glVertex3f(-hs , -hs  , -hs );
  glVertex3f(-hs , -hs  ,  hs );
  
  glVertex3f( hs , -hs  , -hs );
  glVertex3f( hs , -hs  ,  hs );
  
  glVertex3f( hs ,  hs  , -hs );
  glVertex3f( hs ,  hs  ,  hs );
  
  glVertex3f(-hs ,  hs  , -hs );
  glVertex3f(-hs ,  hs  ,  hs );
  glEnd();
  glEndList();
}
// ============================================================================
// GLCubeObject::setSquareSize()                                                     
void GLCubeObject::setSquareSize(float _ss)
{
  square_size = _ss;
  buildDisplayList();
}
} // end of namespace
// ============================================================================
