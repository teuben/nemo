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
#include "glcube.h"

// ============================================================================
// Constructor                                                                 
GLCube::GLCube(float _square_size, const QColor &c, bool activated ):GLObject()
{
  square_size = _square_size;
  dplist_index = glGenLists( 1 );
  buildDisplayList();
  setColor(c);
  is_activated=activated;
}

// ============================================================================
// Destructor                                                                  
// Delete display list                                                         
GLCube::~GLCube()
{
  glDeleteLists( dplist_index, 1 );
}

// ============================================================================
// GLCube::buildDisplayList()                                            
// Build Display List                                                          
void GLCube::buildDisplayList()
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
// GLCube::setSquareSize()                                                     
void GLCube::setSquareSize(float _ss)
{
  square_size = _ss;
}
// ============================================================================
