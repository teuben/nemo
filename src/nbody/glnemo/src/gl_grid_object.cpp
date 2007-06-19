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
// GLGridObject class implementation                                           
//                                                                             
// Manage OpenGL Grid Object                                                   
// ============================================================================

#include "gl_grid_object.h"
#include <iostream>

#define LOCAL_DEBUG 0
#include "print_debug.h"

int   GLGridObject::nsquare=27;        // default #nsquare   
float GLGridObject::square_size=2.0;   // default square size

using namespace std;
// ============================================================================
// Constructor                                                                 
GLGridObject::GLGridObject(int axe_parm, const QColor &c, bool activated ):GLObject()
{
  if (axe_parm < 0 || axe_parm > 2 ) 
    axe_parm=0;
  dplist_index = glGenLists( 1 );
  axe = axe_parm;
  buildDisplayList();
  setColor(c);
  is_activated=activated;
  
}

// ============================================================================
// Destructor                                                                  
// Delete display list                                                         
GLGridObject::~GLGridObject()
{
  glDeleteLists( dplist_index, 1 );
}
// ============================================================================
// GLGridObject::buildDisplayList()                                            
// Build Display List                                                          
void GLGridObject::buildDisplayList()
{
  float x=0.,y=0.,z=0.,x1=0.,y1=0.,z1=0.;
  GLfloat 
    inf=((GLfloat) (-nsquare/2)*square_size),
    sup=((GLfloat) ( nsquare/2)*square_size);

  // display list
  glNewList( dplist_index, GL_COMPILE );
  glBegin(GL_LINES);
  // draw all the lines
  for (int i=-nsquare/2; i<=nsquare/2; i++) {
    switch (axe) {
    case 0:                  // X/Y plane
      x = x1 = square_size * (GLfloat) i;
      y = inf; y1 = sup;
      z=0.0; z1=0.0;
      break;
    case 1:                  // Y/Z plane
      y = y1 = square_size * (GLfloat) i;
      z = inf; z1 = sup;
      x=0.0; x1=0.0;
      break;
    case 2:                 // X/Z plane
      z = z1 = square_size * (GLfloat) i;
      x = inf; x1 = sup;
      y=0.0; y1=0.0;
      break;
    }
    // one line
    glVertex3f(x , y  ,z );
    glVertex3f(x1, y1 ,z1);

  }
  // perpendicular line
  for (int i=-nsquare/2; i<=nsquare/2; i++) {
    switch (axe) {
    case 0:                  // X/Y plane
      y = y1 = square_size * (GLfloat) i;
      x = inf; x1 = sup;
      z=0.0; z1=0.0;
      break;
    case 1:                  // Y/Z plane
      z = z1 = square_size * (GLfloat) i;
      y = inf; y1 = sup;
      x=0.0; x1=0.0;
      break;
    case 2:                 // X/Z plane
      x = x1 = square_size * (GLfloat) i;
      z = inf; z1 = sup;
      y=0.0; y1=0.0;
      break;
    }
    // one line
    glVertex3f(x , y  ,z );
    glVertex3f(x1, y1 ,z1);
  }
  glEnd();
  glEndList();
}
// ============================================================================
