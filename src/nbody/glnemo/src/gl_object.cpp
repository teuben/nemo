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
// GLObject class implementation                                               
//                                                                             
// Manage OpenGL  Object                                                       
// ============================================================================

#include "gl_object.h"
#include <qgl.h>

#define LOCAL_DEBUG 0
#include "print_debug.h"

// ============================================================================
// Constructor
GLObject::GLObject()
{
  is_activated=TRUE;
  this->setColor(yellow);
}
// ============================================================================
// Destructor
GLObject::~GLObject()
{
} 
// ============================================================================
// 
void GLObject::display()
{
  if (is_activated) {
    //glColor4ub(255,15,255,255);
    glColor4ub(mycolor.red(), mycolor.green(), mycolor.blue(),255);
    //qglColor(mycolor);
    glCallList(dplist_index);
  }
}
// ============================================================================
// 
void GLObject::toggleActivate()
{
  is_activated = ! is_activated;
}
// ============================================================================
// 
void GLObject::setColor(const QColor &c)
{
  mycolor = c;
}
// ============================================================================
// 
void GLObject::buildDisplayList()
{
}
