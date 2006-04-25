// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2005                                  
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
#include <iostream>
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
  particles_alpha=255;
}
// ============================================================================
// Destructor                                                                  
GLObject::~GLObject()
{
} 
// ============================================================================
// GLObject::display()                                                         
// Display object, via display list, if activated                              
void GLObject::display()
{
  if (is_activated) {
    glColor4ub(mycolor.red(), mycolor.green(), mycolor.blue(),particles_alpha);
    glCallList(dplist_index);
  }
}
// ============================================================================
// GLObject::updateAlphaSlot()                                                 
// update particle alpha color                                                 
void GLObject::updateAlphaSlot(int _alpha)
{ 
  particles_alpha = _alpha;
}
// ============================================================================
// GLObject::toggleActivate()                                                  
// toggle activate                                                             
void GLObject::toggleActivate()
{
  is_activated = ! is_activated;
}
// ============================================================================
// GLObject::setColor()                                                        
// Set object color                                                            
void GLObject::setColor(const QColor &c)
{
  mycolor = c;
}
// ============================================================================
// GLObject::buildDisplayList()                                                
// Build display list                                                          
void GLObject::buildDisplayList()
{
}
// ============================================================================
