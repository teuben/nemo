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
#include "globject.h"
#include <GL/glu.h>

#define LOCAL_DEBUG 0
//#include "print_debug.h"

namespace glnemo {

#define DOF 4000000
//float ortho_left,ortho_right,ortho_bottom,ortho_top;
// ============================================================================
// Constructor                                                                 
GLObject::GLObject()
{
  is_activated=TRUE;
  this->setColor(Qt::yellow);
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
void GLObject::display(int my_list)
{
  if (is_activated) {
    if (my_list == -1) my_list=dplist_index;
    glColor4ub(mycolor.red(), mycolor.green(), mycolor.blue(),particles_alpha);
    glCallList((GLuint) my_list);
  }
}
// ============================================================================
// GLObject::updateAlphaSlot()                                                 
// update particle alpha color                                                 
void GLObject::updateAlphaSlot(const int _alpha)
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
// GLObject::setProjection                                              
//
void GLObject::setProjection(const int x, const int y, const int width, const int height)
{
  glViewport( x, y, width, height);
  double ratio =  ((double )width) / ((double )height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if (1) {
    gluPerspective(45.,ratio,0.0005,(float) DOF);
  }
  else {
#if 1
    //computeOrthoFactor();
    float range=3;
    float fx,fy=fx=1.0;
    float ortho_right = range;
    float ortho_left  =-range;
    float ortho_top   = range;
    float ortho_bottom=-range;
    glOrtho(ortho_left   * fx  * 1.,
            ortho_right  * fx  * 1.,
            ortho_bottom * fy  * 1.,
            ortho_top    * fy  * 1.,
            -1000,1000);
            //(float) -DOF/2.,(float) -DOF/2.);
#endif
  }
}

// ============================================================================

// ============================================================================
// compute some factors for the orthographic projection
void GLObject::computeOrthoFactor()
{
  if (ratio<1.0) {
    fx = 1.0  ; fy = 1./ratio;
  }
  else {
    fx = ratio; fy = 1.0;
  }
}
}
