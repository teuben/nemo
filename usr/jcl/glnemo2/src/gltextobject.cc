// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2010                                  
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
//                                                                             
// GLTextObject class implementation                                           
//                                                                             
// Manage OpenGL Text Object                                                   
// ============================================================================

#include "glwindow.h"
#include "gltextobject.h"

namespace glnemo {
// ============================================================================
// constructor                                                                 
GLTextObject::GLTextObject(const QString text,const QFont &f,
                           const QColor &c, 
                           bool activated):GLObject()
{
  if (text=="") {;};
  setColor(c);
  font = f;
  is_activated = activated;
  x = y = x_text = 0; 
}
// ============================================================================
// constructor                                                                 
GLTextObject::GLTextObject(const QFont &f,const QColor &c, 
                           bool activated):GLObject()
{
  setColor(c);
  is_activated = activated;
  font = f;  
  x = y = x_text = 0;
}
// ============================================================================
// constructor                                                                 
GLTextObject::GLTextObject(bool activated):GLObject()
{
  //dplist_index = glGenLists( 1 ); // create a new display List
  is_activated = activated;
  x = y = x_text = 0;
}
// ============================================================================
// destructor                                                                  
GLTextObject::~GLTextObject()
{
}    
// ============================================================================
// GLTextObject::setText()                                                     
// set label and text                                                          
void GLTextObject::setText(const QString &p_label,const QString& p_text)
{
  label = p_label;
  text  = p_text;
}
// ============================================================================
// GLTextObject::setFont()                                                     
// set Font                                                                    
void GLTextObject::setFont(QFont &f)
{
  font = f;
}
// ============================================================================
// GLTextObject::getLabelWidth()                                               
// return label width in pixels                                                
int GLTextObject::getLabelWidth()
{
  QFontMetrics fm(font);
  return (fm.width(label,label.length()));
}
// ============================================================================
// GLTextObject::getTextWidth()                                                
// return text width in pixels                                                 
int GLTextObject::getTextWidth()
{
  QFontMetrics fm(font);
  return (fm.width(text,text.length()));
}
// ============================================================================
// GLTextObject::getHeight()                                                   
// return font height in pixels                                                
int GLTextObject::getHeight()
{
  QFontMetrics fm(font);
  return (fm.height());
}
// ============================================================================
// GLTextObject::setPos()                                                      
// specify new positions x and y (in pixels) for the given text                
void GLTextObject::setPos(const int new_x, const int new_y, 
                          const int new_x_text)
{
  x = new_x;
  y = new_y;
  x_text = new_x_text;
}
// ============================================================================
// GLTextObject::display()                                                     
// display text object if activated                                            
void GLTextObject::display(GLWindow * glbox)
{
  if (is_activated) {   
    glColor4ub(mycolor.red(), mycolor.green(), mycolor.blue(),255);
    // It's not possible to use renderText method inside a display List,
    // because 'renderText' already use a display list. It drove me nuts
    // during a while :(.
    glbox->renderText(x,y,label,font);
    //std::cerr << x_text << " X "<<y<< " : " << text.toStdString() << "\n";
    glbox->renderText(x_text,y,text,font);
  }
}
} // namespace
// ============================================================================

