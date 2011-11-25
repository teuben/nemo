// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2011                                  
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
GLTextObject::GLTextObject(const QString text,const fntRenderer &f,
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
GLTextObject::GLTextObject(const fntRenderer &f,const QColor &c, 
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
void GLTextObject::setFont(fntRenderer &f)
{
  font = f;
}
// ============================================================================
// GLTextObject::getLabelWidth()                                               
// return label width in pixels                                                
int GLTextObject::getLabelWidth()
{
  float l,r,b,t;
  font.getFont()->getBBox(label.toStdString().c_str(),font.getPointSize(),0,&l,&r,&b,&t);
  return (r-l);
}
// ============================================================================
// GLTextObject::getTextWidth()                                                
// return text width in pixels                                                 
int GLTextObject::getTextWidth()
{
  float l,r,b,t;
  font.getFont()->getBBox(text.toStdString().c_str(),font.getPointSize(),0,&l,&r,&b,&t);
  return (r-l);
}
// ============================================================================
// GLTextObject::getHeight()                                                   
// return font height in pixels                                                
int GLTextObject::getHeight()
{
  float l,r,b,t;
  //font.getFont()->getBBox("this is a test",font.getPointSize(),0,&l,&r,&b,&t);
  font.getFont()->getBBox(text.toStdString().c_str(),font.getPointSize(),0,&l,&r,&b,&t);
  return (t-b);
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
void GLTextObject::display(const int width, const int height)
{
  if (width) {;} // remove compiler warning
  if (is_activated) {   
    glColor4ub(mycolor.red(), mycolor.green(), mycolor.blue(),255);
    //std::cerr << "text and label" << label.toStdString()<< "/" << text.toStdString()<<"\n";
    // label    
    font.start2f(x,height-y);
    font.puts(label.toStdString().c_str());
    // text
    font.start2f(x_text,height-y);
    font.puts(text.toStdString().c_str());
  }
}
} // namespace
// ============================================================================

