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
// GLTextObject class implementation                                           
//                                                                             
// Manage OpenGL Text Object                                                   
// ============================================================================

#include "gl_text_object.h"
#include "glbox.h"


// ============================================================================
// constructor
GLTextObject::GLTextObject(const QString text,const QFont &f,
                           const QColor &c, 
                           bool activated):GLObject()
{
  if (text);  // do nothing (remove compiler warning)
  
  //dplist_index = glGenLists( 1 ); // create a new display List
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
  //dplist_index = glGenLists( 1 ); // create a new display List
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
//  
void GLTextObject::setText(const QString &p_label,const QString& p_text)
{
  label = p_label;
  text  = p_text;
}
// ============================================================================
//  
void GLTextObject::setFont(QFont &f)
{
  font = f;
}
// ============================================================================
//  
int GLTextObject::getLabelWidth()
{
  QFontMetrics fm(font);
  return (fm.width(label,label.length()));
}
// ============================================================================
//  
int GLTextObject::getTextWidth()
{
  QFontMetrics fm(font);
  return (fm.width(text,text.length()));
}
// ============================================================================
//  
int GLTextObject::getHeight()
{
  QFontMetrics fm(font);
  return (fm.height());
}
// ============================================================================
//  
#if 0
void GLTextObject::buildDisplayList(const int x,const int y)
{
  // display list
  glNewList( dplist_index, GL_COMPILE );
    renderText(x,y,text,font);
    cerr << dplist_index << ": [" << text << "]\n";
  glEndList();
}
#endif
// ============================================================================
//
void GLTextObject::setPos(const int new_x, const int new_y, 
                          const int new_x_text)
{
  x = new_x;
  y = new_y;
  x_text = new_x_text;
}
// ============================================================================
//
void GLTextObject::display(GLBox * glbox)
{
  if (is_activated) {   
    glColor4ub(mycolor.red(), mycolor.green(), mycolor.blue(),255);
    // It's not possible to use renderText method inside a display List,
    // because 'renderText' already use a display list. It drove me nuts
    // during a while :(.
    glbox->renderText(x,y,label,font);
    glbox->renderText(x_text,y,text,font);
  }
}
//
