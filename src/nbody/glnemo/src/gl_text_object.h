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
// GLTextObject class definition                                               
//                                                                             
// Manage OpenGL Text Object used in the Head Up Display                       
// ============================================================================
#ifndef GL_TEXT_OBJECT_H
#define GL_TEXT_OBJECT_H

#include <qgl.h>
#include "gl_object.h"
#include <iostream>
class GLBox;
using namespace std;

class GLTextObject : public GLObject {

  public:
  GLTextObject(const QString text,const QFont &f,const QColor &c=green, bool activated=TRUE);
  GLTextObject(const QFont &f,const QColor &c=green, 
               bool activated=TRUE);
  GLTextObject( bool activated=TRUE);
    
  ~GLTextObject();
  
  void setText(const QString &p_label,const QString &p_text);
  void setFont(QFont &f);
  int getLabelWidth();
  int getTextWidth();
  int getHeight();
  void setPos(const int,const  int, const int);
  void display(GLBox * gg);
  private:
  // data
  QString label,text;
  QFont font;
  int x,y;      // xy label text position
  int x_text;   // x offset text position
};
#endif
// ============================================================================
