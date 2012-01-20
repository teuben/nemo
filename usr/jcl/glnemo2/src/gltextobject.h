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
//                                                                             
// GLTextObject class definition                                               
//                                                                             
// Manage OpenGL Text Object used on the On Screen Display Display             
// ============================================================================
#ifndef GL_TEXT_OBJECT_H
#define GL_TEXT_OBJECT_H

#include <qgl.h>
#include "globject.h"
#include <iostream>
#include "fnt.h"
namespace glnemo { 
class GLWindow;

using namespace std;

class GLTextObject : public GLObject {

  public:
  GLTextObject(const QString text,const fntRenderer &f,const QColor &c=Qt::green, bool activated=TRUE);
  GLTextObject(const fntRenderer &f,const QColor &c=Qt::green, 
               bool activated=TRUE);
  GLTextObject( bool activated=TRUE);
    
  ~GLTextObject();
  
  void setText(const QString &p_label,const QString &p_text);
  void setFont(fntRenderer &f);
  int getLabelWidth();
  int getTextWidth();
  int getHeight();
  void setPos(const int,const  int, const int);
  void display(const int width, const int height);
  private:
  // data
  QString label,text;
  fntRenderer font;
  int x,y;      // xy label text position
  int x_text;   // x offset text position
};
} // namespace
#endif
// ============================================================================
