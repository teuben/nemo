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
// GLHudObject class definition                                                
//                                                                             
// Manage OpenGL Hud Object                                                    
// ============================================================================
#ifndef GL_HUD_OBJECT_H
#define GL_HUD_OBJECT_H

#include <qgl.h>
#include "gl_object.h"
class GLBox;
class GLHudObject;
#include "gl_text_object.h"


class GLHudObject : public GLObject {
  Q_OBJECT
  
  public:
  GLHudObject(const int w, const int h,const QFont &f,
              const QColor &c=green, bool activated=TRUE);
  ~GLHudObject();
  enum HudKeys {
    DataType  , 
    Title     , 
    Nbody     , 
    Time      , 
    Getdata   , 
    Zoom      , 
    Rot       , 
    Trans     ,
    Loading   ,
    Projection,
    n_HudKeys
  };
  
  //GLTextObject hud_text[n_HudKeys];
  GLTextObject * hud_text;
  
  public slots:
  void setWH(int width, int height);  
  void setFont(const QFont);
  void setFont(const HudKeys k, const QFont);
  void setText(const HudKeys k, const QString text);
  void setText(const HudKeys k, const int);
  void setText(const HudKeys k, const float);
  void setText(const HudKeys k, const float, const float, const float);
  void setTextColor(const HudKeys k,const QColor );
  void keysToggle(const HudKeys k);
  void keysActivate(const HudKeys k, const bool status);
  void updateDisplay();
  void updateDisplay(const HudKeys k);
  void display(GLBox *);
  void updateColor(const QColor);
  private:
  static char * HUDText[n_HudKeys];
  QFont font;
};
#endif
// ============================================================================

