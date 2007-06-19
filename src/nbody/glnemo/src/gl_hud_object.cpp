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
// GLHudObject class implementation                                            
//                                                                             
// Manage OpenGL HUD (Head Up Display) Object                                  
// ============================================================================


#include "gl_hud_object.h"
#include <iostream>
using namespace std;

// init static data
char * GLHudObject::HUDText[n_HudKeys] = {
    ""                   , // DataType
    ""                   , // Title
    "nbody"              , // Nbody
    "time"               , // Time
    ""                   , // Getdata
    "Zoom"               , // Zoom
    "Rot"                , // Rotation Angle
    "Trans"              , // Tanslation
    ""                   , // Loading
    ""                     // Projection type
};

#define MAX(x,y) ((x) > (y) ? (x) : (y))
// ============================================================================
// constructor                                                                 
GLHudObject::GLHudObject(const int w, const int h,
                         const QFont &f,const QColor &c, 
                         bool activated):GLObject()
{
  mycolor = c;
  is_activated = activated;
  font = f;
  width = w;
  height = h;
  
  //hud_text = new GLTextObject[n_HudKeys](font,c);
  hud_text = new GLTextObject[n_HudKeys];
  for (int i=0; i<n_HudKeys; i++) {
    hud_text[i].setFont(font);
    hud_text[i].setColor(c);
  }
  setTextColor(GLHudObject::Loading,red);
  // initialization
  for (int i=0; i<n_HudKeys; i++) {
    setText((const HudKeys) i,"");
  }
  setText(Loading,"Loading....");
  setText(Trans,0.,0.,0.);
  setText(Rot,0.,0.,0.);
  setText(Zoom,(const float) 0.);
  hud_text[Loading].setActivate(FALSE);
  updateDisplay();
}
// ============================================================================
// destructor                                                                  
GLHudObject::~GLHudObject()
{
}
// ============================================================================
// GLHudObject::setWH()                                                        
// set witdh and height window                                                 
void GLHudObject::setWH(int w, int h)
{
  width  = w;
  height = h;
  updateDisplay();
}
// ============================================================================
// GLHudObject::setText()                                                      
// set Text to the selected HubObject                                          
void GLHudObject::setText(const HudKeys k, const QString text)
{
  hud_text[k].setText(HUDText[k], text);
  updateDisplay((const HudKeys) k);
}
// ============================================================================
// GLHudObject::setText()                                                      
// set int value to the selected HubObject                                     
void GLHudObject::setText(const HudKeys k, int value)
{
  setText(k,QString(": %1").arg(value));  
}
// ============================================================================
// GLHudObject::setText()                                                      
// set float value to the selected HubObject                                   
void GLHudObject::setText(const HudKeys k, float value)
{
  setText(k,QString(": %1").arg(value,0,'f',4));
}  
// ============================================================================
// GLHudObject::setText()                                                      
// set triple float value to the selected HubObject                            
void GLHudObject::setText(const HudKeys k,const float x, const float y,
                          const float z)
{
QString my_text=QString(": %1 %2 %3").arg(x,0,'f',2)
                                     .arg(y,0,'f',2)
                                     .arg(z,0,'f',2);             
  setText(k,my_text);
}  
// ============================================================================
// GLHudObject::setTextColor()                                                 
// set text color to the selected HubObject                                    
void GLHudObject::setTextColor(const HudKeys k, const QColor c)
{
  hud_text[k].setColor(c);
}
// ============================================================================
// GLHudObject::keysToggle()                                                   
// Toggle  the selected HubObject                                              
void GLHudObject::keysToggle(const HudKeys k)
{
  hud_text[k].toggleActivate();
}
// ============================================================================
// GLHudObject::keysActivate()                                                 
// Activate the selected HubObject                                             
void GLHudObject::keysActivate(const HudKeys k, const bool status)
{
  hud_text[k].setActivate(status);
}
// ============================================================================
// GLHudObject::setFont()                                                      
// set global font                                                             
void GLHudObject::setFont(QFont f)
{
  font =f;
}
// ============================================================================
// GLHudObject::setFont()                                                      
// set font to the selected HubObject                                          
void GLHudObject::setFont(const HudKeys k,QFont f)
{
  hud_text[k].setFont(f);
}
// ============================================================================
// GLHudObject::updateDisplay()                                                
// Update display for all HudObject activated                                  
void GLHudObject::updateDisplay()
{
  if (is_activated) {
    for (int i=0; i<n_HudKeys; i++) {
      if (hud_text[i].getActivate()) {
        updateDisplay((const HudKeys) i);
      }
    }
  }
}
// ============================================================================
// GLHudObject::updateDisplay()                                                
// Update display for HudObject selected                                       
void GLHudObject::updateDisplay(const HudKeys k)
{
  int x=0,y=0,x_text=0, max;
  QFontMetrics fm(font);
  
  switch (k) {
    case DataType:
      x=0;
      y=hud_text[k].getHeight();                                   
      break;
    case Title:
      x=0;
      x_text=(width/2)-hud_text[k].getTextWidth()/2;
      y=hud_text[k].getHeight();
      break;
    case Nbody:
      max=MAX(hud_text[Nbody].getLabelWidth(),
              hud_text[Time].getLabelWidth());
      x = 0;
      y = hud_text[Time].getHeight()+2+hud_text[Nbody].getHeight();
      x_text = max;
      break;
    case Time:
      max=MAX(hud_text[Nbody].getLabelWidth(),
              hud_text[Time].getLabelWidth());  
      x = 0;
      y = hud_text[Time].getHeight();
      x_text = max;
      break;
    case Getdata:
      x = 0;
      x_text=width-hud_text[k].getTextWidth()-3;
      y=hud_text[k].getHeight();
      break;
    case Zoom:
      x=0;
      y=height-5-hud_text[Trans].getHeight()-1-
               hud_text[Rot].getHeight()-1;
      x_text = hud_text[Trans].getLabelWidth();
      break;
    case Rot:
      x=0;
      y=height-5-hud_text[Trans].getHeight()-1;
      x_text = hud_text[Trans].getLabelWidth();
      break;
    case Trans:
      x=0;
      y=height-5;
      x_text = hud_text[Trans].getLabelWidth();
      break;
    case Loading:
      x = 0;
      x_text=width-hud_text[k].getTextWidth()-3;
      y=hud_text[Getdata].getHeight()+2+hud_text[k].getHeight();    
      break;  
    case Projection:
      x=0;
      x_text=(width/2)-hud_text[k].getTextWidth()/2;
      y=height-5;
      break;          
    case n_HudKeys:
      break;      
  } 
  if (k < n_HudKeys ) {
    hud_text[k].setPos(x,y, x_text); 
  }
}
// ============================================================================
// GLHudObject::display()                                                      
// render hud text object                                                      
void GLHudObject::display(GLBox * gg)
{
  if (is_activated) {
    for (int i=0; i<n_HudKeys; i++) {
      if (hud_text[i].getActivate()) {
        hud_text[(const HudKeys) i].display(gg);
      }
    }
  }  
}
// ============================================================================
// GLHudObject::updateColor()                                                  
// render hud text object                                                      
void GLHudObject::updateColor(const QColor col)
{
  if (is_activated) {
    for (int i=0; i<n_HudKeys; i++) {
        hud_text[(const HudKeys) i].setColor(col);
    }
  }  
}
// ============================================================================

