// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2015                                  
// e-mail:   Jean-Charles.Lambert@lam.fr                                      
// address:  Centre de donneeS Astrophysique de Marseille (CeSAM)              
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include "globjectosd.h"
#include <GL/glu.h>

namespace glnemo {
// init static data
char * GLObjectOsd::OsdText[n_OsdKeys] = {
 (char *) ""                   , // DataType
 (char *) ""                   , // Title
 (char *) "nbody"              , // Nbody
 (char *) "time"               , // Time
 (char *) ""                   , // Getdata
 (char *) "Zoom"               , // Zoom
 (char *) "Rot"                , // Rotation Angle
 (char *) "Trans"              , // Tanslation
 (char *) ""                   , // Loading
 (char *) ""                     // Projection type
};

#define MAX(x,y) ((x) > (y) ? (x) : (y))

// ============================================================================
// constructor                                                                 
GLObjectOsd::GLObjectOsd(const int w, const int h,
			 const fntRenderer &f,const QColor &c,
			 bool activated):GLObject()
{
  mycolor = c;
  is_activated = activated;
  font = f;
  width = w;
  height = h;
  
  //Osd_text = new GLTextObject[n_OsdKeys](font,c);
  Osd_text = new GLTextObject[n_OsdKeys];
  for (int i=0; i<n_OsdKeys; i++) {
    Osd_text[i].setFont(font);
    Osd_text[i].setColor(c);
    Osd_text[i].setActivate(true);
  }
  setTextColor(GLObjectOsd::Loading,Qt::red);
  // initialization
  for (int i=0; i<n_OsdKeys; i++) {
    setText((const OsdKeys) i,"");
  }
  setText(Loading,"Loading....");
  setText(Trans,0.,0.,0.);
  setText(Rot,0.,0.,0.);
  setText(Zoom,(const float) 0.);
  setText(Projection,"Perspective");
  Osd_text[Loading].setActivate(FALSE);
  updateDisplay();
}
// ============================================================================
// destructor                                                                  
GLObjectOsd::~GLObjectOsd()
{
}
// ============================================================================
// GLObjectOsd::setWH()
// set witdh and height window                                                 
void GLObjectOsd::setWH(int w, int h)
{
  width  = w;
  height = h;
  updateDisplay();
}
// ============================================================================
// GLObjectOsd::setText()
// set Text to the selected HubObject                                          
void GLObjectOsd::setText(const OsdKeys k, const QString text)
{
  Osd_text[k].setText(OsdText[k], text);
  updateDisplay((const OsdKeys) k);
}
// ============================================================================
// GLObjectOsd::setText()
// set int value to the selected HubObject                                     
void GLObjectOsd::setText(const OsdKeys k, int value)
{
  setText(k,QString(": %1").arg(value));
}
// ============================================================================
// GLObjectOsd::setText()
// set float value to the selected HubObject                                   
void GLObjectOsd::setText(const OsdKeys k, float value)
{
  setText(k,QString(": %1").arg(value,0,'f',4));
}
// ============================================================================
// GLObjectOsd::setText()
// set triple float value to the selected HubObject                            
void GLObjectOsd::setText(const OsdKeys k,const float x, const float y,
			  const float z)
{
  QString my_text=QString(": %1 %2 %3").arg(x,0,'f',2)
      .arg(y,0,'f',2)
      .arg(z,0,'f',2);
  setText(k,my_text);
}
// ============================================================================
// GLObjectOsd::setTextColor()
// set text color to the selected HubObject                                    
void GLObjectOsd::setTextColor(const OsdKeys k, const QColor c)
{
  Osd_text[k].setColor(c);
}
// ============================================================================
// GLObjectOsd::setColor()
// set text color to the selected HubObject                                    
void GLObjectOsd::setColor(const QColor c)
{
  for (int i=0; i<n_OsdKeys; i++) {
    Osd_text[i].setColor(c);
  }
}
// ============================================================================
// GLObjectOsd::keysToggle()
// Toggle  the selected HubObject                                              
void GLObjectOsd::keysToggle(const OsdKeys k)
{
  Osd_text[k].toggleActivate();
}
// ============================================================================
// GLObjectOsd::keysActivate()
// Activate the selected HubObject                                             
void GLObjectOsd::keysActivate(const OsdKeys k, const bool status)
{
  Osd_text[k].setActivate(status);
}
// ============================================================================
// GLObjectOsd::setFont()
// set global font                                                             
void GLObjectOsd::setFont(fntRenderer f)
{
  font =f;
  // first we update the fonts
  for (int i=0; i<n_OsdKeys; i++) {
    Osd_text[i].setFont(f);
  }
  // 2nd we update text positions
  // first step must be complete otherwise
  // it crashs the application
  for (int i=0; i<n_OsdKeys; i++) {
    updateDisplay((const OsdKeys) i);
  }
}
// ============================================================================
// GLObjectOsd::setFont()
// set font to the selected HubObject                                          
void GLObjectOsd::setFont(const OsdKeys k,fntRenderer f)
{
  Osd_text[k].setFont(f);
}
// ============================================================================
// GLObjectOsd::updateDisplay()
// Update display for all OsdObject activated
void GLObjectOsd::updateDisplay()
{
  if (is_activated) {
    for (int i=0; i<n_OsdKeys; i++) {
      if (Osd_text[i].getActivate()) {
	updateDisplay((const OsdKeys) i);
      }
    }
  }
}
// ============================================================================
// GLObjectOsd::updateDisplay()
// Update display for OsdObject selected
void GLObjectOsd::updateDisplay(const OsdKeys k)
{
  int x=0,y=0,x_text=0, max;
  
  switch (k) {
    case DataType:
      x=0;
      y=Osd_text[k].getHeight();
      break;
    case Title:
      x=0;
      x_text=(width/2)-Osd_text[k].getTextWidth()/2;
      y=Osd_text[k].getHeight();
      break;
    case Nbody:
      max=MAX(Osd_text[Nbody].getLabelWidth(),
	      Osd_text[Time].getLabelWidth());
      x = 0;
      y = Osd_text[Time].getHeight()+2+Osd_text[Nbody].getHeight();
      x_text = max;
      break;
    case Time:
      max=MAX(Osd_text[Nbody].getLabelWidth(),
	      Osd_text[Time].getLabelWidth());
      x = 0;
      y = Osd_text[Time].getHeight();
      x_text = max;
      break;
    case Getdata:
      x = 0;
      x_text=width-Osd_text[k].getTextWidth()-3;
      y=Osd_text[k].getHeight();
      break;
    case Zoom:
      x=0;
      y=height-5-Osd_text[Trans].getHeight()-1-
	  Osd_text[Rot].getHeight()-1;
      x_text = Osd_text[Trans].getLabelWidth();
      break;
    case Rot:
      x=0;
      y=height-5-Osd_text[Trans].getHeight()-1;
      x_text = Osd_text[Trans].getLabelWidth();
      break;
    case Trans:
      x=0;
      y=height-5;
      x_text = Osd_text[Trans].getLabelWidth();
      break;
    case Loading:
      x = 0;
      x_text=width-Osd_text[k].getTextWidth()-3;
      y=Osd_text[Getdata].getHeight()+2+Osd_text[k].getHeight();
      break;
    case Projection:
      x=0;
      x_text=(width/2)-Osd_text[k].getTextWidth()/2;
      y=height-5;
      break;
    case n_OsdKeys:
      break;
  }
  if (k < n_OsdKeys ) {
    Osd_text[k].setPos(x,y, x_text);
  }
}
// ============================================================================
// GLObjectOsd::display()
// render Osd text object
void GLObjectOsd::display()
{
  if (is_activated) {
    // save OpenGL state
    glDisable( GL_DEPTH_TEST );
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    //glOrtho(0.,width,0.,height,-1,1);
    gluOrtho2D(0.,width,0.,height);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    //glBlendFunc( GL_SRC_ALPHA, GL_ONE ); // original
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  // No Alpha bending accumulation
    glEnable(GL_BLEND);
    
    font.begin();
    
    for (int i=0; i<n_OsdKeys; i++) {
      if (Osd_text[i].getActivate()) {
	Osd_text[(const OsdKeys) i].display(width,height);
      }
    }
    
    font.end();
    glDisable(GL_BLEND);
    // Restore OpenGL state
    glMatrixMode( GL_PROJECTION );
    glPopMatrix();
    glMatrixMode( GL_MODELVIEW );
    glPopMatrix();
    glEnable( GL_DEPTH_TEST );
  }
}
// ============================================================================
// GLObjectOsd::updateColor()
// render Osd text object
void GLObjectOsd::updateColor(const QColor col)
{
  if (is_activated) {
    for (int i=0; i<n_OsdKeys; i++) {
      Osd_text[(const OsdKeys) i].setColor(col);
    }
  }
}
// ============================================================================




}
