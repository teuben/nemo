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
// GLBox class implementation                                                  
//                                                                             
// Manage all OpenGL stuff                                                     
// ============================================================================
#include <iostream>
#include <iomanip>
#include <qfont.h>
#include <math.h> // happy gcc 3.4.x
#include "glbox.h"
#include "globjwin.h"

#define LOCAL_DEBUG 0
#include "print_debug.h"

#define POLY 1
using namespace std;

// Initialize static class variables:
int GLBox::width=894, GLBox::height=633;

/*!
  Create a GLBox widget
*/
// ----------------------------------------------------------------------------
// constructor
GLBox::GLBox( QWidget* parent, const char* name, const QGLWidget* shareWidget,
            const bool _blending, const bool _dbuffer, const bool _grid, 
            const float _psize)
  : QGLWidget( parent, name, shareWidget )
{
  xRot = yRot = zRot = 0.0;       // default object rotation
  xTrans = yTrans = zTrans = 0.0; // default object translation
  ixTrans = iyTrans = izTrans = 0;
  scale = 1.0;                    // default object scale
  enable_s = false;
  makeCurrent();   // activate OpenGL context, can build display list by now

  nb_object = 0;   // no object so far

  // grid
  show_grid=_grid;
  gridx = new GLGridObject(0,QColor(136,141,102));  //yellow);
  gridy = new GLGridObject(1,red,FALSE);
  gridz = new GLGridObject(2,blue,FALSE);
  line_aliased = FALSE;
  // blending
  particles_size = _psize;
  blending = _blending;
  depth_buffer = _dbuffer;
  MAX_PARTICLES_SIZE=5.0;
  // poly
  show_poly = FALSE;
  // HUD
  //QFont f=QFont("Courier", 12, QFont::Light);
  QFont f;
  f.setFamily("fixed");
  f.setRawMode(true);
  //f.setStyleStrategy(QFont::OpenGLCompatible);
  f.setPixelSize(10);
  f.setFixedPitch ( true )  ;
  f.setStyleHint(QFont::AnyStyle, QFont::PreferBitmap);  
  
  hud = new GLHudObject(width,height,f,yellow);
}
// ----------------------------------------------------------------------------
// Release allocated resources
GLBox::~GLBox()
{
  //makeCurrent();				// We're going to do gl calls
  delete gridx;
  delete gridy;
  delete gridz;
}
// ----------------------------------------------------------------------------
//   Set up the OpenGL rendering state. Robustly access shared display list.
void GLBox::initializeGL()
{
  // Let OpenGL clear to black
  qglClearColor( black ); 

  glEnable(GL_DEPTH_TEST);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  zoom = -10.;
  
  // Enable GL textures
  //glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
  //glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );

  // Nice texture coordinate interpolation
  glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );
  std::cerr << "-- Initialize GL --\n";
#if POLY  
  //glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
  loadImage();
#endif  
}
// ----------------------------------------------------------------------------
// collect Data and fill DisplayList
void GLBox::getData(const int * nbody,const float * pos,const ParticlesRangeVector  * prv)
{
  // delete previous objects
#if 0
  for (int i=0; i<nb_object; i++) {
    delete vparticles_object[i];
  }
  nb_object=0;
#endif
#if 0
  std::cerr << "void GLBox::getData, nbody = " << *nbody << "\n";
  std::cerr << "void GLBox::getData, prv-size()=" << prv->size() << "\n";
  std::cerr << "void GLBox::getData, nobject = " << nb_object << "\n";
#endif
  
  //vparticles_object = new GLParticlesObject[prv->size()];
  //nb_object=0;
  for (int i=0; i< (int ) prv->size(); i++ ) {
    //(*prv)[i].printRange();
    if (i>=nb_object ) {
      nb_object++;
      GLParticlesObject * pobj;
      pobj = new GLParticlesObject( nbody, pos,
      				 &((*prv)[i]));
     
      //vparticles_object.push_back(*pobj);
      vparticles_object[i] = pobj;
    }
    else { // update the current object
#if 1    
      vparticles_object[i]->updateObject( nbody, pos,
					  &((*prv)[i]));
#endif                                          
    }
  }
  // desactivate not selected objects
  for (int i=(int ) prv->size(); i< nb_object; i++ ) {
    vparticles_object[i]->setActivate(FALSE);
  }
  updateGL();
}
// ----------------------------------------------------------------------------
// setup zoom according to a z value
void GLBox::setZoom(float z)
{
  zoom = z;
  //cerr << "Zoom in setZoom = " << zoom << "\n";
  hud->setText(GLHudObject::Zoom,(const float) zoom);
  hud->updateDisplay();  
  updateGL();
}
// ----------------------------------------------------------------------------
// setup zoom according to wheel mouse
void GLBox::setZoom(int z)
{
  if (z>0) {
    zoom *= 1.025;//1.1;
  } else {
    zoom *= 0.975;//0.9;
  }
  //cerr << "zoom = " << zoom << "\n";
  //setTranslation();
  hud->setText(GLHudObject::Zoom,(const float) zoom);
  hud->updateDisplay();
  updateGL();
}
// ----------------------------------------------------------------------------
// Get Zoom
float GLBox::getZoom()
{
  return zoom;
}
// ----------------------------------------------------------------------------
// set Width and Height window property
void GLBox::setWH(int w, int h)
{
  width  = w;
  height = h;
  hud->setWH(w,h);
  updateGL();
}

// ----------------------------------------------------------------------------
// set projection according to Width and Height window property
void GLBox::setProjection(int width, int height)
{
  glViewport( 0, 0, width, height);
  ratio =  ((double )width) / ((double )height);
  //std::cerr << "ration = " << ratio << "\n";
  glMatrixMode(GL_PROJECTION); 
  glLoadIdentity();
  //glOrtho(-10.0, 10.0, -10.0, 10.0, 1000, -1000);
  gluPerspective(45.,ratio,0.05,10000.0);
  setViewPort();
  //std::cerr << "Width = " << width << " Height = " << height <<"\n";

}
// ----------------------------------------------------------------------------
//Paint the box. The actual openGL commands for drawing the box are
//performed here.
void GLBox::paintGL()
{
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  setProjection(width, height);
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();

  glEnable(GL_DEPTH_TEST);
  //std::cerr << "zoom = " << zoom << "object = " << object << "\n";
  //std::cerr << "In PaintGL\n";
  //glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mMod);
  // translate camera
#if 1  
  glTranslatef( 0.0, 0.0, zoom);  
#else  
  gluLookAt(0.0,  0.0,  zoom,
            0.0,  0.0,  0.0 ,
            0.0,  -1.0,  0.0);
  //glScalef( scale, scale, scale );
#endif

  glRotatef( xRot, 1.0, 0.0, 0.0 ); 
  glRotatef( yRot, 0.0, 1.0, 0.0 ); 
  glRotatef( zRot, 0.0, 0.0, 1.0 );
  
  // Grid Anti aliasing
  if (line_aliased) {
    glEnable(GL_LINE_SMOOTH);
  } else {
    glDisable(GL_LINE_SMOOTH);
  }
  // Display grid
  if (show_grid) {
    gridx->display();
    gridy->display();
    gridz->display();
  }  
  // Translate particles
  glTranslatef( xTrans, yTrans, zTrans);  
  
  setModelMatrix(); // save ModelView Matrix
  setProjMatrix();  // save Projection Matrix

  //bestFit();
  
  // display particles object
  //glDisable(GL_DEPTH_TEST);

  glPointSize(particles_size);
  glEnable(GL_POINT_SMOOTH);

  // control blending on particles
  if (blending) {
    glBlendFunc( GL_SRC_ALPHA, GL_ONE );  
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //glBlendFunc(GL_DST_ALPHA, GL_ONE);
    glEnable(GL_BLEND);  
  } else {
    glDisable(GL_BLEND);
  }
  // control depht buffer on particles
  if (depth_buffer) {
    glEnable(GL_DEPTH_TEST);
  } else {
    glDisable(GL_DEPTH_TEST);
  }

  // Display Particles
  for (int i=0; i<nb_object; i++) {
    vparticles_object[i]->display();
  }  

#if POLY  
  if (show_poly) { // Display Polygons
    glLoadIdentity(); // reset opengl state machine
    for (int i=0; i<nb_object; i++) {
      vparticles_object[i]->displayPolygons(mModel,texture[0],u_max,v_max);
    }
  }
#endif  
  glEnable(GL_DEPTH_TEST); // in any case enable depth buffer for grids

  // Try to print some text
  hud->display(this);
#if 0
  if (enable_s)  
    screenshot();
#endif  
}
// ----------------------------------------------------------------------------
// Set up the OpenGL view port, matrix mode, etc.
void GLBox::resizeGL( int w, int h )
{  
  glViewport( 0, 0, (GLint)w, (GLint)h );
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glFrustum(-1.0, 1.0, -1.0, 1.0, 1.0, 1000.0);
  hud->setWH(w,h);
}
// ----------------------------------------------------------------------------
//  Set the rotation angle of the object to n degrees around the X,Y,Z axis.
void GLBox::setRotation( int x, int y, int z )
{
  xRot = (GLfloat)(x % 360);
  yRot = (GLfloat)(y % 360);
  zRot = (GLfloat)(z % 360);
  //cerr << "xrot=" << xRot << " yrot=" << yRot << "\n";
  hud->setText(GLHudObject::Rot,xRot,yRot,zRot);
  hud->updateDisplay();
  updateGL();
}
// ----------------------------------------------------------------------------
//  Set the rotation angle of the object to \e degrees around the X axis.
void GLBox::setXRotation( int degrees )
{
  xRot = (GLfloat)(degrees % 360);
  updateGL();
}
// ----------------------------------------------------------------------------
//  Set the rotation angle of the object to \e degrees around the Y axis.
void GLBox::setYRotation( int degrees )
{
  yRot = (GLfloat)(degrees % 360);
  updateGL();
}
// ----------------------------------------------------------------------------
//  Set the rotation angle of the object to \e degrees around the Z axis.
void GLBox::setZRotation( int degrees )
{
  zRot = (GLfloat)(degrees % 360);
  updateGL();
}
// ----------------------------------------------------------------------------
// return x,y,z pixel translation according to current translation
// viewport coordinates and zoom
void GLBox::getPixelTranslation(int *x, int *y, int *z)
{
  GLint	Viewport[4];
  glGetIntegerv(GL_VIEWPORT,Viewport);
  *x = (int) (-xTrans * Viewport[2]/zoom);
  *y = (int) ( yTrans * Viewport[3]/zoom);
  *z = (int) ( zTrans * Viewport[2]/zoom);
  //cerr << "GPT : " << *x << " " << *y << " " << *z << "\n";
}
// ----------------------------------------------------------------------------
//  Set the translation angle 
void GLBox::setTranslation()
{
  setTranslation(ixTrans, iyTrans, izTrans);
}
// ----------------------------------------------------------------------------
//  Set the translation angle 
void GLBox::setTranslation( int x, int y, int z )
{
  GLint	Viewport[4];
  glGetIntegerv(GL_VIEWPORT,Viewport);
  
  ixTrans = x;
  iyTrans = y;
  izTrans = z;
  
  xTrans = (GLfloat)(-x*zoom/(Viewport[2]));
  yTrans = (GLfloat)(y*zoom/(Viewport[3])); //Viewport[3]*ratio));
  zTrans = (GLfloat)(z*zoom/(Viewport[2]));
  
  //cerr << "xrot=" << xRot << " yrot=" << yRot << "\n";
  hud->setText(GLHudObject::Trans,xTrans,yTrans,zTrans);
  hud->updateDisplay();
  updateGL();
}
// ----------------------------------------------------------------------------
//  toggle grid
void GLBox::toggleGrid( )
{
  show_grid = ! show_grid;
  updateGL();
}
// ----------------------------------------------------------------------------
//  setParticlesSize
void GLBox::setParticlesSize(int value)
{
  particles_size = 1.0 + ((GLfloat) value) * MAX_PARTICLES_SIZE / 100.0;
  //cerr << "Particles size =" << particles_size << "\n"; 
  updateGL();
}
// ----------------------------------------------------------------------------
//  
void GLBox::setHud(const GLHudObject::HudKeys k, const QString text)
{
  hud->setText(k,text);
  updateGL();
}
// ----------------------------------------------------------------------------
//  
void GLBox::setHud(const GLHudObject::HudKeys k, const float value)
{
  hud->setText(k,(const float) value);
  updateGL();
}
// ----------------------------------------------------------------------------
//  
void GLBox::setHud(const GLHudObject::HudKeys k, const int value)
{
  hud->setText(k,(const int) value);
  updateGL();
}
// ----------------------------------------------------------------------------
//  
void GLBox::setHudToggle(const GLHudObject::HudKeys k)
{
  hud->keysToggle(k);
  updateGL();
}
// ----------------------------------------------------------------------------
//  
void GLBox::setHudActivate(const GLHudObject::HudKeys k, const bool status)
{
  hud->keysActivate(k,status);
  updateGL();
}
void GLBox::bestFit()
{
  
  GLdouble mProj[16],mMod[16];
  glGetDoublev(GL_PROJECTION_MATRIX, (GLdouble *) mProj);
  glGetDoublev(GL_MODELVIEW_MATRIX, (GLdouble *) mMod);
  //float cc[4] = { 9.87903 ,10.5217 ,10.9075 , 1};
  float cc[4] = { -56.20 ,13.0 ,-28.0 , 1};
  
  float T1[4] = { 0.,0.,0.,0.};
  float T2[4] = { 0.,0.,0.,0.};
#define MT(row,col)  mTrans[col*4+row]
#define MP(row,col)  mProj[col*4+row]
#define MM(row,col)  mMod[col*4+row]

#if 1
  std::cerr << "ModelView Matrix:\n";
  for (int i=0;i<4;i++) {
    for (int j=0;j<4;j++) {
       T1[i] += (MM(i,j) * cc[j]);
       std::cerr << setw(10) << MM(i,j) <<"\t";
    }
    std::cerr << "\n";
  }
  std::cerr << "Vector T1:\n";
  for (int j=0;j<4;j++) {
      std::cerr << setw(10) << T1[j] <<"\t";
  }
  std::cerr << "\n";
  
  std::cerr << "Projection Matrix:\n";
  for (int i=0;i<4;i++) {
    for (int j=0;j<4;j++) {
       T2[i] += (MP(i,j) * T1[j]);
       std::cerr << setw(10) << MP(i,j) <<"\t";
    }
    std::cerr << "\n";
  } 
  std::cerr << "Vector T2:\n";
  for (int j=0;j<4;j++) {
      std::cerr << setw(10) << T2[j] <<"\t";
  }
  std::cerr << "\n";  
#endif
  T2[0] /= T2[3];
  T2[1] /= T2[3];
  T2[2] /= T2[3];
  std::cerr << "\n$$$$$$$$$$$$$\n";
  
  int view[4];
  glGetIntegerv(GL_VIEWPORT,view);
  std::cerr << "ViewPort = " << view[0] << " " << view[1] << "\n";
  std::cerr << "ViewPort = " << view[2] << " " << view[3] << "\n";
  GLdouble winX,winY,winZ;
  gluProject(cc[0],cc[1],cc[2],mMod,mProj,view,&winX,&winY,&winZ);
  float winX2=view[0] + view[2] * (T2[0]+1.)/2.;
  float winY2=view[1] + view[3] * (T2[1]+1.)/2.;
  std::cerr << "Winx = "<< winX << " Winy = "<<winY<<"\n";
  std::cerr << "Winx2 = "<< winX2 << " Winy2 = "<<winY2<<"\n";
  std::cerr << "zoomx = " << -cc[2] + MP(0,0) * cc[0] * view[2]/ (MP(3,2)*(2.*(0.-view[0])-view[2])) <<"\n";
    //std::cerr << cc[i]*zoom/(zoom-cc[2])<<"\t";
//  }
  std::cerr << "\n";
    
}
void GLBox::screenshot() {
    //updateGL();
    swapBuffers();
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
  
    QImage img(viewport[2],viewport[3],32);
    GLfloat *r=new GLfloat[viewport[2]];
    GLfloat *g=new GLfloat[viewport[2]];
    GLfloat *b=new GLfloat[viewport[2]];
    
    for(int y=0;y<viewport[3];y++) {
      glReadPixels(viewport[0],y,viewport[2],1,GL_RED,GL_FLOAT,(void *) r);
      glReadPixels(viewport[0],y,viewport[2],1,GL_GREEN,GL_FLOAT,(void *) g);
      glReadPixels(viewport[0],y,viewport[2],1,GL_BLUE,GL_FLOAT,(void *) b);
      
      for(int x=0;x<viewport[2];x++) {
        uint *p = (uint *)img.scanLine(viewport[3]-1-y)+x;
        *p = qRgb((uint) (r[x]*255),(uint) (g[x]*255),(uint) (b[x]*255));    
      }
  }
    QImageIO iio;
  iio.setImage(img);
    
  static int no=0;
  char shot[100];
  sprintf(shot,"frame.%05d.png",no);
  //takeScreenshot(shot);
  std::cerr <<  ">> shot = " << shot << "\n";
  no++;
  
  //iio.setFileName("screenshot.bmp");
  iio.setFileName(shot);
  iio.setFormat("PNG");
  iio.write();
  
  delete[] r;
  delete[] g;
  delete[] b;
}
// ============================================================================
// 
void GLBox::loadImage()
{
  //QString name = QFileDialog::getOpenFileName(".", "Images (*.png *.xpm *.jpg)", this, "Choose", "Select an image");
  
  // In case of Cancel
  //if (name.isEmpty())
  //  return;
  //QString name="/home/jcl/download/PEngine/images/particle.png";
  QString name="/home/jcl/download/qq/smoke.png";
  
  QImage img(name);
  
  if (img.isNull())
    {
      std::cerr << "Unable to load " << name << ", unsupported file format" << std::endl;
      return;
    }
    

  std::cout << "Loading " << name << ", " << img.width() << "x" << img.height() << " pixels" << std::endl;

  // 1E-3 needed. Just try with width=128 and see !
  int newWidth  = 1<<(int)(1+log(img.width() -1+1E-3) / log(2.0));
  int newHeight = 1<<(int)(1+log(img.height()-1+1E-3) / log(2.0));

  u_max = img.width()  / (float)newWidth;
  v_max = img.height() / (float)newHeight;
  
  if ((img.width()!=newWidth) || (img.height()!=newHeight))
    {
      std::cout << "Image size set to " << newWidth << "x" << newHeight << " pixels" << std::endl;
      img = img.copy(0, 0, newWidth, newHeight);
    }

  tratio = newWidth / (float)newHeight;
  
  QImage glImg = QGLWidget::convertToGLFormat(img);  // flipped 32bit RGBA
  glEnable(GL_TEXTURE_2D);  
  glGenTextures(1, &texture[0]);					// Create The Texture
  
  // Typical Texture Generation Using Data From The Bitmap
  glBindTexture(GL_TEXTURE_2D, texture[0]);
  // Bind the img texture...
  glTexImage2D(GL_TEXTURE_2D, 0, 4, glImg.width(), glImg.height(), 0,
	       GL_RGBA, GL_UNSIGNED_BYTE, glImg.bits());
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
  //glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
  //glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
  glDisable(GL_TEXTURE_2D);
}
//
