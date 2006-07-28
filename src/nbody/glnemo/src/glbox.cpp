// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2005                                  
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
#include <string.h>
#include <qdatetime.h>
#include "glbox.h"

#include "smoke.h"
#define LOCAL_DEBUG 0
#include "print_debug.h"

#define POLY 1
using namespace std;

#define DOF 4000000
// Initialize static class variables:
int GLBox::width=894, GLBox::height=633;
float ortho_left,ortho_right,ortho_bottom,ortho_top;

// ============================================================================
// constructor                                                                 
GLBox::GLBox( QWidget* parent, const char* name,
            GlobalOptions * _options, AnimationEngine * _anim_engine,const QGLWidget* shareWidget)
  : QGLWidget( parent, name, shareWidget )
{
  xRot = yRot = zRot = 0.0;       // default object rotation
  xTrans = yTrans = zTrans = 0.0; // default object translation
  ixTrans = iyTrans = izTrans = 0;
  scale = 1.0;                    // default object scale
  enable_s = false;

  // get options
  store_options = _options;
  anim_engine   = _anim_engine; 
  makeCurrent();   // activate OpenGL context, can build display list by now

  nb_object = 0;   // no object so far

  // grid
  gridx = new GLGridObject(0,QColor(136,141,102),store_options->xy_grid);  //yellow);
  gridy = new GLGridObject(1,QColor(136,141,102),store_options->yz_grid);
  gridz = new GLGridObject(2,QColor(136,141,102),store_options->xz_grid);
  cube  = new GLCube(store_options->mesh_length*store_options->nb_meshs,
                     QColor(136,141,102),store_options->show_cube);
  resizeGrid(store_options->mesh_length,store_options->nb_meshs);
  line_aliased = FALSE;
  MAX_PARTICLES_SIZE=5.0;
  // poly
  show_poly = store_options->show_poly;
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
  hudProjection();
  tree = new GLOctree(_options);
}

// ============================================================================
// Destructor                                                                  
// Release allocated resources                                                 
GLBox::~GLBox()
{
  //makeCurrent();				// We're going to do gl calls
  delete gridx;
  delete gridy;
  delete gridz;
  delete cube;
}

// ============================================================================
// GLBox::initializeGL()                                                       
// Set up the OpenGL rendering state. Robustly access shared display list.     
void GLBox::initializeGL()
{
  qglClearColor( black );		// Let OpenGL clear to black
  glEnable(GL_DEPTH_TEST);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  store_options->zoom = -10.;
  
  // Enable GL textures
  //glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
  //glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );

  // Nice texture coordinate interpolation
  glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );
  PRINT_D std::cerr << "-- Initialize GL --\n";
#if POLY  
  //glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
  loadImage();
#endif  
}
// ============================================================================
// GLBox::paintGL()                                                            
// Paint the box. The actual openGL commands for drawing the box are           
// performed here.                                                             
void GLBox::paintGL()
{
#if 0
  static QTime t;
  static bool first=true;
  if (first) {
    first=false;
    t.start();
  } else {
    //t.restart();
    std::cerr << t.elapsed() << endl;
  }
#endif
  // record frame if activates
  anim_engine->record->beginFrame();
  
  qglClearColor( store_options->background_color);
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );  
  setProjection(width, height);
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();

  glEnable(GL_DEPTH_TEST);
  // translate camera
#if 1  
  if (store_options->perspective) {
    glTranslatef( 0.0, 0.0, store_options->zoom);  
  }

#else  
  gluLookAt(0.0,  0.0,  store_options->zoom,
            0.,0.,0.,
            0.,  1.0,  0.);
  //glScalef( scale, scale, scale );
#endif
  glRotatef( store_options->xrot, 1.0, 0.0, 0.0 ); 
  glRotatef( store_options->yrot, 0.0, 1.0, 0.0 ); 
  glRotatef( store_options->zrot, 0.0, 0.0, 1.0 );


#if 1  
  // Grid Anti aliasing
  if (line_aliased) {
    glEnable(GL_LINE_SMOOTH);
  } else {
    glDisable(GL_LINE_SMOOTH);
  }
  // Display grid
  if (store_options->show_grid) {
    gridx->display();
    gridy->display();
    gridz->display();
    cube->display();
  }  
#endif
  // Translate particles
  glTranslatef( store_options->xtrans, store_options->ytrans, store_options->ztrans);  
  
  setModelMatrix(); // save ModelView Matrix
  setProjMatrix();  // save Projection Matrix

#if 0  
  GLdouble xx,yy,zz;
  gluUnProject(width,0,0.5,mModel,mProj,viewport,&xx,&yy,&zz);
  std::cerr << "xx ="<<xx<<" yy="<<yy<<" zz="<<zz<<"\n";
#endif
  //bestFit();
  
  // display particles object
  //glDisable(GL_DEPTH_TEST);

  glPointSize(store_options->psize);
  glEnable(GL_POINT_SMOOTH);

  // octree display
  //treeUpdate();
  tree->display();
  
  // control blending on particles
  if (store_options->blending) {
    glBlendFunc( GL_SRC_ALPHA, GL_ONE );  
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //glBlendFunc(GL_DST_ALPHA, GL_ONE);
    glEnable(GL_BLEND);  
  } else {
    glDisable(GL_BLEND);
  }
#if 1
  // control depht buffer on particles
  if (store_options->dbuffer) {
    glEnable(GL_DEPTH_TEST);
  } else {
    glDisable(GL_DEPTH_TEST);
  }
#endif
  // Display Particles
  if (store_options->show_part) {
    for (int i=0; i<nb_object; i++) {
      vparticles_object[i]->updateAlphaSlot(store_options->particles_alpha);
      vparticles_object[i]->display();
    }  
  }

// display sprites
#if GL_EXT_ENABLE
if (0) {
  for (int i=0; i<nb_object; i++) {
      vparticles_object[i]->displaySprites(texture[0]);
  }
} 
else  
#endif
#if POLY  
  if (store_options->show_poly) { // Display Polygons          
    glLoadIdentity();             // reset opengl state machine
    for (int i=0; i<nb_object; i++) {
      vparticles_object[i]->displayPolygons(mModel,texture[0],u_max,v_max);
    }
  }
#endif  
  glEnable(GL_DEPTH_TEST); // in any case enable depth buffer for grids
  glDisable(GL_BLEND);     // disable blending for HUD                 
  // Head Up Display
  hud->display(this);
#if 0
  if (enable_s)  
    screenshot();
#endif  
  // record end of frame if activated
  anim_engine->record->endFrame(store_options);
}
// ============================================================================
// GLBox::resizeGL()                                                           
// Set up the OpenGL view port, matrix mode, etc.                              
void GLBox::resizeGL( int w, int h )
{  
  glViewport( 0, 0, (GLint)w, (GLint)h );
  ratio =  ((double )w) / ((double )h);
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glFrustum(-1.0, 1.0, -1.0, 1.0, 1.0, (float) DOF);
  hud->setWH(w,h);
  computeOrthoFactorRatio(); // recalculate Ortho projection
}
// ============================================================================
// GLBox::hudProjection()                                                      
// toogle projection type on the HUD                                           
void GLBox::hudProjection()
{
  if (store_options->perspective) {
    hud->setText(GLHudObject::Projection,"Perspective");  
  } 
  else {
    hud->setText(GLHudObject::Projection,"Orthographic");
  }
}
// ============================================================================
// GLBox::getData()                                                            
// collect Data and fill DisplayList                                           
void GLBox::getData(const int * nbody               ,  
                    const float * pos_              ,
                    ParticlesSelectVector  * psv)
{
  static int nbody_save=0;
  static float * pos=NULL;
  
  // copy positions into working array to prevent
  // bad display during "playing" snapshot       
  if (*nbody > nbody_save) {
    nbody_save=*nbody;
    if (pos) delete pos;
    pos = new float[nbody_save*3];
  }
  memcpy((float *) pos, (float *) pos_, sizeof(float)* (*nbody) * 3);
  
  
  for (int i=0; i< (int ) psv->size(); i++ ) {
    (*psv)[i].vps->defaultIndexTab();  // reset index, in case of loading
  }
  // create octree
  tree->update(nbody,pos,psv); 
  // delete previous objects
#if 0
  for (int i=0; i<nb_object; i++) {
    delete vparticles_object[i];
  }
  nb_object=0;
#endif
  
  // loop on all the objects stored in psv
  for (int i=0; i< (int ) psv->size(); i++ ) {
    (*psv)[i].vps->printRange();
    if (i>=nb_object ) { // create a new object
      nb_object++;
      GLParticlesObject * pobj;
      pobj = new GLParticlesObject( nbody, pos,
      				 ((*psv)[i].vps));
     
      //vparticles_object.push_back(*pobj);
      vparticles_object[i] = pobj;
    }
    else { // update the current object    
    #if 0
      (*psv)[i].vps->defaultIndexTab();  // reset index, in case of loading
    #endif  
      vparticles_object[i]->updateObject( nbody, pos,
					  ((*psv)[i].vps));        
    }
  }
  // desactivate not selected objects
  for (int i=(int ) psv->size(); i< nb_object; i++ ) {
    vparticles_object[i]->setActivate(FALSE);
  }
  setParticlesSize((int) store_options->psize,false);
  setTextureSize(store_options->texture_size,false);
  changeTextureAlphaColor(store_options->texture_alpha_color,false);
  // create octree
  //tree->update(nbody,pos,psv);  
#if 0  
  if (store_options->octree_enable) {
    for (int i=0; i<nb_object; i++) {
      vparticles_object[i]->rebuildDisplayList();
    }
  }
#endif  
  updateGL();
}
// ============================================================================
// GLBox::setZoom()                                                            
// setup zoom according to a z value                                           
void GLBox::setZoom(float z)
{
  zoom = z;
  //cerr << "Zoom in setZoom = " << zoom << "\n";
  hud->setText(GLHudObject::Zoom,(const float) zoom);
  hud->updateDisplay();  
  store_options->zoom=zoom;
  //store_options->zoomo=zoom;
  updateGL();
}
// ============================================================================
// GLBox::setZoom()                                                            
// setup zoom according to wheel mouse
void GLBox::setZoom(int z)
{
  if (z>0) {
    store_options->zoom *= 1.025; //1.1;
    store_options->zoomo *= 1.025;//1.1;
  } else {
    store_options->zoom *= 0.975; //0.9;
    store_options->zoomo *= 0.975;//0.9;
  }
  hud->setText(GLHudObject::Zoom,(const float) store_options->zoom);
  hud->updateDisplay();
  updateGL();
}
// ============================================================================
// GLBox::getZoom()                                                            
// return zoom value                                                           
float GLBox::getZoom()
{
  return store_options->zoom;
}
// ============================================================================
// GLBox::setWH()                                                              
// set Width and Height opengl window                                          
void GLBox::setWH(int w, int h)
{
  width  = w;
  height = h;
  hud->setWH(w,h);
  updateGL();
}
// ============================================================================
// GLBox::setRotation()                                                        
//  Set the rotation angle of the object to n degrees around the X,Y,Z axis.   
void GLBox::setRotation( int x, int y, int z )
{
  // rotate angles
  xRot = (GLfloat)(x % 360);
  yRot = (GLfloat)(y % 360);
  zRot = (GLfloat)(z % 360);
  // hud display
  hud->setText(GLHudObject::Rot,xRot,yRot,zRot);
  hud->updateDisplay();
  // save values
  store_options->xrot =xRot;
  store_options->yrot =yRot;
  store_options->zrot =zRot;
  updateGL();
}
// ============================================================================
// GLBox::setXRotation()                                                       
// Set the rotation angle of the object to n degrees around the X axis.        
void GLBox::setXRotation( int degrees )
{
  xRot = (GLfloat)(degrees % 360);
  store_options->xrot =xRot;
  updateGL();
}
// ============================================================================
// GLBox::setYRotation()                                                       
// Set the rotation angle of the object to n degrees around the Y axis.        
void GLBox::setYRotation( int degrees )
{
  yRot = (GLfloat)(degrees % 360);
  store_options->yrot =yRot;
  updateGL();
}
// ============================================================================
// GLBox::setZRotation()                                                       
// Set the rotation angle of the object to n degrees around the Z axis.        
void GLBox::setZRotation( int degrees )
{
  zRot = (GLfloat)(degrees % 360);
  store_options->zrot =zRot;
  updateGL();
}
// ============================================================================
// GLBox::getPixelTranslation()                                                
// return x,y,z pixel translation according to current translation viewport    
// coordinates and zoom                                                        
void GLBox::getPixelTranslation(int *x, int *y, int *z)
{
  GLint	Viewport[4];
  glGetIntegerv(GL_VIEWPORT,Viewport);
  *x = (int) (-store_options->xtrans * Viewport[2]/zoom);
  *y = (int) ( store_options->ytrans * Viewport[3]/zoom);
  *z = (int) ( store_options->ztrans * Viewport[2]/zoom);
}
// ============================================================================
// GLBox::setTranslation()                                                     
// Set the translation                                                         
void GLBox::setTranslation()
{
  setTranslation(ixTrans, iyTrans, izTrans);
}
// ============================================================================
// GLBox::setTranslation()                                                     
// Set the translation angle                                                   
void GLBox::setTranslation( int x, int y, int z )
{
  GLint	Viewport[4];
  glGetIntegerv(GL_VIEWPORT,Viewport);
  
  ixTrans = x;
  iyTrans = y;
  izTrans = z;
  // compute translation
  xTrans = (GLfloat)(-x*zoom/(Viewport[2]));
  yTrans = (GLfloat)(y*zoom/(Viewport[3])); //Viewport[3]*ratio));
  zTrans = (GLfloat)(z*zoom/(Viewport[2]));
  // diplsy on HUD
  hud->setText(GLHudObject::Trans,xTrans,yTrans,zTrans);
  hud->updateDisplay();
  // save
  store_options->xtrans=xTrans;
  store_options->ytrans=yTrans;
  store_options->ztrans=zTrans;
  updateGL();
}
// ============================================================================
// GLBox::toggleGrid()                                                         
// toggle grid                                                                 
void GLBox::toggleGrid( )
{
  store_options->show_grid = ! store_options->show_grid;
  updateGL();
}
// ============================================================================
// GLBox::resizeGrid()                                                         
// resize grid                                                                 
void GLBox::resizeGrid(float square_size, int nb_square )
{
  gridx->setNbSquare(nb_square);
  gridx->setSquareSize(square_size);
  gridx->rebuild();
  gridy->rebuild();
  gridz->rebuild();
  cube->setSquareSize(square_size * nb_square);
  cube->rebuild();
}
// ============================================================================
// GLBox::setParticlesSize()                                                   
// set particles size                                                          
void GLBox::setParticlesSize(int value, const bool ugl)
{
  particles_size = 1.0 + ((GLfloat) value) * MAX_PARTICLES_SIZE / 100.0;
  if (ugl) updateGL();
}
// ============================================================================
// GLBox::changeColorHUD()                                                     
// Set Text value to the specified HudObject                                   
void GLBox::changeColorHUD(const QColor color)
{
  hud->updateColor(color);
  updateGL();
}
// ============================================================================
// GLBox::setHud()                                                             
// Set Text value to the specified HudObject                                   
void GLBox::setHud(const GLHudObject::HudKeys k, const QString text)
{
  hud->setText(k,text);
  updateGL();
}
// ============================================================================
// GLBox::setHud()                                                             
// Set Float value to the specified HudObject                                  
void GLBox::setHud(const GLHudObject::HudKeys k, const float value)
{
  hud->setText(k,(const float) value);
  updateGL();
}
// ============================================================================
// GLBox::setHud()                                                             
// Set Int value to the specified HudObject                                    
void GLBox::setHud(const GLHudObject::HudKeys k, const int value)
{
  hud->setText(k,(const int) value);
  updateGL();
}
// ============================================================================
// GLBox::setHudToggle                                                         
// Toggle Hud object selected                                                  
void GLBox::setHudToggle(const GLHudObject::HudKeys k)
{
  hud->keysToggle(k);
  updateGL();
}
// ============================================================================
// GLBox::setHudActivate()                                                     
// Activate HudObject selected                                                 
void GLBox::setHudActivate(const GLHudObject::HudKeys k, const bool status)
{
  hud->keysActivate(k,status);
  updateGL();
}
// ============================================================================
// GLBox::setHudActivateNoGL()                                                 
// Activate all HUD object according to 'store_options->hud' status            
void GLBox::setHudActivateNoGL()
{
  if (store_options->hud) {  // HUD enable
    hud->setActivate(true);       
  } 
  else {                  // HUD disable
    hud->setActivate(false);
  }
  hud->keysActivate(GLHudObject::Getdata,store_options->hud_data_type);      
  hud->keysActivate(GLHudObject::Title,store_options->hud_title);      
  hud->keysActivate(GLHudObject::Time,store_options->hud_time);      
  hud->keysActivate(GLHudObject::Zoom,store_options->hud_zoom);      
  hud->keysActivate(GLHudObject::Rot,store_options->hud_rot);      
  hud->keysActivate(GLHudObject::Trans,store_options->hud_trans);      
  hud->keysActivate(GLHudObject::Nbody,store_options->hud_nbody);      
  hud->keysActivate(GLHudObject::Projection,store_options->hud_projection);  
  
  GlobalOptions * so = store_options;
  hud->setText(GLHudObject::Rot,so->xrot,so->yrot,so->zrot);
  hud->setText(GLHudObject::Trans,so->xtrans,so->ytrans,so->ztrans);
  hud->setText(GLHudObject::Zoom,so->zoom);
}
// ============================================================================
// GLBox::setHudActivate()                                                     
// Activate all HUD object according to 'store_options->hud' status            
void GLBox::setHudActivate()
{
  setHudActivateNoGL(); 
  updateGL();
}
// ============================================================================
// GLBox::setTextureSize()                                                     
// Change texture size for gaz like particles                                  
void GLBox::setTextureSize(const float ts, const bool ugl)
{
    for (int i=0; i<nb_object; i++) {
      vparticles_object[i]->setTextureSize(ts);
    }
    if (ugl) updateGL();
}
// ============================================================================
// GLBox::changeTextureAlphaColor()                                            
// Change Texture Alpha color for gaz like particles                           
void GLBox::changeTextureAlphaColor(const int alpha, const bool ugl)
{
    for (int i=0; i<nb_object; i++) {
      vparticles_object[i]->setTextureAlphaColor(alpha);
    }
    if (ugl) updateGL();
}
// ============================================================================
// GLBox::setProjection()                                                      
// set projection according to Width and Height window property                
float fx,fy;
void GLBox::setProjection(int width, int height)
{
  glViewport( 0, 0, width, height);
  ratio =  ((double )width) / ((double )height);
  glMatrixMode(GL_PROJECTION); 
  glLoadIdentity();
  
  // 
  if (store_options->perspective) {
    gluPerspective(45.,ratio,0.0005,(float) DOF);
  } 
  else {
    computeOrthoFactorRatio();
    glOrtho(ortho_left   * fx  * store_options->zoomo,
            ortho_right  * fx  * store_options->zoomo,
            ortho_bottom * fy  * store_options->zoomo,
            ortho_top    * fy  * store_options->zoomo,
            -1000,1000);
            //(float) -DOF/2.,(float) -DOF/2.);            
  }
  setViewPort();
}
// ============================================================================
// GLBox::setOrthoProjection()                                                 
// set Orthographic projection according to the range specified                
void GLBox::setOrthoProjection(float range)
{
  if (range == 0.0) {
    if (nb_object) { // we need at least one object
      float absxmax=fabs(vparticles_object[0]->coo_max[0]);
      float absymax=fabs(vparticles_object[0]->coo_max[1]);
      for (int i=0; i< nb_object; i++ ) {
        if (vparticles_object[i]->getActivate()) {
          float x=fabs(vparticles_object[i]->coo_max[0]);
          float y=fabs(vparticles_object[i]->coo_max[1]);
          if (x>absxmax) {
            absxmax = x;
          }
          if (y>absymax) {
            absymax = y;
          }
        }
      }
      //float range;
      if (absxmax >= absymax) {
        range = absxmax;
      }
      else {
        range = absymax;
      }
    }
    else {
      range=6.0; // default range if NO objecta
    }
  }
  //range=6.;
  ortho_right = range;
  ortho_left  =-range;
  ortho_top   = range;
  ortho_bottom=-range;
  computeOrthoFactorRatio();
  updateGL();
}
// ============================================================================
// GLBox::computeOrthoFactorRatio()                                            
// compute Orthographic factor ratio                                           
void GLBox::computeOrthoFactorRatio()
{
  if (ratio<1.0) {
    fx = 1.0;
    fy = 1./ratio;
  } 
  else {
    fx = ratio;
    fy = 1.0;
  }
}
// ============================================================================
// GLBox::loadImage()                                                          
// load embeded texture image                                                  
void GLBox::loadImage()
{
  //QString name = QFileDialog::getOpenFileName(".", "Images (*.png *.xpm *.jpg)", 
  //                                            this, "Choose", "Select an image");
  
  // In case of Cancel
  //if (name.isEmpty())
  //  return;
  //QString name="/home/jcl/download/PEngine/images/particle.png";
#if 0  
  QString name="/home/jcl/works/glnemo/src/images/smoke.png";
  QImage img(name);
#else  
  QString name="smoke";  
  QImage img=qembed_findImage(name);
#endif  
  if (img.isNull())
    {
      std::cerr << "Unable to load " << name << ", unsupported file format" << std::endl;
      return;
    }
    

  PRINT_D std::cout << "Loading " << name << ", " << img.width() 
                    << "x" << img.height() << " pixels" << std::endl;

  // 1E-3 needed. Just try with width=128 and see !
  int newWidth  = 1<<(int)(1+log(img.width() -1+1E-3) / log(2.0));
  int newHeight = 1<<(int)(1+log(img.height()-1+1E-3) / log(2.0));

  u_max = img.width()  / (float)newWidth;
  v_max = img.height() / (float)newHeight;
  
  if ((img.width()!=newWidth) || (img.height()!=newHeight))
    {
      PRINT_D std::cout << "Image size set to " << newWidth 
                        << "x" << newHeight << " pixels" << std::endl;
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
// ============================================================================
// GLBox::bestFit()                                                            
// find out the "best" zoom to fit all the particles in the drawing windows    
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
#if 0  
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
#endif  
  delete[] r;
  delete[] g;
  delete[] b;
}
// ============================================================================
// GLBox::updateOptions()                                                      
void GLBox::updateOptions(GlobalOptions * go, const bool only_transform)
{
  if (only_transform) {               // only                    
    store_options->copyTransform(*go); // copy transformation data
  } else {                            // else                    
    *store_options = *go;             // copy all                
  }
  setHudActivateNoGL(); 
  updateGL();
}
// ============================================================================
// GLBox::updateOptions()                                                      
void GLBox::takeScreenshot(QImage & img)
{
  img=grabFrameBuffer();
}
// ============================================================================
// GLBox::treeUpdate()                                                         
void GLBox::treeUpdate() 
{ 
  tree->update(); 
  for (int i=0; i<nb_object; i++) {
    vparticles_object[i]->rebuildDisplayList();
  }
 updateGL();
}
// ============================================================================
