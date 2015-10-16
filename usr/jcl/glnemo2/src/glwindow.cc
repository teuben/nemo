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
#include <QtGlobal>
#if QT_VERSION >= QT_VERSION_CHECK(5, 0, 0)
#include <GL/glew.h>
#include <QtGui>
#else
#include <QtGui>
#include <GL/glew.h>
#endif
#include <QtOpenGL>
#include <QMutex>
#include <assert.h>
#include <limits>
#include <math.h>
#include "glwindow.h"
#include "glgridobject.h"
#include "glcubeobject.h"
#include "globaloptions.h"
#include "particlesdata.h"
#include "particlesobject.h"
#include "tools3d.h"
#include "fnt.h"
//#include "cshader.h"

namespace glnemo {
#define DOF 4000000
  
  bool GLWindow::GLSL_support = false;
  GLuint framebuffer, renderbuffer;
  GLdouble GLWindow::mIdentity[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  //float store_options->ortho_range;
  
// ============================================================================
// Constructor                                                                 
// BEWARE when parent constructor QGLWidget(QGLFormat(QGL::SampleBuffers),_parent)
// is called, we get antialiasing during screenshot capture but we can loose    
// performance during rendering. You have been warned !!!!!                     
GLWindow::GLWindow(QWidget * _parent, GlobalOptions*_go,QMutex * _mutex, Camera *_camera) //:QGLWidget(QGLFormat(QGL::SampleBuffers),_parent)
{
  // copy parameters
  parent        = _parent;
  store_options = _go;
  camera        = _camera;
  //setAttribute(Qt::WA_NoSystemBackground);
  // reset coordinates
  resetEvents(true);
  is_mouse_pressed   = FALSE;
  is_mouse_zoom      = FALSE;
  is_key_pressed     = FALSE;
  is_shift_pressed   = FALSE;
  is_pressed_ctrl    = FALSE;
  gpv.clear();
  pov = NULL;
  p_data = NULL;
  mutex_data = _mutex;
  nframe = 0; // used during bench
  p_data = new ParticlesData();
  pov    = new ParticlesObjectVector();
  gl_select = new GLSelection();
  // rotatation mode is screen by default
  store_options->rotate_screen = true;
  // reset rotation matrixes
  resetMatScreen();
  resetMatScene();
  reset_screen_rotation = false;
  reset_scene_rotation  = false;
  // initialyse rotation variables
  last_xrot = last_urot = 0.;
  last_yrot = last_vrot = 0.;
  last_zrot = last_wrot = 0.;
  // SCENE/object matrix index
  i_umat = 0;
  i_vmat = 1;
  i_wmat = 2;
  
  connect(gl_select, SIGNAL(updateGL()), this, SLOT(updateGL()));
  //connect(gl_select, SIGNAL(updateZoom()), this, SLOT(osdZoom()));
  connect(gl_select, SIGNAL(updateZoom()), this, SLOT(updateOsdZrt()));
  connect(this,SIGNAL(sigScreenshot()),parent,SLOT(startAutoScreenshot()));
  setFocus();
  wwidth=894;wheight=633;
  zoom_dynam = 0;
  // leave events : reset event when we leave opengl windows
  connect(this,SIGNAL(leaveEvent()),this,SLOT(resetEvents()));
  
  initializeGL();
  checkGLErrors("initializeGL");
  shader = NULL;
  vel_shader = NULL;
  initShader();
  checkGLErrors("initShader");
  ////////
  
  // grid
  GLGridObject::nsquare = store_options->nb_meshs;
  GLGridObject::square_size = store_options->mesh_length;
  gridx = new GLGridObject(0,store_options->col_x_grid,store_options->xy_grid);
  gridy = new GLGridObject(1,store_options->col_y_grid,store_options->yz_grid);
  gridz = new GLGridObject(2,store_options->col_z_grid,store_options->xz_grid);

  // axes
  axes = new GLAxesObject();
  
  // cube
  cube  = new GLCubeObject(store_options->mesh_length*store_options->nb_meshs,store_options->col_cube,store_options->show_cube);
  // load texture
  GLTexture::loadTextureVector(gtv);
  
  // build display list in case of screenshot
  if (store_options->show_part && pov ) {
    //std::cerr << "GLWindow::initializeGL() => build display list\n";
    for (int i=0; i<(int)pov->size(); i++) {
      // !!!! DEACTIVATE gpv[i].buildDisplayList();;
      gpv[i].buildVelDisplayList();;
      gpv[i].setTexture();
      //gpv[i].buildVboPos();
    }
  }
  
  // Osd
  fntRenderer text;
  font = new fntTexFont(store_options->osd_font_name.toStdString().c_str());
  text.setFont(font);
  text.setPointSize(store_options->osd_font_size );
  osd = new GLObjectOsd(wwidth,wheight,text,store_options->osd_color);
  // colorbar
  gl_colorbar = new GLColorbar(store_options,true);
  
  ////////
  // FBO
  // Set the width and height appropriately for you image
  fbo = false;
  //Set up a FBO with one renderbuffer attachment
  // init octree
  tree = new GLOctree(store_options);
  tree->setActivate(true);
  if (GLWindow::GLSL_support) {
    glGenFramebuffersEXT(1, &framebuffer);
    glGenRenderbuffersEXT(1, &renderbuffer);
  }
  checkGLErrors("GLWindow constructor");
}

// ============================================================================
// Destructor
GLWindow::~GLWindow()
{
  delete gridx;
  delete gridy;
  delete gridz;
  delete gl_select;
  delete cube;
  delete tree;
  delete axes;
  if (GLWindow::GLSL_support) {
    glDeleteRenderbuffersEXT(1, &renderbuffer);
    glDeleteRenderbuffersEXT(1, &framebuffer);
    if (shader) delete shader;
    if (vel_shader) delete vel_shader;
  }
  std::cerr << "Destructor GLWindow::~GLWindow()\n";
}
#define COPY 0
// ============================================================================
// update
void GLWindow::updateGL()
{
  if ( !store_options->duplicate_mem) {
    mutex_data->lock();
    if (store_options->new_frame) update();
    else QGLWidget::update();
    mutex_data->unlock();
  }
  else QGLWidget::update();
}
//QMutex mutex1;

// ============================================================================
// update
void GLWindow::update(ParticlesData   * _p_data,
                      ParticlesObjectVector * _pov,
                      GlobalOptions         * _go,
                      const bool update_old_obj)
{
  store_options = _go;
  mutex_data->lock();

  if (store_options->duplicate_mem) {
    //pov->clear();
    //*pov    = *_pov;
  }
  //else pov    = _pov;
  pov    = _pov;
  
  if (p_data!=_p_data) {
    if (store_options->duplicate_mem)  *p_data = *_p_data;
    else                                p_data = _p_data;
  }
  // octree
  store_options->octree_enable = true;
  store_options->octree_display = true;
  store_options->octree_level = 0;
  //tree->update(p_data, _pov);
  gl_colorbar->update(&gpv,p_data->getPhysData(),store_options,mutex_data);
  
  for (unsigned int i=0; i<pov->size() ;i++) {
    if (i>=gpv.size()) {
      GLObjectParticles * gp = new GLObjectParticles(p_data,&((*pov)[i]),
                                                     store_options,&gtv,shader,vel_shader);
      //GLObjectParticles * gp = new GLObjectParticles(&p_data,pov[i],store_options);
      gpv.push_back(*gp);
      delete gp;
    } else {      
      gpv[i].update(p_data,&((*pov)[i]),store_options, update_old_obj);
      //gpv[i].update(&p_data ,pov[i],store_options);
        
    }
  }

  store_options->new_frame=false;
  gl_select->update(&gpv,store_options,mutex_data);
  mutex_data->unlock();
  updateGL();
  
}
// ============================================================================
// update
void GLWindow::update(ParticlesObjectVector * _pov)
{
  if (store_options->duplicate_mem) {
    pov->clear();
    *pov    = *_pov;
  }
  else pov    = _pov;

}
// ============================================================================
// updateBondaryPhys
void GLWindow::updateBoundaryPhys(const int i_obj, const bool ugl)
{
  assert(i_obj < (int) gpv.size());
  gpv[i_obj].updateBoundaryPhys();
  if (ugl) updateGL();
}
// ============================================================================
// updateVbo
void GLWindow::updateVbo(const int i_obj)
{
  assert(i_obj < (int) gpv.size());
  gpv[i_obj].updateVbo();
}
// ============================================================================
// updateColorVbo
void GLWindow::updateColorVbo(const int i_obj)
{
  assert(i_obj < (int) gpv.size());
  gpv[i_obj].updateColorVbo();
}
// ============================================================================
// update
void GLWindow::update()
{
  if (pov) {
    update(p_data,pov,store_options);
  }
}
// ============================================================================
// update velocity vectors for the [index] object
void GLWindow::updateVel(const int index)
{
  if (pov && index < (int) pov->size()) {
    gpv[index].updateVel();
  }
}
// ============================================================================
// changeColorMap                                                              
void GLWindow::changeColorMap()
{
  for (unsigned int i=0; i<pov->size() ;i++) {
     //gpv[i].buildVboColor();
    gpv[i].updateColormap();
   }
  updateGL();
}
// ============================================================================
// reverseColorMap                                                             
void GLWindow::reverseColorMap()
{
  //store_options->reverse_cmap = !store_options->reverse_cmap;
//  for (unsigned int i=0; i<pov->size() ;i++) {
//    gpv[i].buildVboColor();
//  }
  updateGL();
}
// ============================================================================
// rebuildGrid                                                             
void GLWindow::rebuildGrid(bool ugl)
{
  GLGridObject::nsquare = store_options->nb_meshs;
  GLGridObject::square_size = store_options->mesh_length;
  gridx->rebuild();
  gridy->rebuild();
  gridz->rebuild();
  cube->setSquareSize(store_options->nb_meshs*store_options->mesh_length);
  if (ugl) updateGL();
}
// ============================================================================
// updatedGrid                                                             
void GLWindow::updateGrid(bool ugl)
{
  gridx->setActivate(store_options->xy_grid);
  gridx->setColor(store_options->col_x_grid);
  
  gridy->setActivate(store_options->yz_grid);
  gridy->setColor(store_options->col_y_grid);
  
  gridz->setActivate(store_options->xz_grid);
  gridz->setColor(store_options->col_z_grid);
  
  cube->setActivate(store_options->show_cube);
  cube->setColor(store_options->col_cube);
  
  if (ugl) updateGL();
}
// ============================================================================
// init Light                                                                  
void GLWindow::initLight()
{
}
#define OSD 0
// ============================================================================
// move, translate and re-draw the whole scene according to the objects and
// features selected
long int CPT=0;
void GLWindow::paintGL()
{
  CPT++; 
  //std::cerr << "GLWindow::paintGL() --> "<<CPT<<"\n";
  if (store_options->auto_gl_screenshot) {
    store_options->auto_gl_screenshot = false;
    emit sigScreenshot();
    store_options->auto_gl_screenshot = true;
  }
  if ( !store_options->duplicate_mem)
    mutex_data->lock();
  if (fbo && GLWindow::GLSL_support) {
    std::cerr << "FBO GLWindow::paintGL() --> "<<CPT<<"\n";
    //glGenFramebuffersEXT(1, &framebuffer);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, framebuffer);
    //glGenRenderbuffersEXT(1, &renderbuffer);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, renderbuffer);
    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_RGBA8, texWidth, texHeight);
    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                  GL_RENDERBUFFER_EXT, renderbuffer);
    GLuint status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    if (status != GL_FRAMEBUFFER_COMPLETE_EXT) {
    }
  } 
  //setFocus();
  
  qglClearColor(store_options->background_color);
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  // set projection
  setProjection(0, 0,  wwidth, wheight);
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();
  //glEnable(GL_DEPTH_TEST);
  
  // rotation around scene/object axes
  float ru=store_options->urot-last_urot;
  float rv=store_options->vrot-last_vrot;
  float rw=store_options->wrot-last_wrot;
  
  // the following code compute OpenGL rotation 
  // around UVW scene/object axes
  if (ru!=0 ||
      rv!=0 ||
      rw!=0) {
    glLoadIdentity();
    if (ru!=0)
      glRotatef(ru, mScene[0],mScene[1], mScene[2] );
    if (rv!=0)
      glRotatef(rv, mScene[4],mScene[5], mScene[6] );
    if (rw!=0)
      glRotatef(rw, mScene[8],mScene[9], mScene[10]);
    
    last_urot = store_options->urot;
    last_vrot = store_options->vrot;
    last_wrot = store_options->wrot;
    glMultMatrixd (mScene);
    glGetDoublev (GL_MODELVIEW_MATRIX, mScene);
  }
 
  // rotation around screen axes
  float rx=store_options->xrot-last_xrot;
  float ry=store_options->yrot-last_yrot;
  float rz=store_options->zrot-last_zrot;

  // the following code compute OpenGL rotation 
  // around XYZ screen axes
  if (rx!=0 ||
      ry!=0 ||
      rz!=0) {
    glLoadIdentity();
    // rotate only around the screen axes about the delta angle from the previous
    // rotation, otherwise it mess up the rotation
    glRotatef( rx, 1.0, 0.0, 0.0 );
    glRotatef( ry, 0.0, 1.0, 0.0 );
    glRotatef( rz, 0.0, 0.0, 1.0 );
    last_xrot = store_options->xrot;
    last_yrot = store_options->yrot;
    last_zrot = store_options->zrot;
    
    glMultMatrixd (mScreen); // apply previous rotations on the current one
    glGetDoublev (GL_MODELVIEW_MATRIX, mScreen); // save screen rotation matrix
  }
  if (reset_screen_rotation) { 
    glLoadIdentity ();
    glGetDoublev (GL_MODELVIEW_MATRIX, mScreen); // set to Identity
    reset_screen_rotation=false;
  }
  if (reset_scene_rotation) { 
    glLoadIdentity ();
    glGetDoublev (GL_MODELVIEW_MATRIX, mScene); // set to Identity
    reset_scene_rotation=false;
    last_urot = last_vrot = last_wrot = 0.0;
  }  

  glLoadIdentity (); // reset OGL rotations
  // set camera
  if ( store_options->perspective) {
    camera->setEye(0.0,  0.0,  -store_options->zoom);
    camera->moveTo();
  }
  glGetDoublev(GL_MODELVIEW_MATRIX, (GLdouble *) mRot);
  
  // apply screen rotation on the whole system
  glMultMatrixd (mScreen);   
  // apply scene/world rotation on the whole system
  glMultMatrixd (mScene);   
  
  // Grid Anti aliasing
#ifdef GL_MULTISAMPLE
  glEnable(GL_MULTISAMPLE);
#endif
  if (1) { //line_aliased) {
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);    
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    //glLineWidth (0.61);
    glLineWidth (1.0);
  } else {
    glDisable(GL_LINE_SMOOTH);
  }

  // grid display
  if (store_options->show_grid) {
    //glEnable( GL_DEPTH_TEST );
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    gridx->display();
    gridy->display();
    gridz->display();
    cube->display();
    glDisable(GL_BLEND);
  }

  // sphere display
  if (0) {
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    GLUquadricObj *quadric=gluNewQuadric();
    gluQuadricDrawStyle(quadric,GLU_LINE);
    gluQuadricNormals(quadric, GLU_SMOOTH);
    GLdouble radius=GLdouble(store_options->mesh_length*store_options->nb_meshs/2.0);
    GLint subdivisions=16;
    gluSphere(quadric, radius, subdivisions,subdivisions);
    gluDeleteQuadric(quadric);
    glDisable(GL_BLEND);

  }
  // camera display path and control points
  camera->display();

  setModelMatrix(); // save ModelView  Matrix
  setProjMatrix();  // save Projection Matrix
  // move the scene
  glTranslatef( store_options->xtrans, store_options->ytrans, store_options->ztrans);
  glGetDoublev(GL_MODELVIEW_MATRIX, (GLdouble *) mModel2);  
  
  // nice points display
  glEnable(GL_POINT_SMOOTH);
  
  // control blending on particles
  if (store_options->blending) {
    glEnable(GL_BLEND);
    glBlendFunc( GL_SRC_ALPHA, GL_ONE ); // original
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //glDepthFunc(GL_LESS);
  }
  else
    glDisable(GL_BLEND);
  // control depht buffer on particles
  if (store_options->dbuffer) glEnable (GL_DEPTH_TEST);
  else                        glDisable(GL_DEPTH_TEST);
  //glDepthFunc(GL_LESS);
  // Display objects (particles and velocity vectors)
  if (store_options->show_part && pov ) {
    //mutex_data->lock();
    bool first=true;
    bool obj_has_physic=false;
    for (int i=0; i<(int)pov->size(); i++) {
      gpv[i].display(mModel2,wheight);

      if (first) {
        const ParticlesObject * po = gpv[i].getPartObj();
        if (po->hasPhysic()) { //store_options->phys_min_glob!=-1 && store_options->phys_max_glob!=-1) {
          obj_has_physic=true;
          first=false;
        }
      }
    }
    if (obj_has_physic) {
      if (fbo) // offscreen rendering activated
        gl_colorbar->display(texWidth,texHeight);
      else
        gl_colorbar->display(QGLWidget::width(),QGLWidget::height());
    }

    //mutex_data->unlock();
  }
  
  // octree
  if (store_options->octree_display || 1) {
    tree->display();
  }
  
  // On Screen Display
  if (store_options->show_osd) osd->display();
    
  // display selected area
  gl_select->display(QGLWidget::width(),QGLWidget::height());

  // draw axes
  if (store_options->axes_enable)
    axes->display(mScreen, mScene, wwidth,wheight,
                  store_options->axes_loc,store_options->axes_psize, store_options->perspective);

  // reset viewport to the windows size because axes object modidy it
  glViewport(0, 0,  wwidth, wheight);

  if (fbo && GLWindow::GLSL_support) {
    fbo = false;
    //imgFBO = grabFrameBuffer();
    imgFBO = QImage( texWidth, texHeight,QImage::Format_RGB32);
    glReadPixels( 0, 0, texWidth, texHeight, GL_RGBA, GL_UNSIGNED_BYTE, imgFBO.bits() );
    // Make the window the target
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

   // Delete the renderbuffer attachment
   //glDeleteRenderbuffersEXT(1, &renderbuffer);
   //glDeleteRenderbuffersEXT(1, &framebuffer);
  } 
  if ( !store_options->duplicate_mem) mutex_data->unlock();

  nframe++; // count frames
  //glDrawPixels(gldata.width(), gldata.height(), GL_RGBA, GL_UNSIGNED_BYTE, gldata.bits());
  emit doneRendering();
}
// ============================================================================
void GLWindow::initShader()
{
  if (store_options->init_glsl) {
    const GLubyte* gl_version=glGetString ( GL_VERSION );
    std::cerr << "OpenGL version : ["<< gl_version << "]\n";
    int major = 0;
    int minor = 0;
    glGetIntegerv(GL_MAJOR_VERSION, &major);
    glGetIntegerv(GL_MINOR_VERSION, &minor);
    std::cerr << "OpenGL :"<< major << "." << minor << "\n";
    GLSL_support = true;
    std::cerr << "begining init shader\n";
    int err=glewInit();
    if (err==GLEW_OK && GLEW_ARB_multitexture && GLEW_ARB_vertex_shader && GLEW_ARB_fragment_shader && GL_VERSION_2_0)
      qDebug() << "Ready for GLSL\n";
    else {
      qDebug() << "BE CAREFULL : No GLSL support\n";
      GLSL_support = false;
      //exit(1);
    }

    if (GLSL_support ) {
      // check GLSL version supported
      const GLubyte* glsl_version=glGetString ( GL_SHADING_LANGUAGE_VERSION );
      std::cerr << "GLSL version supported : ["<< glsl_version << "]\n";
      //GLuint glsl_num;
      //glGetStringi(GL_SHADING_LANGUAGE_VERSION,glsl_num);
      //std::cerr << "GLSL version NUM : ["<< glsl_num << "]\n";
      // particles shader
      shader = new CShader(GlobalOptions::RESPATH.toStdString()+"/shaders/particles.vert.cc",
                           GlobalOptions::RESPATH.toStdString()+"/shaders/particles.frag.cc");
      shader->init();
      // velocity shader
      if (1) {

#if 0
          // Geometry shader OpenGL 3.30 and above only
          vel_shader = new CShader(GlobalOptions::RESPATH.toStdString()+"/shaders/velocity.vert330.cc",
                                   GlobalOptions::RESPATH.toStdString()+"/shaders/velocity.frag330.cc",
                                   GlobalOptions::RESPATH.toStdString()+"/shaders/velocity.geom330.cc");

#else
          vel_shader = new CShader(GlobalOptions::RESPATH.toStdString()+"/shaders/velocity.vert.cc",
                                   GlobalOptions::RESPATH.toStdString()+"/shaders/velocity.frag.cc");

#endif
          if (!vel_shader->init() ) {
              delete vel_shader;
              vel_shader=NULL;
          }
      }
    }

  }
  else { // Initialisation of GLSL not requested
    std::cerr << "GLSL desactivated from user request, slow rendering ...\n";
    GLSL_support = false;
  }
  std::cerr << "END OF INITSHADER \n";
}
// ============================================================================
// check OpenGL error message                                                  
void GLWindow::checkGLErrors(std::string s) 
{
  GLenum error;
  while ((error = glGetError()) != GL_NO_ERROR) {
    std::cerr << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n";
    std::cerr << s << ": error - " << (char *) gluErrorString(error)<<"\n";
  }
}
// ============================================================================
// initialyse the openGl engine
void GLWindow::initializeGL()
{
  std::cerr << "\n>>>>>>>>> initializeGL()\n\n";
#if 0
  qglClearColor( Qt::black );		// Let OpenGL clear to black
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LINE_SMOOTH);
#ifdef GL_MULTISAMPLE
  glEnable(GL_MULTISAMPLE);
#endif
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  //store_options->zoom = -10.;
  // Nice texture coordinate interpolation
  glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );
  //PRINT_D std::cerr << "-- Initialize GL --\n";

#endif

  makeCurrent();   // activate OpenGL context, can build display list by now
  
}
// ============================================================================
// resize the opengl viewport according to the new window size
void GLWindow::resizeGL(int w, int h)
{
  wwidth = w;
  wheight= h;
  glViewport( 0, 0, (GLint)w, (GLint)h );
  //float ratio =  ((double )w) / ((double )h);
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glFrustum(-1.0, 1.0, -1.0, 1.0, 1.0, (float) DOF);
  //gluPerspective(45.,ratio,0.0005,(float) DOF);
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();
  osd->setWH(w,h);
}
// ============================================================================
// set up the projection according to the width and height of the windows
void GLWindow::setProjection(const int x, const int y, const int width, const int height)
{
  glViewport( x, y, width, height);
  ratio =  ((double )width) / ((double )height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if (store_options->perspective) {
    gluPerspective(45.,ratio,0.0005,(float) DOF);
//    double mp[16];
//    glGetDoublev(GL_PROJECTION_MATRIX, (GLdouble *) mp);
//    for (int i=0;i<16;i++) std::cerr << "// "<< mp[i];
//    std::cerr << "\n";
  }
  else {
    computeOrthoFactor();    
    //std::cerr << "RANGE="<<store_options->ortho_range<<" zoom="<<store_options->zoom<<" zoomo="<<store_options->zoomo<<"\n";
    //std::cerr << "fx = "<< fx << " fy=" << fy <<  " range*fx*zoom0=" <<  store_options->zoomo*fx*store_options->ortho_range <<  "\n";
    ortho_right = store_options->ortho_range;
    ortho_left  =-store_options->ortho_range;
    ortho_top   = store_options->ortho_range;
    ortho_bottom=-store_options->ortho_range;
    //std::cerr << "zoom0 = " << store_options->zoomo << "\n";
    glOrtho(ortho_left   * fx  * store_options->zoomo,
            ortho_right  * fx  * store_options->zoomo,
            ortho_bottom * fy  * store_options->zoomo,
            ortho_top    * fy  * store_options->zoomo,
            -100000,100000);
            //(float) -DOF/2.,(float) -DOF/2.);
  }
  glGetIntegerv(GL_VIEWPORT,viewport);
}
// ============================================================================
// compute some factors for the orthographic projection
void GLWindow::computeOrthoFactor()
{
  if (ratio<1.0) {
    fx = 1.0  ; fy = 1./ratio;
  }
  else {
    fx = ratio; fy = 1.0;
  }
}
// ============================================================================
// reset rotation and translation to 0,0,0 coordinates
void GLWindow::resetEvents(bool pos)
{
  if (pos) {
    x_mouse= y_mouse= z_mouse=0;
    tx_mouse=ty_mouse=tz_mouse=0;
  }
  is_pressed_left_button  =FALSE;
  is_pressed_right_button =FALSE;
  is_pressed_middle_button=FALSE;
  is_mouse_pressed        =FALSE;
  is_translation          =FALSE;
  is_ctrl_pressed         =FALSE;
}
// ============================================================================
// manage rotation/translation according to mousePresssEvent
void GLWindow::mousePressEvent( QMouseEvent *e )
{
  setFocus();
  if ( e->button() == Qt::LeftButton ) {  // left button pressed
    if (is_shift_pressed)
      setCursor(Qt::CrossCursor);
    is_mouse_pressed       = TRUE;
    is_pressed_left_button = TRUE;
    setMouseTracking(TRUE);
    last_posx = e->x();
    last_posy = e->y();
    if  (is_translation) {;} //!parent->statusBar()->message("Translating X/Y");
    else                 {;} //!parent->statusBar()->message("Rotating X/Y");
  }
  if ( e->button() == Qt::RightButton ) { // right button pressed
    is_mouse_pressed        = TRUE;
    is_pressed_right_button = TRUE;
    setMouseTracking(TRUE);
    last_posz = e->x();
    if (is_translation) {;} //!parent->statusBar()->message("Translating Z");
    else                {;} //!parent->statusBar()->message("Rotating Z");
  }
  //if ( e->button() == Qt::MiddleButton ) {
  if ( e->button() == Qt::MidButton ) {
    //std::cerr << "Middle button pressed\n";
    is_mouse_pressed        = TRUE;
    is_pressed_middle_button= TRUE;
    setMouseTracking(TRUE);
    last_posx = e->x();
    last_posy = e->y();
  }
  emit sigKeyMouse( is_key_pressed, is_mouse_pressed);
  //!options_form->downloadOptions(store_options);
}
// ============================================================================
// GLObjectWindow::mouseReleaseEvent()
// manage mouseReleaseEvent
void GLWindow::mouseReleaseEvent( QMouseEvent *e )
{
  if (e) {;}  // do nothing... just to remove the warning :p
  is_pressed_left_button = FALSE;
  is_pressed_right_button = FALSE;
  is_mouse_pressed = FALSE;
  is_pressed_middle_button = FALSE;
  setMouseTracking(FALSE);
  //!statusBar()->message("Ready");
  //!options_form->downloadOptions(store_options);
#if DRAWBOX
  draw_box->draw( glbox,part_data, &psv, store_options,"");
#endif
  if (is_shift_pressed) {
    if ( !store_options->duplicate_mem) mutex_data->lock();
    //JCL 07/21/2015 setPerspectiveMatrix(); // toggle to perspective matrix mode
    gl_select->selectOnArea(pov->size(),mProj,mModel,viewport);
    setPerspectiveMatrix(); // toggle to perspective matrix mode
    gl_select->zoomOnArea(mProj,mModel,viewport);
    osd->setText(GLObjectOsd::Zoom,(const float) store_options->zoom);
    osd->updateDisplay();
    if ( !store_options->duplicate_mem) mutex_data->unlock();
  }
  setCursor(Qt::ArrowCursor);
  emit sigKeyMouse( is_key_pressed, is_mouse_pressed);
  //!draw_box->show();
}

// ============================================================================
// GLObjectWindow::mouseMoveEvent()
// manage mouseMoveEvent
void GLWindow::mouseMoveEvent( QMouseEvent *e )
{
  int dx=0,dy=0,dz=0;
  setFocus();
  if (is_pressed_left_button) {
    // offset displcacement
    dx = e->x()-last_posx;
    dy = e->y()-last_posy;
    //std::cerr << "dxdy="<< dx << " " << dy << "\n";
    // save last position
    last_posx = e->x();
    last_posy = e->y();
    if (is_shift_pressed && !is_mouse_zoom) { // user selection request
      gl_select->getMouse(e);
      updateGL();
    }
    else
      if (is_translation) {
        // total rotation
        tx_mouse+=dx;
        ty_mouse+=dy;
      }
      else {
        if (is_mouse_zoom) {
          setZoom(-dx);
        }
        else {
          // total rotation
          x_mouse+=dx;
          y_mouse+=dy;
        }
      }
  }
  if ( !gl_select->isEnable()) {
    if ( is_pressed_right_button) {
      // offset displcacement
      dz = e->x()-last_posz;
      // save last position
      last_posz = e->x();
      if (is_translation) {
        tz_mouse-=dz; // total rotation
      }
      else {
        z_mouse-=dz;  // total translation
      }
    }
    if (is_translation) {
      setTranslation(tx_mouse,ty_mouse,tz_mouse);
      //setTranslation(dx,dy,-dz);
    }
    else {
      //std::cerr << "xyz mouse="<<x_mouse<<" "<<y_mouse<<" "<< z_mouse<<"\n";
      if (store_options->rotate_screen)
        setRotationScreen(y_mouse,x_mouse,z_mouse);
      else
        setRotationScene(y_mouse,x_mouse,z_mouse);      
    }
  }
  //!options_form->downloadOptions(store_options);
  if (is_pressed_middle_button) {
    dx = e->x()-last_posx;
    dy = e->y()-last_posy;
    // save last position
    last_posx = e->x();
    last_posy = e->y();
    emit sigMouseXY(dx,dy);
    //std::cerr << "dx="<<dx<< "  dy="<<dy<<"\n";
  }
}
// ============================================================================
// manage zoom according to wheel event
void GLWindow::wheelEvent(QWheelEvent * e)
{
  setZoom(e->delta());
  //!options_form->downloadOptions(store_options);
}
// ============================================================================
// manage keyboard press events
void GLWindow::keyPressEvent(QKeyEvent * k)
{
  setFocus();
  if (k->key() == Qt::Key_Control ) {
    is_key_pressed = TRUE;
    is_translation = TRUE;
    is_pressed_ctrl= TRUE;
    getPixelTranslation(&tx_mouse,&ty_mouse,&tz_mouse);
    if (is_shift_pressed) {
      is_translation=FALSE;
      is_mouse_zoom=TRUE;
    }
  }
  if (k->key() == Qt::Key_A) {
    is_key_pressed = TRUE;
    //!glbox->toggleLineAliased();
  }
  if (k->key() == Qt::Key_Plus) {
    is_key_pressed = TRUE;
    setZoom(-1);
    //!statusBar()->message("Zoom IN");
  }
  if (k->key() == Qt::Key_Minus) {
    is_key_pressed = TRUE;
    setZoom(1);
    //!statusBar()->message("Zoom OUT");
  }
  if (k->key() == Qt::Key_Shift) {
    is_shift_pressed = TRUE;
    if (is_pressed_ctrl) {
      is_translation=FALSE;
      is_mouse_zoom=TRUE;
    }
  }
  emit sigKeyMouse( is_key_pressed, is_mouse_pressed);
  //!options_form->downloadOptions(store_options);
}
// ============================================================================
// manage keyboard release events
void GLWindow::keyReleaseEvent(QKeyEvent * k)
{
  if (k->key() == Qt::Key_Control ) {
    is_translation = FALSE;
    is_pressed_ctrl= FALSE;
    is_mouse_zoom  = FALSE;
  }
  if (k->key() == Qt::Key_Shift) {
    is_shift_pressed = FALSE;
    is_mouse_zoom  = FALSE;
    gl_select->reset();
    updateGL();
  }
  is_key_pressed = FALSE;
  emit sigKeyMouse( is_key_pressed, is_mouse_pressed);
  //!options_form->downloadOptions(store_options);
}
// ============================================================================
//  Set the rotation angle of the object to n degrees around the X,Y,Z axis of the screen
void GLWindow::setRotationScreen( const int x, const int y, const int z )
{
  // rotate angles
  GLfloat xRot = (GLfloat)(x % 360);
  GLfloat yRot = (GLfloat)(y % 360);
  GLfloat zRot = (GLfloat)(z % 360);
  // hud display
  osd->setText(GLObjectOsd::Rot,xRot,yRot,zRot);
  osd->updateDisplay();
  // save values
  store_options->xrot =xRot;
  store_options->yrot =yRot;
  store_options->zrot =zRot;
  updateGL();
}
// ============================================================================
//  Set the rotation angle of the object to n degrees around the U,V,W axis of the scene/object
void GLWindow::setRotationScene( const int x, const int y, const int z )
{
  // rotate angles
  GLfloat xRot = (GLfloat)(x % 360);
  GLfloat yRot = (GLfloat)(y % 360);
  GLfloat zRot = (GLfloat)(z % 360);
//  // hud display
//  osd->setText(GLObjectOsd::Rot,xRot,yRot,zRot);
//  osd->updateDisplay();
  // save values
  store_options->urot =xRot;
  store_options->vrot =yRot;
  store_options->wrot =zRot;
  updateGL();
}
// -----------------------------------------------------------------------------
// rotateAroundAxis()
void GLWindow::rotateAroundAxis(const int axis)
{
  if (!is_key_pressed              && // no interactive user request
      !is_mouse_pressed) {

      switch (axis) {
      case 0: // X
              store_options->xrot += store_options->ixrot;
              y_mouse = (int) store_options->xrot;              
              break;
      case 1: // Y
              store_options->yrot += store_options->iyrot;
              x_mouse = (int) store_options->yrot;
              break;
      case 2: // Z
              store_options->zrot += store_options->izrot;
              z_mouse = (int) store_options->zrot;
              break;
      case 3: // U
              store_options->urot += store_options->iurot;
              //y_mouse = (int) store_options->urot;              
              i_umat = 0;
              i_vmat = 1;
              i_wmat = 2;
              break;
      case 4: // V
              store_options->vrot += store_options->ivrot;
              //x_mouse = (int) store_options->vrot;
              i_umat = 4;
              i_vmat = 5;
              i_wmat = 6;              
              break;
      case 5: // W
              store_options->wrot += store_options->iwrot;
              //z_mouse = (int) store_options->wrot;
              i_umat = 8;
              i_vmat = 9;
              i_wmat = 10;              
              break;
      }
      if (axis < 3) {
        setRotationScreen(store_options->xrot,store_options->yrot,store_options->zrot);
      } else {
        setRotationScene(store_options->urot,store_options->vrot,store_options->wrot);
      }
      //setRotation(x,y,z);
      //updateGL();
  }
}
// -----------------------------------------------------------------------------
// translateAlongAxis()
void GLWindow::translateAlongAxis(const int axis)
{
  if (!is_key_pressed              && // no interactive user request
      !is_mouse_pressed) {
      switch (axis) {
      case 0: tx_mouse++;  // X
              setTranslation(tx_mouse,ty_mouse,tz_mouse);
              break;
      case 1: ty_mouse++;  // Y
              setTranslation(tx_mouse,ty_mouse,tz_mouse);
              break;
      case 2: tz_mouse--;  // Z
              setTranslation(tx_mouse,ty_mouse,tz_mouse);
              break;
      }
    updateGL();
  }
}
// ============================================================================
// GLBox::getPixelTranslation()
// return x,y,z pixel translation according to current translation viewport
// coordinates and zoom
void GLWindow::getPixelTranslation(int *x, int *y, int *z)
{
  GLint	Viewport[4];
  glGetIntegerv(GL_VIEWPORT,Viewport);
  *x = (int) (-store_options->xtrans * Viewport[2]/store_options->zoom);
  *y = (int) ( store_options->ytrans * Viewport[3]/store_options->zoom);
  *z = (int) ( store_options->ztrans * Viewport[2]/store_options->zoom);
}
// ============================================================================
// GLBox::setTranslation()
// Set the translation angle
void GLWindow::setTranslation( const int x, const int y, const int z )
{
  GLint	Viewport[4];
  glGetIntegerv(GL_VIEWPORT,Viewport);
#if 0
  ixTrans = x;
  iyTrans = y;
  izTrans = z;
#endif
  // compute translation
  GLfloat xTrans = (GLfloat)(-x*store_options->zoom/(Viewport[2]));
  GLfloat yTrans = (GLfloat)( y*store_options->zoom/(Viewport[3])); //Viewport[3]*ratio));
  GLfloat zTrans = (GLfloat)( z*store_options->zoom/(Viewport[2]));
  // display on HUD
  osd->setText(GLObjectOsd::Trans,-xTrans,-yTrans,-zTrans);
  osd->updateDisplay();
  // save
  store_options->xtrans=xTrans;
  store_options->ytrans=yTrans;
  store_options->ztrans=zTrans;
  updateGL();
}
// -----------------------------------------------------------------------------
// updateOsd()
void GLWindow::updateOsdZrt(bool ugl)
{
  GlobalOptions * g = store_options;
  setOsd(GLObjectOsd::Zoom,(const float) store_options->zoom,
                    g->osd_zoom,false);
  setOsd(GLObjectOsd::Rot,(const float) store_options->xrot,
                    (const float) store_options->yrot,
                    (const float) store_options->zrot, g->osd_rot,false);
  setOsd(GLObjectOsd::Trans,(const float) -store_options->xtrans,
                    (const float) -store_options->ytrans,
                    (const float) -store_options->ztrans,g->osd_trans,false);

  osd->updateDisplay();
  if (ugl) {
   updateGL();
  }

}

// ============================================================================
// setup zoom according to a z value
void GLWindow::setZoom(const int z)
{
  if (z>0) {
    store_options->zoom  *= 1.1;//1.1025; //1.1;
    store_options->zoomo *= 1.1;//1.1025; //1.1;
  } else {
    store_options->zoom  *= 0.9;//0.8075; //0.9;
    store_options->zoomo *= 0.9;//0.8075; //0.9;
  }
  osdZoom();
}
// ============================================================================
// update OSD zoom value
void GLWindow::osdZoom(bool ugl)
{
  osd->setText(GLObjectOsd::Zoom,(const float) store_options->zoom);
  osd->updateDisplay();
  if (ugl) updateGL();
}
// ============================================================================
// GLWindow::setOsd()                                                             
// Set Text value to the specified HudObject                                   
void GLWindow::setOsd(const GLObjectOsd::OsdKeys k, const QString text, bool show, bool ugl)
{
  osd->keysActivate(k,show);
  osd->setText(k,text);
  if (ugl) updateGL();
}
// ============================================================================
// GLWindow::setOsd()                                                             
// Set Float value to the specified HudObject                                  
void GLWindow::setOsd(const GLObjectOsd::OsdKeys k, const float value, bool show, bool ugl)
{
  osd->keysActivate(k,show);
  osd->setText(k,(const float) value);
  if (ugl) updateGL();
}
// ============================================================================
// GLWindow::setOsd()                                                             
// Set Float value to the specified HudObject                                  
void GLWindow::setOsd(const GLObjectOsd::OsdKeys k, const float value1, 
                      const float value2, const float value3, bool show,bool ugl)
{
  osd->keysActivate(k,show);
  osd->setText(k,(const float) value1,(const float) value2,(const float) value3);  
  if (ugl) updateGL();
}
// ============================================================================
// GLWindow::setOsd()                                                             
// Set Int value to the specified HudObject                                    
void GLWindow::setOsd(const GLObjectOsd::OsdKeys k, const int value, bool show, bool ugl)
{
  osd->keysActivate(k,show);
  osd->setText(k,(const int) value);
  if (ugl) updateGL();
}
// ============================================================================
// GLWindow::changeOsdFont()                                                             
// Change OSD font
void GLWindow::changeOsdFont()
{
  fntRenderer text;
  if (font) delete font;
  font = new fntTexFont(store_options->osd_font_name.toStdString().c_str());
  text.setFont(font);
  text.setPointSize(store_options->osd_font_size );
  osd->setFont(text);
  osd->setColor(store_options->osd_color);
  updateGL();
}
// ============================================================================
// set texture on the object
void GLWindow::setTextureObject(const int tex, const int obj)
{
  if (obj < (int) gpv.size() ) {
    gpv[obj].setTexture(tex);
  }
}
// ============================================================================
// setZoom()
// setup zoom according to a z value
void GLWindow::setZoom(const float z)
{
  const float zoom = z;
  store_options->zoom=zoom;
  osdZoom();
}
// ============================================================================
// >> HERE WE FORCE PERSPECTIVE PROJECTION
//    TO COMPUTE BOTH BEST ZOOM FOR ORTHO
//    AND PERSPECTVE PROJECTION
//
//    we must force perspective projection to have
//    te good prjection and modelview matrix
void GLWindow::setPerspectiveMatrix()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.,ratio,0.0005,(float) DOF);
#if 1
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();
  camera->setEye(0.0,  0.0,  -store_options->zoom);
  camera->moveTo();
  // apply screen rotation on the whole system
  glMultMatrixd (mScreen);   
  // apply scene/world rotation on the whole system
  glMultMatrixd (mScene);   
  setModelMatrix(); // save ModelView  Matrix
#endif
  setProjMatrix();  // save Projection Matrix
}

// ============================================================================
// Best Zoom fit
// fit all the particles on the screen from perspective view
void GLWindow::bestZoomFit()
{
  if ( !store_options->duplicate_mem) mutex_data->lock();

  glGetIntegerv(GL_VIEWPORT,viewport);
  ratio =  ((double )viewport[2]) / ((double )viewport[3]);

  setPerspectiveMatrix(); // toggle to perspective matric mode


  Tools3D::bestZoomFromObject(mProj,mModel,
                              viewport, pov, p_data, store_options);
    
  ortho_right = store_options->ortho_range;
  ortho_left  =-store_options->ortho_range;
  ortho_top   = store_options->ortho_range;
  ortho_bottom=-store_options->ortho_range;
  //store_options->zoomo = 1.;
  
  osdZoom();
  if ( !store_options->duplicate_mem) mutex_data->unlock();
}
} // namespace glnemo
