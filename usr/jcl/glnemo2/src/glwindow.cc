// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2008                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           P�le de l'Etoile, site de Ch�teau-Gombert                         
//           38, rue Fr�d�ric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================

#include <QtGui>
#include <GL/glew.h>
#include <QtOpenGL>
#include <QMutex>
#include <assert.h>
#include "glwindow.h"
#include "glgridobject.h"
#include "globaloptions.h"
#include "particlesdata.h"
#include "particlesobject.h"
#include "tools3d.h"

namespace glnemo {
#define DOF 4000000

  const char  GLWindow::vertexShader[] = {
        "// with ATI hardware, uniform variable MUST be used by output          \n"
        "// variables. That's why win_height is used by gl_FrontColor           \n"
        "uniform int win_height;                                                \n"
        "uniform float alpha;                                                   \n"
        "attribute float a_sprite_size;                                         \n"
        "void main()                                                            \n"
        "{                                                                      \n"
        "    vec4 vert = gl_Vertex;                                             \n"
        "    float pointSize =  win_height*a_sprite_size*gl_Point.size;  \n"
        "    vec3 pos_eye = vec3 (gl_ModelViewMatrix * vert);                   \n"
        "    gl_PointSize = max(0.00001, pointSize / (1.0 - pos_eye.z));        \n"
        "    gl_TexCoord[0] = gl_MultiTexCoord0;                                \n"
        "    gl_Position = ftransform();                                        \n"
        "    gl_FrontColor =  vec4(gl_Color.r+win_height-win_height,gl_Color.g,gl_Color.b,gl_Color.a*alpha); \n"
        "}                                                                      \n"
  };
  const char  GLWindow::pixelShader[] = {
        "uniform sampler2D splatTexture;                                        \n"

        "void main()                                                            \n"
        "{                                                                      \n"
        "    vec4 color = gl_Color * texture2D(splatTexture, gl_TexCoord[0].st);\n"
        "    gl_FragColor = color ;                                              \n"
        "}                                                                      \n"
  };
  unsigned int GLWindow::m_program = 0;
  bool GLWindow::GLSL_support = false;
GLuint framebuffer, renderbuffer;
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
  resetEvents();
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
  connect(gl_select, SIGNAL(updateGL()), this, SLOT(updateGL()));
  connect(gl_select, SIGNAL(updateZoom()), this, SLOT(osdZoom()));
  connect(this,SIGNAL(sigScreenshot()),parent,SLOT(startAutoScreenshot()));
  setFocus();
  wwidth=894;wheight=633;
  zoom_dynam = 0;
  // OSD
  //QFont f=QFont("Courier", 12, QFont::Light);
  QFont f;
  //f.setFamily("fixed");
  f.setRawMode(true);
  f.setPixelSize(10);
  f.setFixedPitch (true)  ;
  //f.setStyleHint(QFont::AnyStyle, QFont::PreferBitmap);
  f.setStyleHint(QFont::SansSerif, QFont::PreferAntialias);
  osd = new GLObjectOsd(wwidth,wheight,f,Qt::yellow);
  initializeGL();
  checkGLErrors("initializeGL");
  initShader();
  checkGLErrors("initShader");
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
  delete tree;
  if (GLWindow::GLSL_support) {
    glDeleteRenderbuffersEXT(1, &renderbuffer);
    glDeleteRenderbuffersEXT(1, &framebuffer);
  }
}
#define COPY 0
// ============================================================================
// update
void GLWindow::updateGL()
{
  if ( !store_options->duplicate_mem) {
    mutex_data->lock();
    if (store_options->new_frame) update();
    else QGLWidget::updateGL();
    mutex_data->unlock();
  }
  else QGLWidget::updateGL();
}
//QMutex mutex1;

// ============================================================================
// update
void GLWindow::update(ParticlesData   * _p_data,
                      ParticlesObjectVector * _pov,
                      GlobalOptions         * _go)
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

  for (unsigned int i=0; i<pov->size() ;i++) {
    if (i>=gpv.size()) {
      GLObjectParticles * gp = new GLObjectParticles(p_data,&((*pov)[i]),store_options,&gtv);
      //GLObjectParticles * gp = new GLObjectParticles(&p_data,pov[i],store_options);
      gpv.push_back(*gp);
      delete gp;
    } else {
      gpv[i].update(p_data,&((*pov)[i]),store_options);
      //gpv[i].update(&p_data,pov[i],store_options);
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
     gpv[i].buildVboColor();
   }
  updateGL();
}
// ============================================================================
// reverseColorMap                                                             
void GLWindow::reverseColorMap()
{
  store_options->reverse_cmap = !store_options->reverse_cmap;
  for (unsigned int i=0; i<pov->size() ;i++) {
    gpv[i].buildVboColor();
  }
  updateGL();
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
void GLWindow::paintGL()
{
  if (store_options->auto_gl_screenshot) {
    store_options->auto_gl_screenshot = false;
    emit sigScreenshot();
    store_options->auto_gl_screenshot = true;
  }
  if ( !store_options->duplicate_mem)
    mutex_data->lock();
  if (fbo && GLWindow::GLSL_support) {
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
  setProjection( wwidth, wheight);
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();
  //glEnable(GL_DEPTH_TEST);
 
  // set camera
  camera->setEye(0.0,  0.0,  -store_options->zoom);
  camera->moveTo();

  // rotate the scene
  glRotatef( store_options->xrot, 1.0, 0.0, 0.0 );
  glRotatef( store_options->yrot, 0.0, 1.0, 0.0 );
  glRotatef( store_options->zrot, 0.0, 0.0, 1.0 );

  // Grid Anti aliasing
#ifdef GL_MULTISAMPLE
  glEnable(GL_MULTISAMPLE);
#endif
  if (1) { //line_aliased) {
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    //glLineWidth (0.61);
    glLineWidth (1.0);
  } else {
    glDisable(GL_LINE_SMOOTH);
  }

  // grid display
  if (store_options->show_grid) {
    glEnable(GL_BLEND);
    gridx->display();
    gridy->display();
    gridz->display();
    //!cube->display();
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
  //glPointSize(store_options->psize);
  
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
    for (int i=0; i<(int)pov->size(); i++) {
/*      if (i==0) {
        store_options->render_mode = 2;
      }
      else {
        store_options->render_mode = 0;
      }*/
      gpv[i].display(mModel2,wheight);
    }
    //mutex_data->unlock();
  }
  // octree
  if (store_options->octree_display || 1) {
    tree->display();
  }
  
  
  // On Screen Display
  if (store_options->show_osd) osd->display(this);
  
  // display selected area
  gl_select->display(QGLWidget::width(),QGLWidget::height());

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
}
// ============================================================================
void GLWindow::initShader()
{
  if (store_options->init_glsl) {
    GLSL_support = true;
    std::cerr << "begining init shader\n";
    glewInit();
    if (GLEW_ARB_vertex_shader && GLEW_ARB_fragment_shader && GL_VERSION_2_0)
      qDebug() << "Ready for GLSL\n";
    else {
      qDebug() << "BE CAREFULL : No GLSL support\n";
      GLSL_support = false;
      //exit(1);
    }
    if (GLSL_support ) {
      m_vertexShader = glCreateShader(GL_VERTEX_SHADER);
      if (!m_vertexShader) {
        qDebug() << "Unable to create VERTEX SHADER.....\n";
        exit(1);
      }
      m_pixelShader = glCreateShader(GL_FRAGMENT_SHADER);
      if (!m_pixelShader) {
        qDebug() << "Unable to create PIXEL SHADER.....\n";
        exit(1);
      }
      const char* v = vertexShader;
      const char* p = pixelShader;
      glShaderSource(m_vertexShader, 1, &v, NULL);
      glShaderSource(m_pixelShader, 1, &p, NULL);
      
      GLint compile_status;
      glCompileShader(m_vertexShader);
      checkGLErrors("compile Vertex Shader");
      glGetShaderiv(m_vertexShader, GL_COMPILE_STATUS, &compile_status);
      if(compile_status != GL_TRUE) {
        qDebug() << "Unable to COMPILE VERTEX SHADER.....\n";
        exit(1);
      }
      
      glCompileShader(m_pixelShader);
      checkGLErrors("compile Pixel Shader");
      glGetShaderiv(m_pixelShader, GL_COMPILE_STATUS, &compile_status);
      if(compile_status != GL_TRUE) {
        qDebug() << "Unable to COMPILE PIXEL SHADER.....\n";
        exit(1);
      }
      
      m_program = glCreateProgram();
      
      glAttachShader(m_program, m_vertexShader);
      glAttachShader(m_program, m_pixelShader);
      
      // bind attribute
      //glBindAttribLocation(m_program, 100, "a_sprite_size");
      glLinkProgram(m_program);
      checkGLErrors("link Shader program");
      int  link_status;
      glGetProgramiv(m_program, GL_LINK_STATUS, &link_status);
      if(link_status != GL_TRUE) {
        qDebug() << "Unable to LINK Shader program.....\n";
        exit(1);
      }
      glDeleteShader(m_vertexShader);
      glDeleteShader(m_pixelShader);
      std::cerr << "ending init shader\n";
    }
  }
  else { // Initialisation of GLSL not requested
    std::cerr << "GLSL desactivated from user request, slow rendering ...\n";
    GLSL_support = false;
  }
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
  std::cerr << ">>>>>>>>> initializeGL()\n";
#if 1
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
  //width = height = 0;

  makeCurrent();   // activate OpenGL context, can build display list by now
  //initShader();
  
  // grid
  gridx = new GLGridObject(0,QColor(136,141,102),store_options->xy_grid);
  gridy = new GLGridObject(1,QColor(136,141,102),store_options->yz_grid);
  gridz = new GLGridObject(2,QColor(136,141,102),store_options->xz_grid);

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
  
  //formatInstructions(w,h);
}
// ============================================================================
//
void GLWindow::formatInstructions(int width, int height)
{
  QString text = tr("- Glnemo 2 Demonstration -");
  QFontMetrics metrics = QFontMetrics(font());
  int border = qMax(4, metrics.leading());

  QRect rect = metrics.boundingRect(0, 0, width - 2*border, int(height*0.125),
      Qt::AlignCenter | Qt::TextWordWrap, text);
  
  image = QImage(width, rect.height() + 2*border, QImage::Format_ARGB32_Premultiplied);
  image.fill(qRgba(0, 0, 0, 127));

  QPainter painter;
  painter.begin(&image);
  painter.setRenderHint(QPainter::TextAntialiasing);
  //painter.setRenderHint(QPainter::HighQualityAntialiasing);
  painter.setPen(Qt::green);
  painter.drawText((width - rect.width())/2, border,
                    rect.width(), rect.height(),
                    Qt::AlignCenter | Qt::TextWordWrap, text);
  painter.end();
  gldata = QGLWidget::convertToGLFormat(image);
}
// ============================================================================
// set up the projection according to the width and height of the windows
void GLWindow::setProjection(const int width, const int height)
{
  glViewport( 0, 0, width, height);
  ratio =  ((double )width) / ((double )height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if (store_options->perspective) {
    gluPerspective(45.,ratio,0.0005,(float) DOF);
  }
  else {
#if 0
    computeOrthoFactorRatio();
    glOrtho(ortho_left   * fx  * store_options->zoomo,
            ortho_right  * fx  * store_options->zoomo,
            ortho_bottom * fy  * store_options->zoomo,
            ortho_top    * fy  * store_options->zoomo,
            -1000,1000);
            //(float) -DOF/2.,(float) -DOF/2.);
#endif
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
void GLWindow::resetEvents()
{
   x_mouse= y_mouse= z_mouse=0;
  tx_mouse=ty_mouse=tz_mouse=0;
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
  setMouseTracking(FALSE);
  //!statusBar()->message("Ready");
  //!options_form->downloadOptions(store_options);
#if DRAWBOX
  draw_box->draw( glbox,part_data, &psv, store_options,"");
#endif
  if (is_shift_pressed) {
    if ( !store_options->duplicate_mem) mutex_data->lock();
    gl_select->zoomOnArea(pov->size(),mProj,mModel,viewport);
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
    }
    else {
      setRotation(y_mouse,x_mouse,z_mouse);
    }
  }
  //!options_form->downloadOptions(store_options);

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
//  Set the rotation angle of the object to n degrees around the X,Y,Z axis.
void GLWindow::setRotation( const int x, const int y, const int z )
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
// -----------------------------------------------------------------------------
// rotateAroundAxis()
void GLWindow::rotateAroundAxis(const int axis)
{
  if (!is_key_pressed              && // no interactive user request
      !is_mouse_pressed) {
      switch (axis) {
      case 0: store_options->xrot++;  // X
              y_mouse = (int) store_options->xrot;
              break;
      case 1: store_options->yrot++;  // Y
              x_mouse = (int) store_options->yrot;
              break;
      case 2: store_options->zrot++;  // Z
              z_mouse = (int) store_options->zrot;
              break;
      }
    updateGL();
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
  osd->setText(GLObjectOsd::Trans,xTrans,yTrans,zTrans);
  osd->updateDisplay();
  // save
  store_options->xtrans=xTrans;
  store_options->ytrans=yTrans;
  store_options->ztrans=zTrans;
  updateGL();
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
void GLWindow::setOsd(const GLObjectOsd::OsdKeys k, const QString text, bool ugl)
{
  osd->setText(k,text);
  if (ugl) updateGL();
}
// ============================================================================
// GLWindow::setOsd()                                                             
// Set Float value to the specified HudObject                                  
void GLWindow::setOsd(const GLObjectOsd::OsdKeys k, const float value, bool ugl)
{
  osd->setText(k,(const float) value);
  if (ugl) updateGL();
}
// ============================================================================
// GLWindow::setOsd()                                                             
// Set Int value to the specified HudObject                                    
void GLWindow::setOsd(const GLObjectOsd::OsdKeys k, const int value, bool ugl)
{
  osd->setText(k,(const int) value);
  if (ugl) updateGL();
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
// Best Zoom fit
// fit all the particles on the screen from perspective view
void GLWindow::bestZoomFit()
{
  if ( !store_options->duplicate_mem) mutex_data->lock();
  Tools3D::bestZoomFromObject(mProj,mModel,
                              viewport, pov, p_data, store_options);
  osdZoom();
  if ( !store_options->duplicate_mem) mutex_data->unlock();
}
} // namespace glnemo
