// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2009                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
/**
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
*/
#ifndef GLNEMOGLWINDOW_H
#define GLNEMOGLWINDOW_H
#include <QGLWidget>
#include <QImage>
#include <QMutex>
#include "particlesobject.h"
#include "globjectparticles.h"
#include "glselection.h"
#include "globjectosd.h"
#include "gltexture.h"
#include "gloctree.h"
#include "camera.h"
namespace glnemo {

class GLGridObject;
class GlobalOptions;

class GLWindow : public QGLWidget {
  Q_OBJECT
public:
    GLWindow(QWidget *, GlobalOptions *, QMutex * , Camera * );
    ~GLWindow();

    void bestZoomFit();
    void resize(const int w, const int h ) { resizeGL(w,h);}
    void resizeOsd(const int w, const int h ) { osd->setWH(w,h);};
    void resetView() {    // reset view to initial
      setRotation(0,0,0);
      setTranslation(0,0,0);
      resetEvents();
    }
    static GLuint m_program;    
    static bool GLSL_support;
    void setFBO(bool _b) { fbo = _b; };
    void setFBOSize(GLuint w, GLuint h) { texWidth=w; texHeight=h;};
    QImage grabFrameBufferObject() { return imgFBO;};
    void rotateAroundAxis(const int); 
  
    // select area
    GLSelection * gl_select;
    static void checkGLErrors(std::string s);
    
  signals:
    void sigKeyMouse(const bool, const bool);
    void sigScreenshot();
public slots:
   void  update(ParticlesData   * ,
                ParticlesObjectVector * ,
                GlobalOptions         * );
   void  update(ParticlesObjectVector * );
   void  update();
   void  updateVbo(const int);
   void  updateColorVbo(const int);
   void  changeColorMap();
   void  reverseColorMap();
   void  updateGL();
   void  osdZoom(bool ugl=true);
   void setOsd(const GLObjectOsd::OsdKeys k,const QString text, bool b=true);
   void setOsd(const GLObjectOsd::OsdKeys k,const int value, bool b=true);
   void setOsd(const GLObjectOsd::OsdKeys k,const float value, bool b=true);
   void resetFrame() { nframe=0; };
   int getFrame() { return nframe;};

protected:
  void	initializeGL();
  void	paintGL();
  void	resizeGL( int w, int h );

private slots:
  void updateVel(const  int); // update velocity vector
  void rotateAroundX() { rotateAroundAxis(0);};
  void rotateAroundY() { rotateAroundAxis(1);};
  void rotateAroundZ() { rotateAroundAxis(2);};
  void translateX()    { translateAlongAxis(0); };
  void translateY()    { translateAlongAxis(1); };
  void translateZ()    { translateAlongAxis(2); };
  void setTextureObject(const int, const int);
private:
  // my parent
  QWidget * parent;
  // global options
  GlobalOptions * store_options;
  // grid variables
  GLGridObject * gridx, * gridy, * gridz;
  // Vectors
  GLObjectParticlesVector gpv;
  ParticlesObjectVector * pov;
  ParticlesData   * p_data;
  // projections
  void setProjection(const int w, const int h);
  void computeOrthoFactor();
  float ratio, fx,fy;
  int wwidth, wheight;
  GLuint texWidth, texHeight;
  QImage imgFBO; 
  bool fbo;
  // events
  void mousePressEvent  ( QMouseEvent *e );
  void mouseReleaseEvent( QMouseEvent *e );
  void mouseMoveEvent   ( QMouseEvent *e );
  void wheelEvent       ( QWheelEvent *e );
  void keyReleaseEvent  ( QKeyEvent   *k );
  void keyPressEvent    ( QKeyEvent   *k );
  void resetEvents();
  bool is_pressed_left_button;
  bool is_pressed_right_button;
  bool is_pressed_middle_button;
  bool is_pressed_ctrl;
  bool is_mouse_pressed;
  bool is_mouse_zoom;
  bool is_key_pressed;
  bool is_ctrl_pressed;

  void translateAlongAxis(const int);
  // transformations (rotation, translation)
  bool is_translation;
  int  x_mouse, y_mouse, z_mouse,
      tx_mouse,ty_mouse,tz_mouse,
      last_posx, last_posy, last_posz;
  void setRotation( const int x, const int y, const int z );
  void getPixelTranslation(int *x, int *y, int *z);
  void setTranslation( const int x, const int y, const int z );
  void setZoom(const int z);
  void setZoom(const float z);
  int zoom_dynam;
  // gl matrix
  GLdouble mProj[16], mModel[16], mModel2[16];
  int viewport[4];
  void setProjMatrix()  {
    glGetDoublev(GL_PROJECTION_MATRIX, (GLdouble *) mProj);
  };
  void setModelMatrix() {
    glGetDoublev(GL_MODELVIEW_MATRIX, (GLdouble *) mModel);
  };
  void setViewPort() {
    glGetIntegerv(GL_VIEWPORT,viewport);
  };
  // OSD
  GLObjectOsd * osd;
  QImage image,gldata;
  void formatInstructions(int width, int height);
  // Texture vector
  GLTextureVector gtv;
  // Thread
  QMutex * mutex_data;
  
  bool is_shift_pressed;
  // bench
  int nframe;
  // Shaders
  void initShader();
  unsigned int m_vertexShader, m_pixelShader;

  static const char vertexShader[];
  static const char pixelShader[];
  // Lights
  void initLight();
  // octree
  GLOctree * tree;
  // Camera
  Camera * camera;
};
} // namespace glnemo

#endif
