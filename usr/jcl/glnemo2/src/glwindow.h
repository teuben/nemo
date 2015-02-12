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
/**
	@author Jean-Charles Lambert <Jean-Charles.Lambert@lam.fr>
*/
#ifndef GLNEMOGLWINDOW_H
#define GLNEMOGLWINDOW_H

#include  "cshader.h"
#include <QGLWidget>
#include <QImage>
#include <QMutex>
#include "particlesobject.h"
#include "globjectparticles.h"
#include "glcubeobject.h"
#include "glselection.h"
#include "globjectosd.h"
#include "gltexture.h"
#include "gloctree.h"
#include "glcolorbar.h"
#include "glaxesobject.h"
#include "camera.h"


class fntTexFont;

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
    void resizeOsd(const int w, const int h ) { osd->setWH(w,h);}
    void resetView() {    // reset view to initial
      resetMatScreen();
      resetMatScene();
      reset_screen_rotation = true;
      reset_scene_rotation  = true;
      setRotationScreen(0,0,0);
      setRotationScene(0,0,0);
      setTranslation(0,0,0);
      resetEvents(true);
    }
    static bool GLSL_support;
    void setFBO(bool _b) { fbo = _b; }
    void setFBOSize(GLuint w, GLuint h) { texWidth=w; texHeight=h;}
    QImage grabFrameBufferObject() { return imgFBO;}
    void rotateAroundAxis(const int); 
    void setMouseRot(const float x,const float y, const float z) {
      x_mouse = (int) y;
      y_mouse = (int) x;
      z_mouse = (int) z;
    }
    // select area
    GLSelection * gl_select;
    static void checkGLErrors(std::string s);
    // color bar
    GLColorbar * gl_colorbar;
    
  signals:
    void sigKeyMouse(const bool, const bool);
    void sigScreenshot();
    void leaveEvent();
    void sigMouseXY(const int x, const int y);
public slots:
   void  update(ParticlesData   * ,
                ParticlesObjectVector * ,
                GlobalOptions         * ,
                const bool update_old_obj=true);
   void  update(ParticlesObjectVector * );
   void  update();
   void  updateVbo(const int);   
   void  updateBoundaryPhys(const int, const bool);
   void  updateColorVbo(const int);
   void  changeColorMap();
   void  reverseColorMap();
   void  rebuildGrid(bool ugl=true);
   void  updateGrid(bool ugl=true);
   void  updateGL();
   void forcePaintGL() {
       makeCurrent();
       paintGL();
   }

   void  osdZoom(bool ugl=true);
   void setOsd(const GLObjectOsd::OsdKeys k,const QString text, bool show,bool b=true);
   void setOsd(const GLObjectOsd::OsdKeys k,const int value, bool show, bool b=true);
   void setOsd(const GLObjectOsd::OsdKeys k,const float value, bool show, bool b=true);
   void setOsd(const GLObjectOsd::OsdKeys k, const float value1, 
                      const float value2, const float value3,  bool show,bool b=true);
   void changeOsdFont();     
   void toggleRotateScreen() {
     
     last_posx = last_posy = last_posz =0;
     if (store_options->rotate_screen) {
       y_mouse = store_options->xrot;
       x_mouse = store_options->yrot;
       z_mouse = store_options->zrot;
     } else {
       y_mouse = store_options->urot;
       x_mouse = store_options->vrot;
       z_mouse = store_options->wrot;     
     }
   }

   void resetFrame() { nframe=0; }
   int getFrame() { return nframe;}

protected:
  void	initializeGL();
  void	paintGL();
  void	resizeGL( int w, int h );

private slots:
  void updateVel(const  int); // update velocity vector
  void updateIpvs(const int ipvs=-1) {
    p_data->setIpvs(ipvs);
    gl_colorbar->update(&gpv,p_data->getPhysData(),store_options,mutex_data);
  }
  void leaveEvent ( QEvent * event ) {
      if (event) {;}
      emit leaveEvent();
    }
  void rotateAroundX() { rotateAroundAxis(0);}
  void rotateAroundY() { rotateAroundAxis(1);}
  void rotateAroundZ() { rotateAroundAxis(2);}
  void rotateAroundU() { rotateAroundAxis(3);}
  void rotateAroundV() { rotateAroundAxis(4);}
  void rotateAroundW() { rotateAroundAxis(5);}
  
  void translateX()    { translateAlongAxis(0); }
  void translateY()    { translateAlongAxis(1); }
  void translateZ()    { translateAlongAxis(2); }
  void setTextureObject(const int, const int);
  void resetEvents(bool pos=false);
  void resetMatScreen() {
    memcpy(mScreen,mIdentity, 16*sizeof(GLdouble));
  }
  void resetMatScene() {
    memcpy(mScene,mIdentity, 16*sizeof(GLdouble));
  }

private:
  // my parent
  QWidget * parent;
  // global options
  GlobalOptions * store_options;
  // grid variables
  GLGridObject * gridx, * gridy, * gridz;
  GLCubeObject * cube;
  // axes
  GLAxesObject * axes;
  // Vectors
  GLObjectParticlesVector gpv;
  ParticlesObjectVector * pov;
  ParticlesData   * p_data;
  // projections
  void setProjection(const int x, const int y, const int w, const int h );
  void computeOrthoFactor();
  float ortho_left,ortho_right,ortho_bottom,ortho_top;
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
  float last_xrot, last_yrot, last_zrot;
  float last_urot, last_vrot, last_wrot;
  int   i_umat, i_vmat, i_wmat; // index of the SCENE/Object rotation matrix
  
  void setRotationScreen( const int x, const int y, const int z );
  void setRotationScene ( const int u, const int v, const int w );
  void getPixelTranslation(int *x, int *y, int *z);
  void setTranslation( const int x, const int y, const int z );
  void setZoom(const int z);
  void setZoom(const float z);
  int zoom_dynam;
  // gl matrix
  GLdouble mProj[16], mModel[16], mModel2[16],
  mScreen[16], mScene[16], mRot[16];
  GLdouble static mIdentity[16];
  int viewport[4];
  bool reset_screen_rotation, reset_scene_rotation;
  void setProjMatrix()  {
    glGetDoublev(GL_PROJECTION_MATRIX, (GLdouble *) mProj);
  }
  void setModelMatrix() {
    glGetDoublev(GL_MODELVIEW_MATRIX, (GLdouble *) mModel);
  }
  void setViewPort() {
    glGetIntegerv(GL_VIEWPORT,viewport);
  }
  void setPerspectiveMatrix();
  // OSD
  GLObjectOsd * osd;
  QImage image,gldata;
  // Font
  fntTexFont * font;
  // Font
  //fntRenderer * text;
  // Texture vector
  GLTextureVector gtv;
  // Thread
  QMutex * mutex_data;
  
  bool is_shift_pressed;
  // bench
  int nframe;
  // Shaders
  CShader * shader;
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
