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
//                                                                             
// GLBox class definition                                                      
//                                                                             
// Manage all OpenGL stuff                                                     
// ============================================================================

#ifndef GLBOX_H
#define GLBOX_H

#include <qgl.h>
#include <qimage.h>
#include "gl_grid_object.h"
#include "gl_particles_object.h"
class GLBox;
class GLHudObject;

#include "gl_hud_object.h"

class GLBox : public QGLWidget
{
  Q_OBJECT

    public:

  GLBox( QWidget* parent, const char* name, const QGLWidget* shareWidget=0,
        const bool _blending=TRUE, const bool _dbuffer=TRUE, const bool _grid=FALSE, 
        const float _psize=1.5);
  ~GLBox();
  float getXtrans() { return xTrans;};
  float getYtrans() { return yTrans;};
  float getZtrans() { return zTrans;};
  void printViewport() {
    GLint	Viewport[4];
    glGetIntegerv(GL_VIEWPORT,Viewport);
    cerr << "Viewport:" << " " << Viewport[0]
                        << " " << Viewport[1]
                        << " " << Viewport[2]
                        << " " << Viewport[3] << "\n";
    
  }
  void myResizeGL(int w, int h) {
    resizeGL(w,h);
    updateGL();
    swapBuffers();
    printViewport();
  }
  void screenshot();

  bool enable_s;


// -------------

  public slots:

  void setRotation( int x, int y, int z );
  void setTranslation( int x, int y, int z ); 
  void setTranslation();
  void getPixelTranslation(int * x, int * y, int * z);
  void setXRotation( int degrees );
  void setYRotation( int degrees );
  void setZRotation( int degrees );
  void setWH(int, int);
  void setZoom(int);
  void setZoom(float);
  float getZoom();
  int getWidth() { return width;};
  int getHeight() { return height;};
  void setProjection(int, int);
  void getData(const int *, const float *, const ParticlesRangeVector*);
  void toggleGrid();
  void toggleBlending() { 
    //cerr << "Blending = " << blending << "\n";
    blending = ! blending;
    updateGL();
      }
  void toggleDepthBuffer() {
    depth_buffer = ! depth_buffer;
    updateGL();
  }
  void toggleGridX() { gridx->toggleActivate(); updateGL();};
  void toggleGridY() { gridy->toggleActivate(); updateGL();};
  void toggleGridZ() { gridz->toggleActivate(); updateGL();};
  void toggleHUD()   { hud->toggleActivate();   updateGL();};
  void togglePoly()  { show_poly =! show_poly;  updateGL();};
  void setParticlesSize(int);
  void bestFit();
  public:
  
  bool statusBlending() { return blending;}
  bool statusDepthBuffer() { return depth_buffer;}
  bool statusHUD() { return hud->getActivate(); };
  
  const GLfloat  getParticlesSize() { return particles_size;}

  bool statusGridX() { return gridx->getActivate();};
  bool statusGridY() { return gridy->getActivate();};
  bool statusGridZ() { return gridz->getActivate();};
  bool statusPoly()  { return show_poly; };
  void toggleLineAliased() { line_aliased =! line_aliased; updateGL();};
  void setHud(const GLHudObject::HudKeys k,const QString text);
  void setHud(const GLHudObject::HudKeys k,const int value);
  void setHud(const GLHudObject::HudKeys k,const float value);
  void setHudToggle(const GLHudObject::HudKeys k);  
  void setHudActivate(const GLHudObject::HudKeys k, const bool status);  
  const double * getProjMatrix() { return ((double *) mProj); };
  const double * getModelMatrix() { return ((double *) mModel); };
  const int    * getViewPort() { return ((int *) viewport); };
  
  GLfloat MAX_PARTICLES_SIZE;
  
  protected:
  
  void	initializeGL();
  void	paintGL();
  void	resizeGL( int w, int h );

  private:
  
  static int  width,height; // frame size

  //
  GLfloat ratio;
  
  // zoom/rotation/translation
  float zoom, zoo_power;
  GLfloat xRot, yRot, zRot, scale;
  GLfloat xTrans, yTrans, zTrans;
  int ixTrans, iyTrans, izTrans;
  
  // grid variable
  GLGridObject * gridx, * gridy, * gridz;
  bool show_grid;
  bool line_aliased;

  // particles object variable
  GLParticlesObject * vparticles_object[50];
  int nb_object;
  GLuint texture[1];// Storage For One Texture
  float u_max,v_max,tratio;
  
  // polygones
  bool show_poly;
  
  // blending options
  GLfloat particles_size;
  bool blending,
       depth_buffer;
       
  // Head Up Display
  GLHudObject * hud;
  
  // gl matrix
  GLdouble mProj[16], mModel[16];
  int viewport[4];
  // Method
  void setProjMatrix()  {
    glGetDoublev(GL_PROJECTION_MATRIX, (GLdouble *) mProj);
  };
  void setModelMatrix() {
    glGetDoublev(GL_MODELVIEW_MATRIX, (GLdouble *) mModel);
  }; 
  void setViewPort() {
    glGetIntegerv(GL_VIEWPORT,viewport);
  };
  void loadImage();
};
#endif // GLBOX_H
