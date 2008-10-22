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
#include "particles_data.h"
class GLBox;
class GLHudObject;

#include "gl_hud_object.h"
#include "global_options.h"
#include "animation_engine.h"
#include "gloctree.h"
#include "glcube.h"

class GLBox : public QGLWidget
{
  Q_OBJECT

    public:

  GLBox( QWidget* parent, const char* name,GlobalOptions * _options,
        AnimationEngine * ,
        const QGLWidget* shareWidget=0);
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
  void getData(const ParticlesData *, ParticlesSelectVector*);
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
  // Grids
  
  void toggleGridX() { gridx->toggleActivate(); updateGL();};
  void toggleGridY() { gridy->toggleActivate(); updateGL();};
  void toggleGridZ() { gridz->toggleActivate(); updateGL();};
  void resizeGrid(const float, const int);
  void changeColorGridX(const QColor color) { gridx->setColor(color); updateGL();};
  void changeColorGridY(const QColor color) { gridy->setColor(color); updateGL();};
  void changeColorGridZ(const QColor color) { gridz->setColor(color); updateGL();};
  // Cube
  void toggleCube() { cube->toggleActivate(); updateGL();};
  void changeColorCube(const QColor color) { cube->setColor(color); updateGL();};
  // Hud
  void toggleHUD()   { hud->toggleActivate();   updateGL();};
  void setHudActivate();
  void setHudActivateNoGL();
  void changeColorHUD(const QColor color);
  // divers
  void togglePoly()  { show_poly =! show_poly;  updateGL();};
  void setParticlesSize(int,const  bool ugl=true);
  void bestFit();
  // textures
  void setTextureSize(const float, bool ugl=true );
  void changeTextureAlphaColor(const int alpha, const bool ugl=true);
  // tree
  void treeUpdate(bool ugl=true);
  private slots:
  void updateOptions(GlobalOptions * , const bool);  
  void takeScreenshot(QImage &);
  void updateVelVectorFactor();
  public:
  
  bool statusBlending() { return blending;}
  bool statusDepthBuffer() { return depth_buffer;}
  bool statusHUD() { return hud->getActivate(); };
  
  GLfloat  getParticlesSize() { return particles_size;}

  // return gridXYZ status
  bool statusGridX() { return gridx->getActivate();};
  bool statusGridY() { return gridy->getActivate();};
  bool statusGridZ() { return gridz->getActivate();};
  // set gridXYZ status
  void setGridXActivate(bool status) { gridx->setActivate(status);};
  void setGridYActivate(bool status) { gridy->setActivate(status);};
  void setGridZActivate(bool status) { gridz->setActivate(status);};
  // rebuild grid
  void gridReBuild();
  
  bool statusPoly()  { return show_poly; };
  void toggleLineAliased() { line_aliased =! line_aliased; updateGL();};
  void setHud(const GLHudObject::HudKeys k,const QString text);
  void setHud(const GLHudObject::HudKeys k,const int value);
  void setHud(const GLHudObject::HudKeys k,const float value);
  void setHudToggle(const GLHudObject::HudKeys k);  
  void setHudActivate(const GLHudObject::HudKeys k, const bool status);  
  void hudProjection();
  const double * getProjMatrix() { return ((double *) mProj); };
  const double * getModelMatrix() { return ((double *) mModel); };
  const int    * getViewPort() { return ((int *) viewport); };
  
  // Orthographic projection
  void setOrthoProjection(float);
  
  GLfloat MAX_PARTICLES_SIZE;
  
  protected:
  
  void	initializeGL();
  void	paintGL();
  void	resizeGL( int w, int h );

  private:
  
  static int  width,height; // frame size

  // GLOctree
  GLOctree * tree;
  //
  GLfloat ratio;
  
  // store options
  GlobalOptions * store_options;
  // Animation stuff
  AnimationEngine * anim_engine;
  // zoom/rotation/translation
  float zoom, zoo_power;
  GLfloat xRot, yRot, zRot, scale;
  GLfloat xTrans, yTrans, zTrans;
  int ixTrans, iyTrans, izTrans;
  
  // grid variable
  GLGridObject * gridx, * gridy, * gridz;
  GLCube * cube;
  bool show_grid;
  bool line_aliased;

  // particles object variable
  GLParticlesObject * vparticles_object[50];
  int nb_object;
  GLuint texture[1];// Storage For One Texture
  float u_max,v_max,tratio;

  // particles data
  ParticlesData * p_data;
  
  // polygones
  bool show_poly;
  
  // velocity vector
  float vel_max_norm;
  
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
  void computeOrthoFactorRatio();
  float getVelMaxNorm(const ParticlesData *, ParticlesSelectVector*);
};
#endif // GLBOX_H
