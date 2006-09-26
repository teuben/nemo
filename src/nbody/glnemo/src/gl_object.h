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
//                                                                             
// GLObject class definition                                                   
//                                                                             
// Manage OpenGL  Object                                                       
// ============================================================================
#ifndef GL_OBJECT_H
#define GL_OBJECT_H

#include <qgl.h>
#include <qcolor.h>


class GLObject : public QGLWidget {

 public:
  GLObject();
  ~GLObject();

  // method
  void display(int my_list=-1);
  void setColor(const QColor&);
  void toggleActivate();
  bool getActivate() { return is_activated; };
  void setActivate(bool status) { is_activated = status ;}; 
  void setWH(int new_w, int new_h) { width=new_w; height=new_h; };
// protected slots:
  void updateAlphaSlot(int);
  
 protected:
  bool is_activated;    // 
  QColor mycolor;
  GLuint dplist_index;
  int width,height;              // Display width and height
  int particles_alpha;  
  // method
  void  buildDisplayList();
};
#endif
// ============================================================================
