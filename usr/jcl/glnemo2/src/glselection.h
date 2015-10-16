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
	@author Jean-Charles Lambert <jean-charles.lambert@lam.fr>
*/
#ifndef GLNEMOGLSELECTION_H
#define GLNEMOGLSELECTION_H
#include <QObject>
#include <QMouseEvent>
#include <QTimer>
#include <QMutex>
#include <globaloptions.h>
#include <particlesdata.h>
#include "globjectparticles.h"
#include "vec3d.h"
#include <QGLWidget>

namespace glnemo {

class GLSelection : public QObject
{
  Q_OBJECT

public:
    GLSelection();
    ~GLSelection();
    void reset();
    void update(const GLObjectParticlesVector *,
                GlobalOptions   *, QMutex * );
    bool isEnable()   { return enable;}
    void setEnable(bool _b) { enable=_b;    }
    void getMouse(QMouseEvent *);
    void display(const int, const int);
    void selectOnArea(const int nobj, double mProj[16],double mModel[16],
                    const int viewport[4], const bool whole_screen=false);
    void zoomOnArea(double mProj[16],double mModel[16],
                    const int viewport[4]);
    float X0() { return x0;}
    float X1() { return x1;}
    float Y0() { return y0;}
    float Y1() { return y1;}
    std::vector <int> * getList() { return &list; }
public slots:
 void setZoom(bool _b)     { zoom      = _b; }
 void setAnimZoom(bool _b) { anim_zoom = _b; }

private slots:
  void playZoomAnim();
signals:
  void updateGL();
  void updateZoom();
  void updatePareticlesSelected(const int);
  void sendList(std::vector <int> &);
private:
  std::vector <int> list;      // to save particles
  float x0,y0,x1,y1;
  bool enable;
  bool zoom;
  bool anim_zoom;
  const GLObjectParticlesVector * gpv;
  GlobalOptions * store_options;
  QTimer * anim_timer;
  Vec3D trans_in, trans_out, comvec;
  double zoom_dynamic, zoomo_dynamic;
  int total_frame, frame_counter;
  QMutex * mutex_data;
  int in_area;
  double com[3];
  const ParticlesData * part_data;
  // METHOD
   
};

}

#endif
