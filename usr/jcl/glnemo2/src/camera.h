// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2011                                  
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

#ifndef CAMERA_H
#define CAMERA_H
#include <QObject>
#include <QTimer>
#include <fstream>

#include "catmull_rom_spline.h"
#include "globject.h"

namespace glnemo {

/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
*/
  class Camera: public GLObject{
    Q_OBJECT
  public:
    Camera();
    ~Camera();
    // init
    void init(std::string filename="", const int _p=2000, const float _s=1.0); // MUST BE CALLED AFTER MAKECURRENT()
    
    // positions and mouvement
    void setEye(const float, const float, const float);
    void setCenter(const float, const float, const float);
    void setUp(const float, const float, const float);
    void moveTo();
    // spline
    int loadSplinePoints(std::string file);
    
    // SLOTS
  public slots:
    void setSplineParam(const int ,const double, const bool ugl=true);  // set spline parameters
    void setCamDisplay(const bool , const bool, const bool ugl=true);  // toggle camera display
    void startStopPlay();
    void display();      
    
  signals:
    void updateGL();
    
  private slots:
    void playGL() { emit updateGL();}

  private:
    float
        ex, ey, ez,
        cx, cy, cz,
        ux, uy, uz;
    
    CRSpline * spline;
    // playing timer
    QTimer * play_timer;
    int npoints; // # interpolated points
    float scale; // Scale factor applied on control points
    int index_frame; 
    bool display_ctrl, display_path;
    bool play;   // true if playing
    void buildDisplayList();
    void displayCameraPath();
    
  };
}
#endif // CAMERA_H
