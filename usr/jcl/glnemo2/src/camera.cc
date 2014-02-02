// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2014                                  
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
#include <QtOpenGL>
#include <QFile>
#include <QString>

#include <sstream>
#include "camera.h"
#include <GL/glu.h>

namespace glnemo {
  // ============================================================================
  // constructor                                                                 
  Camera::Camera()
  {  
    ex=ey=ez=0.0;     // eyes
    cx=cy=cz=0.0;     // center 
    ux=uz=0.0;uy=1.0; // up vectors
    
    spline = new CRSpline();
    display_ctrl = false; // toggle display ctrl  
    display_path = false; // toggle display path  
    play         = false; // toggle animation     
    play_timer   = new QTimer(this);
    connect(play_timer,SIGNAL(timeout()),this,SLOT(playGL())); // update GL at every timeout()
    index_frame  = 0;
  }
  // ============================================================================
  // destructor                                                                 
  Camera::~Camera()
  {
    if (spline) delete spline;
    glDeleteLists( dplist_index, 1 );
  }
  // ============================================================================
  // destructor                                                                 
  void Camera::init(std::string filename, const int _p, const float _s)
  {
    npoints = _p; // #interpolated points
    scale   = _s; // scale factor applied on top of ctrl points
    dplist_index = glGenLists( 1 );    // get a new display list index
    if (filename != "") {
      loadSplinePoints(filename);
    }
  }
  // ============================================================================
  // setEye                                                                      
  void Camera::setEye(const float x, const float y, const float z)
  {
    ex=x;ey=y;ez=z;
  }
  // ============================================================================
  //  setCenter                                                                  
  void Camera::setCenter(const float x, const float y, const float z)
  {
    cx=x;cy=y;cz=z;
  }
  // ============================================================================
  //  setUp                                                                      
  void Camera::setUp(const float x, const float y, const float z)
  {
    ux=x;uy=y;uz=z;
  }
  // ============================================================================
  //  moveTo                                                                     
  void Camera::moveTo()
  {
    if (!play) {
      gluLookAt(ex, ey, ez,
                cx, cy, cz,
                ux, uy, uz);
    }
    else {
      index_frame = index_frame%npoints;
      //std::cerr << "frame : "<< index_frame << "\n";
      float  t=(float)index_frame / (float)npoints;
      Vec3D rv = spline->GetInterpolatedSplinePoint(t)*scale;   
      gluLookAt(rv.x, rv.y, rv.z,
                cx, cy, cz,
                ux, uy, ez); // ez, why ????!!!!!!
      index_frame++;
    }
  }
  // ============================================================================
  //  loadSplinePoints                                                           
  int Camera::loadSplinePoints(std::string filen)
  {
    QFile infile(QString(filen.c_str()));
    if (!infile.open(QIODevice::ReadOnly | QIODevice::Text))
      return 0;
    // clear ctrl points
    spline->clearCPoints();
        
    QTextStream in(&infile);
    QString line;
    do {
      line = in.readLine();
      if (!line.isNull()) {
        //std::cerr << "line :" << line.toStdString() <<"\n";
        std::istringstream ss(line.toStdString());
        float x,y,z;
        ss >> x;
        ss >> y;
        ss >> z;
        Vec3D v(x,y,z);
        spline->AddSplinePoint(v);      
      }
    } while (!line.isNull());
    infile.close();
    buildDisplayList();
    return 1;
  }
  // ============================================================================
  //  buildDisplayList                                                           
  void Camera::buildDisplayList()
  {
    glNewList( dplist_index, GL_COMPILE );
    glBegin(GL_LINE_STRIP);
    for (int i=0; i<npoints; i++) {
      float  t=(float)i / (float)npoints;
      Vec3D rv = spline->GetInterpolatedSplinePoint(t)*scale;    
      glVertex3f(rv.x, rv.y, rv.z);  
    }
    glEnd();
    glEndList();
  }
  // ============================================================================
  // displayCameraPath                                                           
  void Camera::displayCameraPath()
  {
    if (spline->GetNumPoints() > 0 ) {
      glEnable (GL_LINE_SMOOTH);
      glEnable (GL_BLEND);
      glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
      glLineWidth (1.5);
      GLObject::display(dplist_index);
      glDisable(GL_BLEND);
    }
  }
  // ============================================================================
  // display                                                                     
  // display ctrl and camera path                                                
  void Camera::display()
  {
    if (display_path) displayCameraPath();
    if (display_ctrl) {;}
  }
  // ============================================================================
  // setSplineParam                                                              
  // set Spline Parameters                                                       
  // npoints : # interpolated points in the spline                               
  // scale   : scale factor applied on top of each control points                
  void Camera::setSplineParam(const int _p, const double _s, const bool ugl)
  { 
    npoints=_p; 
    scale  =_s;
    buildDisplayList();
    if (ugl) { // updateGL required
      emit updateGL();
    }
  }
  // ============================================================================
  // setCamDisplay                                                               
  // toggle ctrl and path camera display                                         
  void Camera::setCamDisplay(const bool ctrl, const bool path, const bool ugl)
  {
    display_ctrl = ctrl;
    display_path = path;
    if (ugl) { // updateGL required
      emit updateGL();
    }
  }
  // ============================================================================
  // startStopPlay                                                               
  void Camera::startStopPlay()
  {
    if (!play) {
      play_timer->start(10);
    }
    else {
      play_timer->stop();
    }
    play = !play;
  }
}
