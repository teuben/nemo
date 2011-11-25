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
/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
*/
#ifndef GLCOLORBAR_H
#define GLCOLORBAR_H
#include <QObject>
#include <QMouseEvent>
#include <QTimer>
#include <QMutex>
#include <globaloptions.h>
#include <particlesdata.h>
#include "globjectparticles.h"
#include "vec3d.h"
#include "gltextobject.h"

#include <QGLWidget>

namespace glnemo {


class GLColorbar : public GLObject
{
  Q_OBJECT
public:
    GLColorbar(const GlobalOptions *, 
               bool activated=true);
    ~GLColorbar();
    bool isEnable()   { return is_activated;}
    void setEnable(bool _b) { is_activated=_b;    }
    void update( GLObjectParticlesVector *,PhysicalData * phys_select,
                GlobalOptions   *, QMutex * );      
    void display(const int, const int);  
public slots:
    void updateFont();
private:
    const GLObjectParticlesVector * gpv;
    const GlobalOptions * go;
    QMutex * mutex_data;
    bool is_activated;
    int width,height;
    void drawBox  ();
    void drawColor();
    void drawLegend();
    void drawText(float value, int fac);
    // font stuffs
    GLTextObject * legend;    
    fntTexFont * font;
    int x[4][2];
        
    PhysicalData * phys_select;

};


}
#endif // GLCOLORBAR_H
