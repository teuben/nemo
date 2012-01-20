// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#ifndef GLNEMODENSITYCOLORBAR_H
#define GLNEMODENSITYCOLORBAR_H
#include <QGraphicsScene>
#include <QWidget>

#include "globaloptions.h"
/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
 */
namespace glnemo {


class DensityColorBar: public QGraphicsScene {
  
Q_OBJECT
  public:
    DensityColorBar(GlobalOptions * ,QWidget *parent = 0);

    ~DensityColorBar();
  public slots:
    void draw();
    void draw(const int min, const int max);
    void setGo(GlobalOptions * _go) { go = _go; };
    void resizeEvent ( QResizeEvent * event );
    void clearScene();
  private:
    GlobalOptions * go;
    const QWidget * parent;
};

}

#endif
