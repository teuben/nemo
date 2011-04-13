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
#ifndef GLNEMODENSITYHISTO_H
#define GLNEMODENSITYHISTO_H
#include <QGraphicsScene>
#include <QWidget>
#include <QResizeEvent>
/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
*/

namespace glnemo {

class DensityHisto : public QGraphicsScene {

Q_OBJECT

  public:
    DensityHisto(QWidget *parent = 0);
    ~DensityHisto();
    void drawDensity(const int _d[100]);
    void resizeEvent ( QResizeEvent * event );
  public slots:
    void drawDensity(int _min=0, int _max=100);
    void clearScene();
  private:
    void drawGrid();

    const QWidget * parent;
    const int border;
    int density_histo[100];
    int maxhisto;
    
    
};

}

#endif
