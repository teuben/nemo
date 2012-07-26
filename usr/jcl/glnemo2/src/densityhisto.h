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
#ifndef GLNEMODENSITYHISTO_H
#define GLNEMODENSITYHISTO_H
#include <QGraphicsScene>
#include <QWidget>
#include <QResizeEvent>
#include <vector>
/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
*/

namespace glnemo {

class DensityHisto : public QGraphicsScene {

Q_OBJECT

  public:
    DensityHisto(QWidget *parent = 0);
    ~DensityHisto();
    void drawDensity(const std::vector<int> _d); // _d[100]);
    void resizeEvent ( QResizeEvent * event );
  public slots:
    void drawDensity(int _min=0, int _max=100);
    void clearScene();
  private:
    const int nhisto;
    void drawGrid();

    const QWidget * parent;
    const int border;
    std::vector <int> density_histo;
    int maxhisto;
    
    
};

}

#endif
