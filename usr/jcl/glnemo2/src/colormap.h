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
#ifndef GLNEMOCOLORMAP_H
#define GLNEMOCOLORMAP_H

#include <QObject>
#include <globaloptions.h>
#include <vector>

namespace glnemo {

/**
	@author Jean-Charles Lambert <jean-charles.lambert@lam.fr>
*/
class Colormap: public QObject {
  Q_OBJECT
  public:
    Colormap(GlobalOptions *);

   ~Colormap();
    int load(const int);
    int load();
    void getRGB(std::vector <float> * r,std::vector <float> * g,std::vector <float> *b) {
      r=&R; g=&G; b=&B;
    }
    int getIndexMap() { return cmap; }
  public slots:
    int next();
    int prev();
    void constant() { go->dynamic_cmap = !go->dynamic_cmap;
      emit newColorMap(); }
    void constant(bool b) { go->dynamic_cmap = b;
      emit newColorMap(); }
        void reverse() { go->reverse_cmap = !go->reverse_cmap;
      emit newColorMap(); }
    void reverse(bool b) { go->reverse_cmap = b;
      emit newColorMap(); }
    
  signals:
    void newColorMap();
  private:
    GlobalOptions * go;
    std::vector <float> R,G,B;
    int cmap;
};

}

#endif
