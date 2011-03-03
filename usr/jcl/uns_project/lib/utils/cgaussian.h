// ============================================================================
// Copyright Jean-Charles LAMBERT - 2009                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

/*
  @author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
*/

#ifndef CGAUSSIAN_H
#define CGAUSSIAN_H

namespace jclut {

  template <class T> class CGaussian {
  public:
    CGaussian(const int pixel , const T g );
    ~CGaussian() { delete [] gaussian;}
    T * data() { return gaussian;}
    void applyOnArrayXY(T * tab, const int dimx,
                        const int dimy, const int x, const int y,
                        const T weight=1.0);
  private:
    int pixel;
    T g;
    T * gaussian;
  };
}
#endif
//
  
