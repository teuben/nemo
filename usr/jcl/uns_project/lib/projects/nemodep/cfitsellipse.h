// ============================================================================
// Copyright Jean-Charles LAMBERT - 2011                                       
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
#ifndef CFITSELLIPSE_H
#define CFITSELLIPSE_H
#include <string>
namespace uns_proj {
  class CFitsEllipse {

  public:
    CFitsEllipse(const int xaxis, const int yaxis, int nm=1000, float tmax=50.0);
    ~CFitsEllipse();
    
    void buildGrid(const int nbody, const float * x, const float * val);
    void displayGrid();
    void saveGrid(std::string);
    float intensity(const float x, const float y);
  private:
    float tmax;     // half_size of te grid
    int   nmesh;    // #mesh in the grid
    float * grid;   // Ze grid
    int xaxis,yaxis;
  };
  

  
}

#endif // CFITSELLIPSE_H
