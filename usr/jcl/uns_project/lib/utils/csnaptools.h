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

#ifndef CSNAPTOOLS_H
#define CSNAPTOOLS_H
#include <string>
#include <vector>
namespace jclut {

    class CSnaptools {
    public:
      CSnaptools() {;};
      
      template <class T> static void moveToCod(const int nbody,T * pos, T * Vel, T * mass, T * rho, double cod[6], bool move);
      template <class T> static void moveToCom(const int nbody,T * pos, T * mass);
      static std::string basename(const std::string);
      static std::string dirname(const std::string);
      static std::string parseString(std::string & next_string, const std::string sep=",");
      template <class T> static std::vector<T> stringToVector(const std::string s, const int min, T val, std::string sep=",");
    };
}
#endif
//
  
