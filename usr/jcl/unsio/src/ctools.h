// ============================================================================
// Copyright Jean-Charles LAMBERT - 2009-2010                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
#ifndef CTOOLS_H
#define CTOOLS_H
/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
*/
#include <cctype>
#include <cstring>
#include <string>

namespace tools {

// Physical bits
#define TIME_BIT          (1 <<  1)
#define KEYS_BIT          (1 <<  2)
#define HEADER_BIT        (1 <<  3)
#define MASS_BIT          (1 <<  4)
#define POS_BIT           (1 <<  5)
#define VEL_BIT           (1 <<  6)
#define EPS_BIT           (1 <<  7)
#define RHO_BIT           (1 <<  8)
#define HSML_BIT          (1 <<  9)
#define U_BIT             (1 << 10)
#define ID_BIT            (1 << 11)
#define METAL_BIT         (1 << 12)
#define AGE_BIT           (1 << 13)
#define AUX_BIT           (1 << 14)
#define POT_BIT           (1 << 15)
#define ACC_BIT           (1 << 16)
#define TEMP_BIT          (1 << 17)
  
// Component BITS
#define ALL_BIT           (1 <<  1)
#define GAS_BIT           (1 <<  2)
#define HALO_BIT          (1 <<  3)
#define DISK_BIT          (1 <<  4)
#define BULGE_BIT         (1 <<  5)
#define STARS_BIT         (1 <<  6)
#define BNDRY_BIT         (1 <<  7)

  

  class Ctools {
  public:
    Ctools() {;}
    static std::string fixFortran(const char *);
    static bool isFileExist(std::string);
    static std::string toupper(std::string);
    static std::string tolower(std::string);
    static int compBits(std::string s) {
      int ret=0;
      if (s=="all"   ) ret=ALL_BIT|GAS_BIT|HALO_BIT|DISK_BIT|STARS_BIT|BNDRY_BIT;
      if (s=="gas"   ) ret=GAS_BIT;
      if (s=="halo"  ) ret=HALO_BIT;
      if (s=="disk"  ) ret=DISK_BIT;
      if (s=="stars" ) ret=STARS_BIT;
      if (s=="bndry" ) ret=BNDRY_BIT;
      return ret;
    }
  };
}

#endif
