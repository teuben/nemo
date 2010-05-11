// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008-2010                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
#ifndef UNSENGINE_H
#define UNSENGINE_H
/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
*/

#include <cstring>
#include <string>
#include "snapshotinterface.h"
#include <map>

namespace uns {
     enum StringData {
      // data
      Time      ,
      Pos       ,
      Vel       ,
      Mass      ,
      Id        ,
      Rho       ,
      Hsml      ,
      U         ,  // internal Energy
      Keys      ,
      Aux       ,
      Pot       ,
      Acc       ,
      Age       ,
      Temp      ,  // temperature
      Metal     ,  // total metal : gas+stars
      GasMetal  ,  // metal for gas
      StarsMetal,  // metal for stars
      
      // Nbodies  per component
      Nsel      ,
      Nbody     ,
      Ngas      ,
      Nhalo     ,
      Ndisk     ,
      Nbulge    ,
      Nstars    ,
      Nbndry    ,

      // Components
      Gas       ,
      Halo      ,
      Disk      ,
      Bulge     ,
      Stars     ,
      Bndry     ,
      All       ,
      
      GasMPV    ,
      HaloMPV   ,
      DiskMPV   ,
      BulgeMPV  ,
      StarsMPV  ,
      BndryMPV
  };

  // class Cuns                                              
  // manage Unified Nbody Snapshot Input operations            
  class CunsIn {
  public:
    // constructor for READING opertaions
    CunsIn(const char * ,const char * , const char *, const bool verb=false );
    CunsIn(const std::string ,const std::string,const std::string, const bool verb=false );
    ~CunsIn();
    bool isValid() { return valid;};
    CSnapshotInterfaceIn * snapshot; // object to store data

  private:
    void init(const std::string ,const std::string,const std::string, const bool verb=false );
    std::string simname, sel_comp, sel_time; // IN
    void tryGadget();
    void tryNemo();
    void trySimDB();
    void trySnapList();
    
    //bool findSim();
    bool valid;
    bool verbose;
  };

  // class CunsOut                                              
  // manage Unified Nbody Snapshot Output operations            
  class CunsOut {
  public:

    // constructor for WRITING operations
    CunsOut(const std::string, const std::string, const bool verb=false);
    ~CunsOut();
    bool isValid() { return valid;};
    CSnapshotInterfaceOut * snapshot; // object to store data


    // Map to associate the strings with the enum values
    static std::map<std::string, StringData> s_mapStringValues;

    static  void initializeStringMap(const bool);
  private:
    std::string simname, simtype;           // OUT
    
    //bool findSim();
    bool valid;
    bool verbose;
  };

    
}

#endif
