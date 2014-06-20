// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008-2014
//           Centre de donneeS Astrophysiques de Marseille (CeSAM)              
// e-mail:   Jean-Charles.Lambert@lam.fr                                      
// address:  Aix Marseille Universite, CNRS, LAM 
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS UMR 7326                                       
// ============================================================================
#ifndef UNSENGINE_H
#define UNSENGINE_H
/**
	@author Jean-Charles Lambert <jean-charles.lambert@lam.fr>
*/

#include <cstring>
#include <string>
#include "snapshotinterface.h"
#include <map>

namespace uns {

const std::string VERSION="1.0.1"; // UNSIO version

inline std::string getVersion() { return uns::VERSION; }

enum StringData {
      // data
      Time      ,
      Redshift  ,
      Pos       ,
      Vel       ,
      Mass      ,
      Id        ,
      Rho       ,
      Hsml      ,
      U         ,  // internal Energy
      Keys      ,
      Aux       ,
      Eps       ,
      Pot       ,
      Acc       ,
      Age       ,
      Temp      ,  // temperature
      Metal     ,  // total metal : gas+stars
      GasMetal  ,  // metal for gas
      StarsMetal,  // metal for stars
      Zs        ,
      ZSMT      ,
      Im        ,
      Cm        ,
      Czs       ,
      Czsmt     ,
      Ssl       ,
      
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
    bool isValid() { return valid;}
    uns::CSnapshotInterfaceIn * snapshot; // object to store data

    // Map to associate component with a type
    static std::map<std::string, int> s_mapCompInt;
    //
    // py wrapper
    //
    int nextFrame(const char *  _bits);
    // float
    bool getData(const std::string  comp,const std::string  prop,
                 int * size,float ** farray);
    bool getData(const std::string  prop,
                 int * size,float ** farray);
    bool getData(const std::string  prop,float * fvalue);

    // int
    bool getData(const std::string  comp,const std::string  prop,
                 int * size,int ** iarray);
    bool getData(const std::string  prop,
                 int * size,int ** iarray);
    bool getData(const std::string  prop,int * ivalue);

    
    // py wrapper

  private:
    void init(const std::string ,const std::string,const std::string, const bool verb=false );
    std::string simname, sel_comp, sel_time; // IN
    void tryGadget();
    void tryNemo();
    void trySimDB();
    void trySnapList();
    void tryRamses();

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
    bool isValid() { return valid;}
    uns::CSnapshotInterfaceOut * snapshot; // object to store data


    // Map to associate the strings with the enum values
    static std::map<std::string, StringData> s_mapStringValues;

    // py wrapper
    // setData FLOAT
    int setData(const std::string  comp,const std::string  prop,
                int  size,float * farray, const bool _addr=false);
    int setData(const std::string  prop,
                int  size,float * farray, const bool _addr=false);
    int setData(const std::string  prop,
                float fvalue);
    // setData INT
    int setData(const std::string  comp,const std::string  prop,
                int  size,int * iarray, const bool _addr=false);
    int setData(const std::string  prop,
                int  size,int * iarray, const bool _addr=false);
    int setData(const std::string  prop,
                int ivalue);

    //
    int save();
    // py wrapper

    static  void initializeStringMap(const bool);
  private:
    std::string simname, simtype;           // OUT
    
    //bool findSim();
    bool valid;
    bool verbose;
  };

    
}

#endif
