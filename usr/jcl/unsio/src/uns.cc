// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008-2013
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

#include "uns.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <assert.h>
#include "snapshotinterface.h"
#include "snapshotgadget.h"
#include "snapshotramses.h"
#include "snapshotnemo.h"
#include "snapshotsim.h"
#include "snapshotlist.h"
#include "userselection.h"
#include "ctools.h"

#define DEBUG 0
#include "unsdebug.h"

namespace uns {
  // static variable to store DATA string
  std::map<std::string, StringData> CunsOut::s_mapStringValues;
  std::map<std::string, int> CunsIn::s_mapCompInt;
  // ----------------------------------------------------------------------------
  // READING OPERATIONS
  CunsIn::CunsIn(const std::string _name ,const std::string _comp ,const std::string _time, const bool verb)
  {
    init(_name,_comp,_time,verb);
  }

  // ----------------------------------------------------------------------------
  // constructor for READING operations
  CunsIn::CunsIn(const char * _name ,const char * _comp, const char * _time,
                 const bool verb)
  {
    init(_name,_comp,_time,verb);
  }
  // ----------------------------------------------------------------------------
  //
  void CunsIn::init(const std::string _name ,const std::string _comp ,const std::string _time, const bool verb )
  {
    valid = false;
    simname  = tools::Ctools::fixFortran(_name.c_str(),false);
    sel_comp = tools::Ctools::fixFortran(_comp.c_str(),false);
    sel_time = tools::Ctools::fixFortran(_time.c_str(),false);
    
    // initialise some maps
    CunsIn::s_mapCompInt["gas"        ] = 0;
    CunsIn::s_mapCompInt["halo"       ] = 1;
    CunsIn::s_mapCompInt["dm"         ] = 1;
    CunsIn::s_mapCompInt["disk"       ] = 2;
    CunsIn::s_mapCompInt["bulge"      ] = 3;
    CunsIn::s_mapCompInt["stars"      ] = 4;
    CunsIn::s_mapCompInt["bndry"      ] = 5;
    CunsIn::s_mapCompInt["all"        ] =-1;

    // to lower
    //simname  = tools::Ctools::tolower(simname);
    //sel_comp = tools::Ctools::tolower(sel_comp);
    //sel_time = tools::Ctools::tolower(sel_time);
    
    verbose=verb;
    snapshot = NULL;
    PRINT("name    ["<< simname  <<"]\n");
    PRINT("sel_comp["<< sel_comp <<"]\n");
    PRINT("sel_time["<< sel_time <<"]\n");
    CunsOut::initializeStringMap(verbose);
    if (simname == "-") { // we assume here that "-"
      tryNemo();          // is standard input and a NEMO stream...      
    } else {
      if (tools::Ctools::isFileExist(simname)) {  // file exist 
        tryGadget();               // try gadget 
        if (!valid) {              // gadget failed
          tryRamses();             // try ramses
        }
        if (!valid) {              // gadget failed
          tryNemo();               // try nemo   
        }
        if (!valid) {              // nemo 
          trySnapList();           // try snapshotlist   
        }
        if (!valid) {
          trySimDB();              // try DataBase       
        }
      }
      else {                       // file does not exist
        tryGadget();               // try gadget parallel output
        if (!valid) {
          trySimDB();              // try DataBase       
        }
      }
    }
    if (valid && verb) {
      std::cerr << "File      : " << snapshot->getFileName() << "\n";
      std::cerr << "Interface : " << snapshot->getInterfaceType() << "\n";
    }
  }
  // ----------------------------------------------------------------------------
  // destructor for READING operations  
  CunsIn::~CunsIn()
  {
    if (snapshot) delete snapshot;
  }
  
  // ----------------------------------------------------------------------------
  // tryGadget                                                                   
  void CunsIn::tryGadget()
  {
    PRINT("tryGadget("<<simname<<")\n");
    snapshot = new CSnapshotGadgetIn(simname, sel_comp, sel_time, verbose);
    valid = snapshot->isValidData();
  }
  // ----------------------------------------------------------------------------
  // tryRamses
  void CunsIn::tryRamses()
  {
    PRINT("tryRamses("<<simname<<")\n");
    snapshot = new CSnapshotRamsesIn(simname, sel_comp, sel_time, verbose);
    valid = snapshot->isValidData();
  }
  // ----------------------------------------------------------------------------
  // tryNemo                                                                     
  void CunsIn::tryNemo()
  {
    PRINT("tryNemo("<<simname<<")\n");
    snapshot = new CSnapshotNemoIn(simname, sel_comp, sel_time, verbose);
    valid = snapshot->isValidData();
  }
  // ----------------------------------------------------------------------------
  // trySim                                                                      
  void CunsIn::trySimDB()
  {
#ifndef NOSQLITE3  
    snapshot = new CSnapshotSimIn(simname, sel_comp, sel_time, verbose);
    valid = snapshot->isValidData();
    if (valid && verbose) {
      std::cerr << "CunsIn::trySimDB() It's a simulation...\n";
    }
#else
    valid = false;
#endif
  }
  // ----------------------------------------------------------------------------
  // trySnapList                                                                      
  void CunsIn::trySnapList()
  {
    snapshot = new CSnapshotList(simname, sel_comp, sel_time, verbose);
    valid = snapshot->isValidData();
  }
  // ----------------------------------------------------------------------------
  // WRITING OPERATIONS
  //  ---------------------------------------------------------------------------
  // constructor
  CunsOut::CunsOut(const std::string _name, const std::string _type, const bool _verb )
  {
    simname  = tools::Ctools::fixFortran(_name.c_str(),false);
    simtype  = tools::Ctools::fixFortran(_type.c_str(),false);
    verbose = _verb;
    snapshot= NULL;
    initializeStringMap(verbose);
    simtype = tools::Ctools::tolower(simtype);
    if (simtype == "gadget2") {    
      snapshot = new CSnapshotGadgetOut(simname,simtype,verbose);
    } else {
      if (simtype == "nemo") {
        snapshot = new CSnapshotNemoOut(simname,simtype,verbose);
      }
      else {
        std::cerr << "Unkonwn UNS output file format => ["<<simtype<<"]"
            << " aborting program...... \n\n";
        std::exit(1);      
      }
    }
  }
  // ----------------------------------------------------------------------------
  // destructor for READING operations  
  CunsOut::~CunsOut()
  {
    if (snapshot) delete snapshot;
  }
  // ----------------------------------------------------------------------------
  // initializeStringMap
  void  CunsOut::initializeStringMap(const bool verbose)
  {
    
    CunsOut::s_mapStringValues["time"       ] = uns::Time;
    CunsOut::s_mapStringValues["redshift"   ] = uns::Redshift;
    CunsOut::s_mapStringValues["pos"        ] = uns::Pos;
    CunsOut::s_mapStringValues["vel"        ] = uns::Vel;
    CunsOut::s_mapStringValues["mass"       ] = uns::Mass;
    CunsOut::s_mapStringValues["id"         ] = uns::Id;
    CunsOut::s_mapStringValues["rho"        ] = uns::Rho;
    CunsOut::s_mapStringValues["hsml"       ] = uns::Hsml;
    CunsOut::s_mapStringValues["u"          ] = uns::U;
    CunsOut::s_mapStringValues["aux"        ] = uns::Aux;
    CunsOut::s_mapStringValues["acc"        ] = uns::Acc;
    CunsOut::s_mapStringValues["pot"        ] = uns::Pot;
    CunsOut::s_mapStringValues["keys"       ] = uns::Keys;
    CunsOut::s_mapStringValues["age"        ] = uns::Age;
    CunsOut::s_mapStringValues["temp"       ] = uns::Temp;
    CunsOut::s_mapStringValues["metal"      ] = uns::Metal;
    CunsOut::s_mapStringValues["gas_metal"  ] = uns::GasMetal;
    CunsOut::s_mapStringValues["stars_metal"] = uns::StarsMetal;
    CunsOut::s_mapStringValues["nsel"       ] = uns::Nsel;
    CunsOut::s_mapStringValues["nbody"      ] = uns::Nbody;
    CunsOut::s_mapStringValues["ngas"       ] = uns::Ngas;
    CunsOut::s_mapStringValues["nhalo"      ] = uns::Nhalo;
    CunsOut::s_mapStringValues["ndisk"      ] = uns::Ndisk;
    CunsOut::s_mapStringValues["nbulge"     ] = uns::Nbulge;
    CunsOut::s_mapStringValues["nstars"     ] = uns::Nstars;
    CunsOut::s_mapStringValues["nbndry"     ] = uns::Nbndry;
    CunsOut::s_mapStringValues["gas"        ] = uns::Gas;
    CunsOut::s_mapStringValues["halo"       ] = uns::Halo;
    CunsOut::s_mapStringValues["dm"         ] = uns::Halo;
    CunsOut::s_mapStringValues["ndm"        ] = uns::Halo;
    CunsOut::s_mapStringValues["bulge"      ] = uns::Bulge;
    CunsOut::s_mapStringValues["disk"       ] = uns::Disk;
    CunsOut::s_mapStringValues["stars"      ] = uns::Stars;
    CunsOut::s_mapStringValues["bndry"      ] = uns::Bndry;
    CunsOut::s_mapStringValues["all"        ] = uns::All;
    CunsOut::s_mapStringValues["gas_mpv"    ] = uns::GasMPV;
    CunsOut::s_mapStringValues["halo_mpv"   ] = uns::HaloMPV;
    CunsOut::s_mapStringValues["bulge_mpv"  ] = uns::BulgeMPV;
    CunsOut::s_mapStringValues["disk_mpv"   ] = uns::DiskMPV;
    CunsOut::s_mapStringValues["stars_mpv"  ] = uns::StarsMPV;
    CunsOut::s_mapStringValues["bndry_mpv"  ] = uns::BndryMPV;
    CunsOut::s_mapStringValues["zs"         ] = uns::Zs;
    CunsOut::s_mapStringValues["zsmt"       ] = uns::ZSMT;
    CunsOut::s_mapStringValues["im"         ] = uns::Im;
    CunsOut::s_mapStringValues["ssl"        ] = uns::Ssl;
    CunsOut::s_mapStringValues["cm"         ] = uns::Cm;
    CunsOut::s_mapStringValues["czs"        ] = uns::Czs;
    CunsOut::s_mapStringValues["czsmt"      ] = uns::Czsmt;
    if (verbose) {
      std::cout << "CunsOut::initializeStringMap s_mapStringValues contains "
          << CunsOut::s_mapStringValues.size() << " entries." << std::endl;
    }
  }
  
  
}
// ---------------------------------------------------------------------------- 
// End of file                                                                 
// ---------------------------------------------------------------------------- 
