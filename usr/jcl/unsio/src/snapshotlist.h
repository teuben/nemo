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

/**
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */

#ifndef CSNAPHOTLIST_H
#define CSNAPHOTLIST_H

#include <string>
#include <fstream>
#include "snapshotinterface.h"
#ifndef NOSQLITE3
#include "sqlite_tools.h"
#endif

namespace uns {
  class CunsIn;
  class CSnapshotList : public CSnapshotInterfaceIn {

  public:
    CSnapshotList(const std::string, const std::string, const std::string,
		 const bool verb=false);
    ~CSnapshotList();
    int nextFrame(uns::UserSelection &);
    int close() { return 1;}
    ComponentRangeVector * getSnapshotRange();
    bool getData(const std::string name,int *n,float **f) { return snapshot->getData(name,n,f); }
    bool getData(const std::string name,       float * f) { return snapshot->getData(name,  f); }
    bool getData(const std::string name,int *n,int   **i) { return snapshot->getData(name,n,i); }
    bool getData(const std::string name,       int   * i) { return snapshot->getData(name,  i); }
    bool getData(const std::string comp, const std::string name,int *n,float **f) {
      return snapshot->getData(comp,name,n,f);
    }
    bool getData(const std::string comp, const std::string name,int *n,int **i) {
      return snapshot->getData(comp,name,n,i);
    }
    //float    getEps(const std::string);
    std::string getFileName() { 
      if (snapshot) return snapshot->getFileName();
      else return CSnapshotInterfaceIn::getFileName();
    }
    bool isNewFrame();
    bool shift(std::string name, const float x, const float y, const float z) {
      return snapshot->shift(name,x,y,z);
    }

    virtual ComponentRangeVector * getCrvFromSelection() { return snapshot->user_select.getCrvFromSelection();}
    //bool     isNewFrame();
    std::string getFileStructure() {
      if (snapshot) return snapshot->getFileStructure();
      std::cerr << "Algo error : snapshot not defined...\n";
      assert(0);
      return "";
    }

  private:
    // from ascii database
    bool getLine(const bool force=false);
    bool openFileList();
    bool findSim();
    bool fillNemoRange();
    // from SQLite database

    std::string simname, snapname;
    CunsIn *  unsin;
    CSnapshotInterfaceIn * snapshot;
    std::ifstream fi;
    std::string simtype; // gadget, nemo, ftm
    std::string dirname; // sim's dirname   
    std::string basename;// sim's basename  
    int nframe;          // #frames read
    bool buildGadgetFile();
    bool buildNemoFile();
    virtual int nextFrameSelect(ComponentRangeVector * crvs);
    int addNemoComponent(int&,std::string,std::string );
    std::string nemosim;
    ComponentRangeVector crv;
  };
  
} // namespace

#endif
