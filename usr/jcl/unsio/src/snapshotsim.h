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

/**
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */

#ifndef SNAPSHOTSIM_H
#define SNAPSHOTSIM_H
#ifndef NOSQLITE3
#include <string>
#include <fstream>
#include "snapshotinterface.h"
#include "sqlite_tools.h"

namespace uns {

  class CSnapshotSimIn : public CSnapshotInterfaceIn {

  public:
    CSnapshotSimIn(const std::string, const std::string, const std::string,
		 const bool verb=false);
    ~CSnapshotSimIn();
    int nextFrame(const uns::t_indexes_tab * index_tab, const int nsel);
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
    float    getEps(const std::string);
    std::string getFileName() { 
      if (snapshot) return snapshot->getFileName();
      else return CSnapshotInterfaceIn::getFileName();
    }
    std::string getFileStructure() {
      if (snapshot) return snapshot->getFileStructure();
      std::cerr << "Algo error : snapshot not defined...\n";
      assert(0);
    }
    bool shift(std::string name, const float x, const float y, const float z) {
      return snapshot->shift(name,x,y,z);
    }
    bool     isNewFrame();
    virtual ComponentRangeVector * getCrvFromSelection() { return snapshot->user_select.getCrvFromSelection();}
    virtual int     getCod(const std::string select, const float time, 
			   float * tcod, const std::string base="ANALYSIS/cod",
			   const std::string ext="cod");
    virtual std::string getSimDir() { return dirname; }
  private:
    // from ascii database
    bool getLine();
    bool openDbFile();
    bool readEpsFile();
    bool eps_exist;
    bool findSim();
    bool fillNemoRange();
    // from SQLite database
    jclt::CSQLite3 * sql;  // sqlite3 object
    bool openSqlDb(std::string db="/pil/programs/DB/simulation.dbl");
    bool readSqlEps();
    bool findSqlSim();
    bool fillSqlNemoRange();
    std::string sqlite_db;

    std::string simname;
    CSnapshotInterfaceIn * snapshot;
    std::ifstream fi;
    std::string simtype; // gadget, nemo, ftm
    std::string dirname; // sim's dirname   
    std::string basename;// sim's basename  
    int nframe;          // #frames read
    bool buildGadgetFile();
    bool buildNemoFile();
    int addNemoComponent(int&,std::string,std::string );
    virtual int nextFrameSelect(ComponentRangeVector * crvs);
    std::string nemosim;
    ComponentRangeVector crv;
  };
  
} // namespace
#endif
#endif
