// ============================================================================
// Copyright Jean-Charles LAMBERT - 2009-2014
//           Centre de donneeS Astrophysiques de Marseille (CeSAM)              
// e-mail:   Jean-Charles.Lambert@lam.fr                                      
// address:  Aix Marseille Universite, CNRS, LAM 
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS UMR 7326                                       
// ============================================================================

/**
   @author Jean-Charles Lambert <Jean-Charles.Lambert@lam.fr>
*/

#ifndef SQLITE_TOOLS_H
#define SQLITE_TOOLS_H
#ifndef NOSQLITE3
#include <sqlite3.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <iterator>

namespace jclt {


  class CSQLite3 {
  private:
    sqlite3 *db;
    char *zErrMsg;
    char **result;
    int rc;
    int nrow,ncol;
    bool db_open;
    
  public:
    
    std::vector<std::string> vcol_head;
    std::vector<std::string> vdata;

    
    CSQLite3 (std::string tablename="init.db");
    int exe(std::string s_exe);
    ~CSQLite3();
    bool isOpen() { return db_open;};
    void display();
  };

}
#endif
#endif
