// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2011                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
*/
#ifndef SNAPSHOTPHIBAR_H
#define SNAPSHOTPHIBAR_H
#include <QObject>
#include <string>
#include "snapshotinterface.h"
#include "zlib.h"
#include "globaloptions.h"

namespace glnemo {
  using namespace std;
class SnapshotPhiGrape: public SnapshotInterface
{
    Q_OBJECT
    Q_INTERFACES(glnemo::SnapshotInterface)

public:
    SnapshotPhiGrape();

    ~SnapshotPhiGrape();
    
	SnapshotInterface * newObject(const std::string _filename, const int x=0);
    ComponentRangeVector * getSnapshotRange();
    bool isValidData();
    int initLoading(GlobalOptions * so);
    int nextFrame(const int * index_tab, const int nsel);
    int close();
    QString endOfDataMessage();
 private:
    bool valid;
    gzFile           file;               // file handle for compressed file
    bool detectHeader();
    int full_nbody, frame_number;

  // Buffer stuffs
  char * BUFF;
  char line[300];
  unsigned int size_buff;
  std::string sbuff; // string to store the buffer

  // ============================================================================
  // getLine()
  inline bool getLine()
  {
    bool status = false;
    size_t found=sbuff.find('\n'); // look for \n
    if (found!=string::npos) { // \n found
      memcpy(line,sbuff.c_str(),found); // build new line
      line[found]='\0';       // null character terminaison
      sbuff.erase(0,found+1); // erase from beginning to \n
      status=true;
    }
    return status;
  }
  bool readBuffer();
  bool gzGetLine();
};

}
#endif // SNAPSHOTPHIBAR_H
