// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2005                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
//                                                                             
// SnapshotData Class definition                                               
//                                                                             
//                                                                             
// ============================================================================

#ifndef SNAPSHOT_DATA_H
#define SNAPSHOT_DATA_H

#include <qobject.h>

#include "virtual_data.h"
extern "C" {
  int io_nemo(const char * , const char *, ...);
};

class SnapshotData : public VirtualData
{    
  Q_OBJECT
 public: 
  SnapshotData(const char *, const char *, const char *);
  ~SnapshotData();
  int loadPos(ParticlesSelectVector * );
  int getNbody() { return (*nbody); };
  float * getPos() { return pos; };
  float  getTime() { return (*timu); };
  float * getCooMax() {return coo_max;};
  int * getCooIndexMax() {return i_max;}; 
  char * getDataName() { return nemo_file; };
  char * getDataType() { return "Nemo file"; };
  bool isValidData();  
  void uploadGlData(ParticlesSelectVector *);
  int reload(ParticlesSelectVector *);
  QString endOfDataMessage();
 signals:
  void loadedData(const int *, const float *, ParticlesSelectVector *);

 private:
  const  char * select_part, * select_time;
  char  * nemo_file;
  char  * sel2;
  int   * nemobits;

  bool is_open;   // TRUE if file has been open
  bool is_parsed; // TRUE if particle string has been parsed
  bool is_new_data_loaded; 
  // method
  
  int close(); // close snapshot
};
#endif 
// ============================================================================
