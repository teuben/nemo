// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004                                       
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

#include "particles_range.h"
#include "virtual_data.h"
extern "C" {
  
  int io_nemo(const char * , const char *, ...);
  //extern int bswap (void *, int, int);
  
};

#if 0
//>> From $NEMOSRC/kernel/io/filesecret.h
#define CHKSWAP /* allow mixed endian datasets - this is a bit dangerous */
#define SingMagic  ((011<<8) + 0222)            /* singular items */
#define PlurMagic  ((013<<8) + 0222)            /* plural items */
//<< From $NEMOSRC/kernel/io/filesecret.h
#endif

class SnapshotData : public VirtualData
{    
  Q_OBJECT
 public: 
  SnapshotData(const char *, const char *, const char *);
  ~SnapshotData();
  int loadPos(ParticlesRangeVector * prv);
  int getNbody() { return (*nbody); };
  float * getPos() { return pos; };
  float  getTime() { return (*timu); };
  float * getCooMax() {return coo_max;};
  int * getCooIndexMax() {return i_max;}; 
  char * getDataName() { return nemo_file; };
  char * getDataType() { return "Nemo file"; };
  bool isValidData();  
  void uploadGlData(ParticlesRangeVector *);
  int reload(ParticlesRangeVector *);
  QString endOfDataMessage();
 signals:
  void loadedData(const int *, const float *, const ParticlesRangeVector*);

 private:
  const  char * select_part, * select_time;
  char  * nemo_file;
  char  * sel2;

  bool is_open;   // TRUE if file has been open
  bool is_parsed; // TRUE if particle string has been parsed
  bool is_new_data_loaded; 
  
  //ParticlesRangeVector prv; // store Particles Range

  // method
  
  int close(); // close snapshot
};
#endif 
