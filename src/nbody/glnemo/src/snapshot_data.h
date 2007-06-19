// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2007                                  
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
#include <pthread.h>

#include "virtual_data.h"
extern "C" {
  int io_nemo(const char * , const char *, ...);
};

class SnapshotData : public VirtualData
{    
  Q_OBJECT
 public: 
  SnapshotData(const char *, const char *, const char *, const bool, pthread_mutex_t *);
  ~SnapshotData();
  int loadPos(ParticlesSelectVector *, const bool );
  int     getNbody()       { return (*(part_data->nbody)); };
  float * getPos()         { return part_data->pos;        };
  float   getTime()        { return (*(part_data->timu));  };
  float * getCooMax()      { return part_data->coo_max;    };
  int   * getCooIndexMax() { return part_data->i_max;      };
  const char  * getDataName()    { return nemo_file;             };
  const char  * getDataType()    { return "Nemo file";           };
  bool isValidData();  
  void uploadGlData(ParticlesSelectVector *);
  int reload(ParticlesSelectVector *, const bool);
  QString endOfDataMessage();
 signals:
  void loadedData(const ParticlesData *, ParticlesSelectVector *);

 private:
  const  char * select_part, * select_time;
  bool load_vel;
  char  * nemo_file;
  char  * sel2;

  bool is_open;   // TRUE if file has been open
  bool is_parsed; // TRUE if particle string has been parsed
  bool is_new_data_loaded; 
  pthread_mutex_t * mutex_data;
  // method
  
  int close(); // close snapshot
};
#endif 
// ============================================================================
