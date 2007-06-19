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
#ifndef FILELISTDATA_H
#define FILELISTDATA_H
#include <fstream>

#include "virtual_data.h"
class FileListData: public VirtualData {
  Q_OBJECT
 public:
  FileListData(const char *, const char *, const char *, const bool, pthread_mutex_t *);
  ~FileListData();
  int loadPos(ParticlesSelectVector *, const bool );
  int     getNbody()       { return (*(virtual_data->part_data->nbody)); };
  float * getPos()         { return virtual_data->part_data->pos;        };
  float   getTime()        { return (*(virtual_data->part_data->timu));  };
  float * getCooMax()      { return virtual_data->part_data->coo_max;    };
  int   * getCooIndexMax() { return virtual_data->part_data->i_max;      };
  const char  * getDataName()    { return virtual_data->getDataName();   };
  const char  * getDataType()    { return "File list";           };

  bool isValidData();  
  void uploadGlData(ParticlesSelectVector *);
  int reload(ParticlesSelectVector *, const bool);
  QString endOfDataMessage();
  ParticlesData * getParticlesData() { return virtual_data->part_data;}
 signals:
  void loadedData(const ParticlesData *, ParticlesSelectVector *);

 private:
  VirtualData * virtual_data;
  const  char * select_part, * select_time;
  const char * filename;
  bool load_vel;
  char  * nemo_file;
  char  * sel2;

  bool is_open;   // TRUE if file has been open
  bool is_parsed; // TRUE if particle string has been parsed
  bool is_new_data_loaded; 
  pthread_mutex_t * mutex_data;
  std::ifstream fi;
  bool valid;
  std::string snapname,snaprange;
  // method
  void openFile();
  bool getLine();
  int close(); // close snapshot
  static const char * magic;
  private slots:
  void getData(const ParticlesData *, ParticlesSelectVector*);
};

#endif
