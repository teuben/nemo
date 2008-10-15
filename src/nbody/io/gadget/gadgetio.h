// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
*/
#ifndef GADGETGADGETIO_H
#define GADGETGADGETIO_H
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include "componentrange.h"
#include "userselection.h"

namespace gadget {

int compare( const void * elem1, const void * elem2 );

enum ioop { READ,WRITE };

typedef struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} t_io_header_1;

typedef struct s_data_type_header {
  int size_bytes;
  int items;
} t_data_type_header;

enum DATA_TYPE {
  INT   =4,
  FLOAT =4,
  DOUBLE=8,
  CHAR  =1
};

const int NB_DATA_HEADER = 14;
const t_data_type_header dth[NB_DATA_HEADER] = {
  {INT      ,6}, // int      npart[6];
  {DOUBLE   ,6}, // double   mass[6];
  {DOUBLE   ,1}, // double   time;
  {DOUBLE   ,1}, // double   redshift;
  {INT      ,1}, // int      flag_sfr;
  {INT      ,1}, // int      flag_feedback;
  {INT      ,6}, // int      npartTotal[6];
  {INT      ,1}, // int      flag_cooling;
  {INT      ,1}, // int      num_files;
  {DOUBLE   ,1}, // double   BoxSize;
  {DOUBLE   ,1}, // double   Omega0;
  {DOUBLE   ,1}, // double   OmegaLambda;
  {DOUBLE   ,1}, // double   HubbleParam; 
  {CHAR     ,96} // char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];
};

typedef struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;
  int    Id;     // added by JCL for fast qsort particles reordering
#if 0
  float  Rho, U, Temp, Ne;
#endif
} t_particle_data;

typedef struct particle_data_lite 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;
  int    Id;     // added by JCL for fast qsort particles reordering
} t_particle_data_lite;

class GadgetIO{
public:

    GadgetIO(const std::string&);

    ~GadgetIO();

    int open(const std::string);
    int close();
    int read(const glnemo::t_indexes_tab *index, const int nsel);
    float * getMass()   const { return mass; };
    float * getPos()    const { return pos; };
    float * getVel()    const { return vel; };
    const float * getTime()   const { return &tframe;}
    const int   * getNtotal() const { return &npartTotal;};
    const int getVersion() const { return version;};
    const glnemo::ComponentRangeVector getCRV() const { return crv;};
private:
  std::string filename,file0;
  std::ifstream in;
  std::ofstream out;
  int multiplefiles;
  bool lonely_file;
  //data
  float * mass, * pos, * vel;
  int * id;
  float tframe;
  float ntotmasses;
  t_io_header_1 header;
  int npartTotal, npart;
  bool isLittleEndian();
  bool swap;
  glnemo::ComponentRangeVector  crv;
  void storeComponents();
  //fortran offset record length
  int frecord_offset;
  //control
  bool is_open, is_read;
  bool status;
  int bytes_counter;
  // method
  void readBlockName();
  std::string block_name;
  int readHeader(const int);
  int ioData(char * ptr,const size_t size_bytes,const  int items,const ioop op);
  bool guessVersion();
  int version;

  // swap bytes
  inline void swapBytes(void * x,const int size) {
    char * p=(char *) x;
    for (int i=0;i<size/2;i++) {
     int t=*(p+i); *(p+i)=*(p+size-i-1); *(p+size-i-1)=t;
    }
  }
  // read fortran record
  inline int readFRecord() {
    int len; in.read((char *) &len,sizeof(int));
    // We SWAP data
    if (swap) { // swapping requested
      swapBytes(&len,sizeof(int));
    }
    assert(in.good());
    return len;
  }
};

}

#endif
