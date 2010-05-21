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

#ifndef SNAPSHOTGADGET_H
#define SNAPSHOTGADGET_H

#include <string>
#include <assert.h>
#include <fstream>
#include <map>
#include "snapshotinterface.h"

namespace uns {

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
  //float  Age;
  //float  Metal;
  int    Type;
  int    Id;     // added by JCL for fast qsort particles reordering
} t_particle_data_lite;


  // ---------------------------------------------------
  // CSnapshotGadgetIn
  // READING class
  // ---------------------------------------------------

  class CSnapshotGadgetIn : public CSnapshotInterfaceIn {
    
  public:
    // READING constrcuctor
    CSnapshotGadgetIn(const std::string, const std::string, const std::string, const bool verb=false);
    ~CSnapshotGadgetIn();
    int nextFrame(const uns::t_indexes_tab * index_tab, const int nsel);
//    int nextFrame();
    ComponentRangeVector * getSnapshotRange();
    int getNbody() { return getNtotal();}
    // virtual function implemented
    bool getData(const std::string,int *n,float **);
    bool getData(const std::string,       float * );
    bool getData(const std::string,int *n,int   **);
    bool getData(const std::string,       int   * );
    bool getData(const std::string, const std::string ,int *,float **);
    bool getData(const std::string, const std::string ,int *,int   **);

   int close();

  private:
   bool first_loc;
   int open(const std::string);

   int read (const uns::t_indexes_tab *index, const int nsel);
   int read2(const uns::t_indexes_tab *index, const int nsel);
   int getVersion() const { return version;};
   const uns::ComponentRangeVector getCRV() const { return crv;};
   int   getNtotal() const { return npartTotal;};
     std::string filename,file0;
  std::ifstream in;

  int multiplefiles;
  bool lonely_file;
  //data
  float * mass, * pos, * vel, * rho, * hsml, * age, * metal, * intenerg, * temp;
  int * id;
  int bits; // to store the bits components
  float tframe;
  float ntotmasses;
  t_io_header_1 header;
  int npartTotal, npart;
  bool isLittleEndian();
  bool swap;
  uns::ComponentRangeVector  crv;
  void storeComponents();
  // member data
  float * getMass()   { return mass; }
  float   getTime()   { return tframe;}
  float * getPos()    { return pos; }
  float * getVel()    { return vel; }
  float * getAge(int & n)   { n=header.npartTotal[4]; return age;}
  float * getMetal(int & n) { n=header.npartTotal[0]+header.npartTotal[4]; return metal;}
  float * getMetalGas(int & n) { n=header.npartTotal[0]; return metal;}
  float * getMetalStars(int & n) { n=header.npartTotal[4]; return metal+n;}
  float * getTemp(int & n) { n=header.npartTotal[0]; return temp;}
  float * getU(int & n) { n=header.npartTotal[0]; return intenerg;}
  float * getRho(int & n) { n=header.npartTotal[0]; return rho;}
  float * getHsml(int & n) { n=header.npartTotal[0]; return hsml;}
  
  //fortran offset record length
  int frecord_offset;
  //control
  bool is_open, is_read;
  bool status;
  int bytes_counter;
  // method
  bool readBlockName();
  std::string block_name;
  int readHeader(const int);

  int readData(char * ptr,const size_t size_bytes,const  int items);
  bool guessVersion();
  int version;
  void unitConversion();
  // skip Block
  inline void skipBlock() {
    int len1 = readFRecord();
    in.seekg(len1,std::ios::cur);
    int len2 = readFRecord();
    assert(in.good() && len1==len2);
    if (block_name == "AGE" || block_name == "Z" ) {
      //std::cerr << "len1 = " << len1 << "\nlen2 = " << len2 << "\n";
    }
  }
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

  }; // end of class CSnapshotGadgetIn  





  //
  //
  //

  // ---------------------------------------------------
  // CSnapshotGadgetout
  // WRITING class
  // ---------------------------------------------------

  class CSnapshotGadgetOut : public CSnapshotInterfaceOut {

  public:
    // WRITING constructor
    CSnapshotGadgetOut(const std::string, const std::string, const bool);
    ~CSnapshotGadgetOut();
    int setHeader(void * );
    int setNbody(const int _n);
    int setData(std::string, float);
    int setData(std::string, const int , float *,const bool _addr=false);
    // array by double keys
    int setData(std::string, std::string, const int , float *,const bool _addr=false);
    int setData(std::string, std::string, const int , int   *,const bool _addr=false);
    
    int setData(std::string, const int , int *,const bool _addr=false);
    int setData(std::string, const int , 
		float *, float *, float *, const bool _addr=false);
 

    int save();
    std::vector<double> moveToCom();    

  private:
    //data
    float * mass[6], * pos[6], * vel[6], * rho, * hsml, * age, * metal, * intenerg, * temp;
    int * id[6];
    int ntot_withmasses;
    std::ofstream out;
    // Map to associate the strings with the bool values
    std::map<std::string, bool> ptrIsAlloc[6];
    t_io_header_1 header;
    int bits, npartTotal;
    int bytes_counter;
    int version;
    int setHeader(t_io_header_1 *);
    int writeData(char * ptr,const size_t size_bytes,const int items);
    void saveFile();
    void setupHeader(bool check=false);
    int writeHeader();
    int write();
    bool writeBlockName(std::string);
    inline void writeFRecord(const int len) {
      out.write((char *) &len,sizeof(int));
      assert(out.good());
    }
    // array by keys
    int setMass(std::string, const int _n, float * _mass, const bool addr);
    int setPos (std::string, const int _n, float * _pos , const bool addr);
    int setVel (std::string, const int _n, float * _vel , const bool addr);
    int setId  (std::string, const int _n, int   * _id  , const bool addr);
    
    int setRho (const int _n, float * _rho , const bool addr);
    int setHsml(const int _n, float * _hsml, const bool addr);
    int setU   (const int _n, float * _U   , const bool addr);
  }; // end of class CSnapshotGadgetOut  
} // namespace

#endif
