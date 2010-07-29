// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008-2010                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

/**
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */
#ifndef UNSSNAPSHOTNEMO_H
#define UNSSNAPSHOTNEMO_H
#include "snapshotinterface.h"
#include <map>

extern "C" {
  int io_nemo(const char * , const char *, ...);
#include <stdinc.h>
#include <filestruct.h>
#include <nemo.h>
#include <snapshot/snapshot.h>
};

namespace uns {

  class CSnapshotNemoIn : public CSnapshotInterfaceIn {
  
  public:
    CSnapshotNemoIn(const std::string, const std::string, const std::string,
		  const bool verb=false);
    ~CSnapshotNemoIn();
    int nextFrame(const uns::t_indexes_tab * index_tab, const int nsel);
    int close();
    ComponentRangeVector * getSnapshotRange();
    // virtual function implemented
    bool getData(const std::string,int *n,float **);
    bool getData(const std::string,       float * );
    bool getData(const std::string,int *n,int   **);
    bool getData(const std::string,       int   * );
    bool getData(const std::string, const std::string ,int *,float **);
    bool getData(const std::string, const std::string ,int *,int   **);
   
private:
    int full_nbody;
    int * nemobits , * ionbody;
    float * iotime, *iopos, *iovel, *iomass, *iorho, *ioaux, *ioacc, *iopot;
    float * pos, *vel, *mass, * rho, *acc, *aux, *pot;
    bool first_stream;
    int status_ionemo;
    int last_nbody;
    void checkBits(std::string,const int);
    bool isValidNemo();
    float *  getPos()  { //checkBits("pos",PosBit); 
                         return pos ;}
    float *  getVel()  { //checkBits("vel",VelBit); 
                         return vel ;}
    float *  getMass() { //checkBits("mass",MassBit); 
                         return mass;}
    float *  getRho()  { return rho ;}
    float *  getAux()  { return aux ;}
    float *  getAcc()  { return acc ;}
    float *  getPot()  { return pot ;}
    float    getTime() { return *iotime; }
    int      getNbody(){ return *ionbody;}
};
  
  class CSnapshotNemoOut : public CSnapshotInterfaceOut {

  public:
    // WRITING constructor
    CSnapshotNemoOut(const std::string, const std::string, const bool);
    ~CSnapshotNemoOut();
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
    int close();
  private:
    // Map to associate the strings with the bool values
    std::map<std::string, bool> ptrIsAlloc;
    float * mass, * pos, * vel, * aux, * acc, * pot, * rho;
    float time;
    int * keys;
    int nbody;
    int bits;
    bool is_saved, is_closed;
    // array
    int setArray(const int _n, const int _d, float * src, float ** dest, const char * name, const int tbits, const bool addr);
    int setArray(const int _n, const int _d, int   * src, int   ** dest, const char * name, const int tbits, const bool addr);
};
}
#endif
