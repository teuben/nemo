// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include "ftmio.h"
#include <assert.h>

namespace ftm {

int FtmIO::frecord_offset=4;
// ============================================================================
// constructor                                                                 
FtmIO::FtmIO(const std::string _f)
{
  filename = _f;
  status=false;
  is_open=false;
  is_read=false;
  mass   = NULL;
  pos    = NULL;
  vel    = NULL;
  classe = NULL;
  index1 = NULL;
  comp   = NULL;
}

// ============================================================================
// Destuctor                                                                   
FtmIO::~FtmIO()
{
  if (mass)   delete [] mass;
  if (pos)    delete [] pos;
  if (vel)    delete [] vel;
  if (classe) delete [] classe;
  if (index1) delete [] index1;
  if (comp)   delete    comp;
  mass = NULL;
  pos  = NULL;
  vel  = NULL;
}
// ============================================================================
// open() :                                                                    
// open file and return :                                                      
// 0 : success                                                                 
// 1 : unable to open                                                          
// 2 : not an FTM file                                                         
int FtmIO::open()
{
  int ret;
  in.open(filename.c_str(),std::ios::in | std::ios::binary);
  if ( ! in.is_open())
    ret = 1;               // unable to open
  else {
    is_open=true;
    ret = readHeader();    // try to read header
    if (!ret) status=true; // valid header      
    else close();          // not valid header  
  }
  return ret;
}
// ============================================================================
int FtmIO::close()
{
  if (is_open) in.close();
  is_open = false;
  return 1;
}
// ============================================================================
int FtmIO::readFRecord()
{
  int len;
  in.read((char *) &len,sizeof(int));
  return len;
}
// ============================================================================
// readHeader()                                                                
// read FTM reader and return :                                                
// 0 : success                                                                 
// 2 : not an FTM file                                                         
int FtmIO::readHeader()
{
  int len1,len2;
  // read record
  len1=readFRecord();
  smartRead((char *) &dmpindx  , sizeof(int   )    );
  smartRead((char *) &dmp_date , sizeof(char  )*24 );
  smartRead((char *) &ndim     , sizeof(int   )    );
  smartRead((char *) &eqnindx  , sizeof(int   )    );
  smartRead((char *) &version  , sizeof(double)    );
  smartRead((char *) &gamma    , sizeof(double)    );
  smartRead((char *) &poly     , sizeof(double)    );
  smartRead((char *) &lsfm     , sizeof(int   )    );
  len2=readFRecord();
  // check record length
  int check_len=sizeof(int)*4+sizeof(double)*3+sizeof(char)*24;
  if (in.bad() || len1!=len2 || len1!=check_len)
          return 2;
  // read record
  len1=readFRecord();
  smartRead((char *) &tframe   , sizeof(double)    );
  smartRead((char *) &n1       , sizeof(int   )    );
  smartRead((char *) &n2       , sizeof(int   )    );
  len2=readFRecord();
  // check record length
  check_len=sizeof(int)*2+sizeof(double);
  if (in.bad() || len1!=len2 || len1!=check_len)
         return 2;
  dmp_date[24]='\0';
  n0 = n1+n2;
  // load classes
  if (classe) delete [] classe;
  if (index1) delete [] index1;
  classe =  new int[n0];
  index1 =  new int[n0];
  len1=readFRecord();
  smartRead((char *) classe   , sizeof(int)*n0);
  len2=readFRecord();
  assert(len1 == len2);
  // manage classe
  classComponents();
  std::cerr << "\ndmpidx = " << dmpindx << "\n";
  return 0;
}
// ============================================================================
// smartRead()
// 
int FtmIO::smartRead(char * p, unsigned int len)
{
  in.read(p,len);
  return in.good();
}
// ============================================================================
int FtmIO::skipBytes(unsigned int bytes)
{
  in.seekg(bytes+(frecord_offset*2),std::ios::cur); // fortran record padding offset
  return in.good();
}
// ============================================================================
int FtmIO::read(glnemo::ParticlesData * part_data, const int *index, const int nsel, bool load_vel)
{
  int ret=0;
  int len1,len2;
  float * tmp3n;
  //if ( ! is_read) {
  if (is_open) {
    if (status) {
      tmp3n  =  new float[n0];
    }
    if ( *part_data->nbody < nsel) {
      *part_data->nbody = nsel;
      if (part_data->pos) delete [] part_data->pos;
      part_data->pos = new float[nsel*3];
      if (load_vel) {
        if (part_data->vel) delete [] part_data->vel;
        part_data->vel = new float[nsel*3];
      }
      if (dmpindx>2 && n2>0) { // gas density and hsml
        if (part_data->rho) delete part_data->rho;
        part_data->rho = new glnemo::PhysicalData(glnemo::PhysicalData::rho,nsel);
        for (int i=0; i<nsel; i++) part_data->rho->data[i]=-1.;
        if (part_data->rneib) delete part_data->rneib;
        part_data->rneib = new glnemo::PhysicalData(glnemo::PhysicalData::neib,nsel);
      }
    }
    *part_data->nbody = nsel;
    *part_data->timu  = tframe;
    // load masses
    len1=readFRecord();
    smartRead((char *) tmp3n     , sizeof(float)*n0);
    len2=readFRecord();
    assert(len1 == len2);
    // load pos
    for (int i=0; i<3; i++) {
      len1=readFRecord();
      smartRead((char *) tmp3n  , sizeof(float)*n0);
      swapArrayIndex3D(tmp3n,part_data->pos,n0, i, index, nsel);
      len2=readFRecord();
      assert(len1 == len2);
    }
    // load vel
    for (int i=0; i<3; i++) {
      len1=readFRecord();
      smartRead((char *) tmp3n  , sizeof(float)*n0);
      if (load_vel)
        swapArrayIndex3D(tmp3n,part_data->vel,n0, i, index,nsel);
      len2=readFRecord();
      assert(len1 == len2);
    }
    if (dmpindx>1) {
      skipBlock(); // eps
      skipBlock(); // phi
      for (int i=0; i<3; i++) {
        skipBlock(); // force
      }
      if (lsfm) {
        skipBlock(); // state
        skipBlock(); // state
        skipBlock(); // sfmt
        skipBlock(); // stmz
      }
    }
    if (dmpindx>2 && n2>0) {
      skipBlock(); // u

      // load gas density
      len1=readFRecord();
      smartRead((char *) tmp3n     , sizeof(float)*n2);
      len2=readFRecord();
      assert(len1 == len2);
      for (int i=0, k=0; i<n0; i++) {
        if (index[i]!=-1) {
          if (i>=n1) {
            part_data->rho->data[k] = tmp3n[k];
            k++;
          }
          else {
            part_data->rho->data[k] = 1;
          }
        }
      }

      skipBlock(); // divv
      skipBlock(); // udot1
      skipBlock(); // udot2
      skipBlock(); // udot3
      skipBlock(); // udot4
      
      // load gas hsml
      len1=readFRecord();
      smartRead((char *) tmp3n     , sizeof(float)*n2);
      len2=readFRecord();
      assert(len1 == len2);
      for (int i=0, k=0; i<n0; i++) {
        if (index[i]!=-1) {
          if (i>=n1) {
            part_data->rneib->data[k] = tmp3n[k];
            k++;
          }
          else {
            part_data->rneib->data[k] = 1;
          }
        }
      }

    }
    // garbage collecting
    delete [] tmp3n;
    is_read=true;
    ret=1;
  }
  return ret;

}
// ============================================================================
int FtmIO::read()
{
  int ret=0;
  int len1,len2;
  float * tmp3n;
  //if ( ! is_read) {
  if ( is_open) {
    if (status) {
      pos    =  new float[3*n0];
      vel    =  new float[3*n0];
      mass   =  new float[n0];
      tmp3n  =  new float[n0];
    }

    // load masses
    len1=readFRecord();
    smartRead((char *) mass     , sizeof(float)*n0);
    len2=readFRecord();
    assert(len1 == len2);
    // load pos
    for (int i=0; i<3; i++) {
      len1=readFRecord();
      smartRead((char *) tmp3n  , sizeof(float)*n0);
      swapArrayIndex3D(tmp3n,pos,n0, i);
      len2=readFRecord();
      assert(len1 == len2);
    }
    // load vel
    for (int i=0; i<3; i++) {
      len1=readFRecord();
      smartRead((char *) tmp3n  , sizeof(float)*n0);
      swapArrayIndex3D(tmp3n,vel,n0, i);
      len2=readFRecord();
      assert(len1 == len2);
    }
    // garbage collecting
    delete [] tmp3n;
    is_read=true;
    ret=1;
  }
  return ret;
}
// ============================================================================
// swapArrayIndex3D                                                            
// reverse column range order                                                  
void FtmIO::swapArrayIndex3D(float * ain, float * aout, int n, int dim,
                             const int * tab, const int nsel)
{
  int inew=0;
  for (int i=0; i<n; i++) {
    if (tab[i] != -1 ) {           // particles must be selected
      aout[3*inew+dim] = ain[i] ;  // swap in[3,n] to out[n,3]  
      inew++;                      // new particle              
      assert(inew<=nsel);
    }
  }
}
// ============================================================================
// swapArrayIndex3D                                                            
// reverse column range order                                                  
void FtmIO::swapArrayIndex3D(float * ain, float * aout, int n, int dim)
{
  for (int i=0; i<n; i++) {
    aout[3*i+dim] = ain[i] ;    // swap in[3,n] to out[n,3]
  }
}
// ============================================================================
// swapArrayIndex3D                                                            
// reverse column range order                                                  
void FtmIO::swapArrayIndex3D(float * ain, float * aout, int n)
{
  for (int i=0; i<n; i++) {
    for (int j=0; j<3; j++) {
      aout[3*i+j] = ain[j*n+i] ; // swap in[3,n] to out[n,3]
    }
  }
}
// ============================================================================
// classComponents                                                             
// class particles in Halo, Disk, Hole and Gas components                      
// class - classification:                                                     
//                        classe < -100,000   Stellar halo particle            
//            -100,000 <= classe < 0          Stellar disk particle            
//                        classe = 0          Central BH                       
//                    0 < classe              SPH particle                     
// Indexes :                                                                   
//    Gas:1, Hole:2, Disk:3, Halo:4                                            
void FtmIO::classComponents()
{
  if (comp) delete comp;
  comp = new FtmComponent(classe,index1,n0,n1,n2);
}
// ============================================================================
//                  Class FtmComponent                                         
// ============================================================================

// ============================================================================
// constructor
FtmComponent::FtmComponent(int * classe, int * index1,int n0, int n1, int n2)
{
  ihole = 0;
  idisk = 1;
  igas  = 1;
  for (int i=0; i<n0; i++) {
    if (classe[i] > 0) index1[i] = 1;          // Gas 
    else
      if (classe[i] == 0) {                    // Hole
        index1[i] = 2;
        ihole = i;
      }
      else
        if (classe[i] > -1000000) {            // Disk
          index1[i] = 3;
          igas = i + 1;
        }
        else {                                 // Halo
          index1[i] = 4;
          idisk = i + 1;
          igas = i + 1;
        }
  }
  if (idisk > 1 )    halo.setData(0,idisk-2,"halo");          // halo exist
  if (idisk != igas) disk.setData(idisk-1,igas-2,"disk");     // disk exist
  if (n2>0)          gas.setData(igas-1, igas-1+n2-1, "gas"); // gas exist
  if (n1) {;} // remove compiler warning
}
} // end of namespace ftm
