// ============================================================================
// Copyright Jean-Charles LAMBERT  2008 / 2009                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pale de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

#include "gadgetio.h"
#include <sstream>
#include <assert.h>
#include <algorithm>
#include <string.h>
namespace gadget {

// ============================================================================
// Constructor
  GadgetIO::GadgetIO(const std::string & _f, const bool _verb)
{
  filename = _f;
  status=false;
  is_open=false;
  is_read=false;
  mass   = NULL;
  pos    = NULL;
  vel    = NULL;
  id     = NULL;
  swap   = false;
  tframe = 0.;
  frecord_offset = 4;
  bytes_counter = 0;
  multiplefiles = 0;
  lonely_file   = true;
  verbose = _verb;
  ntotmasses = 0.0;
}

// ============================================================================
// Destuctor                                                                   
GadgetIO::~GadgetIO()
{
  if (mass)   delete [] mass;
  if (pos)    delete [] pos;
  if (vel)    delete [] vel;
  if (id)     delete [] id;
}
// ============================================================================
// open() :                                                                    
// open file and return :                                                      
// 0 : success                                                                 
// 1 : unable to open                                                          
// 2 : not a GADGET file                                                       
int GadgetIO::open(const std::string myfile)
{
  int fail=0;
  in.clear();
  in.open(myfile.c_str(),std::ios::in | std::ios::binary);
  if ( ! in.is_open()) {
   in.close();
   in.clear(); // mandatory under win32 !
   // try to open a multiple gadget file
   file0 = myfile+".0";
   //std::cerr << "In gadget open file0 : " << file0 << "\n";
   in.open(file0.c_str(),std::ios::in | std::ios::binary);
   if (in.is_open()) {
//     assert(in.good());
     lonely_file=false;
     std::cerr << "It's a multiple gadget file.\n";
   }
  }
  if ( ! in.is_open()) {
    fail = 1;               // unable to open
    std::cerr << "In gadget, failed to open...\n";
  }
  else {
    is_open=true;
    if (guessVersion()) {
      fail = readHeader(0);    // try to read header
      if (!fail) status=true; // valid header
      else close();           // not valid header
    }
    else {
      fail=1;
      close();             // not a valid gadgetfile
    }
  }
  return fail;
}
// ============================================================================
int GadgetIO::close()
{
  if (is_open) in.close();
  is_open = false;
  return 1;
}
// ============================================================================
int GadgetIO::read(const glnemo::t_indexes_tab *index, const int nsel)
{
  if (! is_read ) {
    if (! pos)  pos  = new float[3*nsel];
    if (! vel)  vel  = new float[3*nsel];
    if (! mass) mass = new float[nsel];
    if (! id)   id   = new int[nsel];
    is_read=true;
    // allocate memory
    assert(nsel<=npartTotal);
    t_particle_data_lite * P= new t_particle_data_lite[npartTotal];
  
    int pc_new; // particles index
    // loop on all the files
    for (int i=0, pc=0; i<header.num_files || (i==0 && header.num_files==0);i++,pc=pc_new) {
      std::string infile;
      if (header.num_files > 0 ) { // more than one file
        std::ostringstream stm;
        stm << "." << i;               // add ".XX" extension
        infile = filename + stm.str(); // new filename
        if (i>0) {
          close();                 // close previous file
          int fail=open(infile);   // open new file,read header
          assert(!fail);           // fail is true abort
        }
      }
      else infile=filename;        // lonely file
      int len1=0,len2=0;
      bool stop=false;
      std::string next_block_name;
      if (version==1) next_block_name="POS";

      while (readBlockName() && !stop) {
	if (version==1) block_name=next_block_name;
	bool ok=false;
	if (block_name=="POS") { // Postions block
	  ok=true;
	  bytes_counter=0;
	  len1 = readFRecord();
	  pc_new=pc;
	  for(int k=0;k<6;k++)
	    for(int n=0;n<header.npart[k];n++)
	      ioData((char *) &P[pc_new++].Pos[0], sizeof(float), 3, READ);
	  len2 = readFRecord();
	  assert(in.good() && len1==len2 && len1==bytes_counter);
	  if (version==1) next_block_name="VEL";
	}
    
	if (block_name=="VEL") { // Velocities block
	  ok=true;
	  bytes_counter=0;
	  len1 = readFRecord();
	  pc_new=pc;
	  for(int k=0;k<6;k++)
	    for(int n=0;n<header.npart[k];n++)
	      ioData((char *) &P[pc_new++].Vel[0], sizeof(float), 3, READ);
	  len2 = readFRecord();
	  assert(in.good() && len1==len2 && len1==bytes_counter);
	  if (version==1) next_block_name="ID";
	}
  
	if (block_name=="ID") { // IDs block
	  ok=true;
	  bytes_counter=0;
	  len1 = readFRecord();
	  pc_new=pc;
	  for(int k=0;k<6;k++)
	    for(int n=0;n<header.npart[k];n++)
	      ioData((char *) &P[pc_new++].Id, sizeof(int), 1, READ);
	  len2 = readFRecord();
	  assert(in.good() && len1==len2 && len1==bytes_counter);
	  if (version==1) next_block_name="MASS";
	}

	if (block_name=="MASS") { // MASS block
	  ok=true;
	  bytes_counter=0;
	  if (ntotmasses > 0.)
	    len1 = readFRecord();
	  pc_new=pc;
	  for(int k=0;k<6;k++)
	    for(int n=0;n<header.npart[k];n++) {
	      if (header.mass[k] == 0 ) {  // variable mass
		ioData((char *) &P[pc_new++].Mass, sizeof(float), 1, READ);
	      }
	      else {                       // constant mass
		P[pc_new++].Mass = header.mass[k];
	      }
	    }
	  if (ntotmasses > 0.) {
	    len2 = readFRecord();
	    assert(in.good() && len1==len2 && len1==bytes_counter);
	  }
	}

	if (!ok) {
	  skipBlock();
	}
      }
    } // end of loop on numfiles

    // set masses in case of missing mass block name
    if (ntotmasses == 0) {
      if (verbose) {
	std::cerr << "Adding constant masses...\n";
      }
      pc_new=0;
      for(int k=0;k<6;k++)
	for(int n=0;n<header.npartTotal[k];n++) {
	  P[pc_new++].Mass = header.mass[k];
	}
    
    }
    // sort particles according to their indexes
    //qsort(P,npartTotal,sizeof(t_particle_data_lite),gadget::compare);
    
    // copy data to output vectors
    int ic=0;
    for (int i=0; i<npartTotal; i++) {
      int ivalid=index[i].i;
      if (ivalid != -1) {
        assert(ic<=nsel);
        memcpy(pos+(ic*3),&P[ivalid].Pos[0],sizeof(float)*3);
	memcpy(vel+(ic*3),&P[ivalid].Vel[0],sizeof(float)*3);
	memcpy(mass+ic   ,&P[ivalid].Mass  ,sizeof(float)  );
        ic++;
      }
    }
    // garbage collecting
    delete [] P;
  }
  return 1;
}
// ============================================================================
// compare 2 elements                                                          
int compare( const void * a, const void * b )
{
  t_particle_data_lite 
    * pa = ( t_particle_data_lite * ) a, 
    * pb = ( t_particle_data_lite * ) b;
  return (pa->Id - pb->Id);
}

// ============================================================================
// readHeader():                                                               
// read gadget header structure                                                
// return values : 0 success, >0 not a Gadget header data                      
int GadgetIO::readHeader(const int id)
{
  int len1,len2;

  // Header block
  readBlockName();
  bytes_counter=0;
  len1 = readFRecord();
  ioData((char *)  header.npart         , sizeof(int   ),  6, READ);
  ioData((char *)  header.mass          , sizeof(double),  6, READ);
  ioData((char *) &header.time          , sizeof(double),  1, READ);
  ioData((char *) &header.redshift      , sizeof(double),  1, READ);
  ioData((char *) &header.flag_sfr      , sizeof(int   ),  1, READ);
  ioData((char *) &header.flag_feedback , sizeof(int   ),  1, READ);
  ioData((char *)  header.npartTotal    , sizeof(int   ),  6, READ);
  ioData((char *) &header.flag_cooling  , sizeof(int   ),  1, READ);
  ioData((char *) &header.num_files     , sizeof(int   ),  1, READ);
  ioData((char *) &header.BoxSize       , sizeof(double),  1, READ);
  ioData((char *) &header.Omega0        , sizeof(double),  1, READ);
  ioData((char *) &header.OmegaLambda   , sizeof(double),  1, READ);
  ioData((char *) &header.HubbleParam   , sizeof(double),  1, READ);
  ioData((char *)  header.fill          , sizeof(char  ), 96, READ);
  len2 = readFRecord();
  //std::cerr << "len1="<<len1<< " len2="<<len2<< " bytes_counter="<<bytes_counter<<"\n";
  if (in.bad() || len1!=len2 || len1!=bytes_counter)
    return 2;
  if (id==0) {                         // first time
    tframe = header.time;
    npartTotal = 0;
    for(int k=0; k<6; k++)  {
      npartTotal += header.npartTotal[k];  // count particles
    }
    for(int k=0; k<6; k++) {
      if (header.mass[k] == 0 ) {
	ntotmasses += header.npartTotal[k];
      }
    }
    storeComponents();
  }
  return 0;
}
// ============================================================================
// readBlockName : read Gadget2 file format block                              
bool GadgetIO::readBlockName()
{
  bool status=true;
  if (version == 2 ) { // gadget2 file format
    int dummy,dummy1;
    char name[9];
    ioData((char *) &dummy, sizeof(int),  1, READ); // read
    ioData((char *) name  , sizeof(char), 8, READ); // read
    ioData((char *) &dummy1, sizeof(int),  1, READ); // read
    int i=0; while (isupper(name[i])&& i<4) i++;
    name[i]='\0';
    block_name=name;
    block_name=name;
    status = in.good();
    if (verbose) {
      std::cerr << "Block Name : <" << block_name << ">\n";
    }
  }
  return status;
}
// ============================================================================
// guessVersion()                                                              
// detect Gadget2 file format version ( 1 or 2 )                               
// return true if successfully detected otherwise false                        
bool GadgetIO::guessVersion()
{
  bool status=true;
  // try to read 32bits length fortran record
  int dummy;
  swap = false; // no swapping
  ioData((char *) &dummy, sizeof(int), 1, READ); // read
  //std::cerr << "GadgetIO::guessVersion dummy = "<<dummy<<"\n";
  if( dummy != 256  && dummy != 8  ) {           // unknow number
    swap = true;                                 // need bytes swapping
    swapBytes((int *) &dummy, sizeof(int));      // swap
    //std::cerr << "GadgetIO::guessVersion dummy swapped= "<<dummy<<"\n";
    if( dummy != 256  && dummy != 8  )           // unknow swapped number
     status = false;                             // not a gadget file    
  }
  if (status) {
    if (dummy==256) version=1; // gadget1
    else            version=2; // gadget2
    in.seekg(0,std::ios::beg); // rewind
    //std::cerr << "gadget Version : " <<version << "\n"; 
  }

  return status;
}
// ============================================================================
// ioData:                                                                     
// perform IO (Read,Write) operation on Data                                   
int GadgetIO::ioData(char * ptr,const size_t size_bytes,const int items,const ioop op)
{
  char * pp;
  bytes_counter += (size_bytes*items);
  switch (op) {

  case READ :
    // get data from file
    in.read(ptr,size_bytes*items);
    //assert(in.good());
    if (! in.good()) return 0;
    // We SWAP data
    if (swap && (size_bytes != CHAR)) { // swapping requested
      for (int i=0; i<items; i++) {
	swapBytes(ptr,size_bytes);
	ptr += size_bytes;
      }
    }
    break;

  case WRITE :
    // We must SWAP data first
    if (swap && (size_bytes != CHAR)) { // swapping requested
      pp = ptr;
      for (int i=0; i<items; i++) {
	swapBytes(ptr,size_bytes);
	ptr += size_bytes;
      }
      ptr = pp;
    }
    // Save Data to file
    out.write(ptr,size_bytes*items);
    assert(out.good());
    
    // We have to unswap the data now
    if (swap && (size_bytes != CHAR)) { // swapping requested
      for (int i=0; i<items; i++) {
	swapBytes(ptr,size_bytes);
	ptr += size_bytes;
      }
    }
    break;
  } // end of switch (op) ....

  return 1;
}
// ============================================================================
// isLittleEndian()                                                            
// Test endianess: return true if little endian architecture, false if big endian
bool GadgetIO::isLittleEndian()
{
  bool status=false;
  int one=1;
  char * p = (char *)&one;
  if (p[0]==1 ) status=true;
  return status;
}
// ============================================================================
// storeComponents:                                                            
void GadgetIO::storeComponents()
{
  glnemo::ComponentRange cr;
  // all
  cr.setData(0,npartTotal-1);
  cr.setType("all");
  crv.clear();
  crv.push_back(cr);
  // components
  const char * comp [] = { "gas", "halo", "disk", "bulge", "stars", "bndry"};
  for(int k=0,start=0; k<6; k++)  {
    if (header.npartTotal[k]) {
      cr.setData(start,start+header.npartTotal[k]-1,comp[k]);
      crv.push_back(cr);
      start+=header.npartTotal[k];
    }
  }
}

}
