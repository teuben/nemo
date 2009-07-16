// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2009                                  
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
#include "gadgetio.h"
#include <sstream>
#include <assert.h>
#include <algorithm>
#include <string.h>
#include <cmath>

namespace gadget {

// ============================================================================
// Constructor
GadgetIO::GadgetIO(const std::string & _f)
{
  filename = _f;
  status=false;
  is_open=false;
  is_read=false;
  mass   = NULL;
  pos    = NULL;
  vel    = NULL;
  intenerg=NULL;
  swap   = false;
  tframe = 0.;
  frecord_offset = 4;
  bytes_counter = 0;
  multiplefiles = 0;
  lonely_file   = true;
  ntot_withmasses = 0;
}

// ============================================================================
// Destuctor                                                                   
GadgetIO::~GadgetIO()
{
  if (mass)     delete [] mass;
  if (pos)      delete [] pos;
  if (vel)      delete [] vel;
  if (intenerg) delete [] intenerg;
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
  std::cerr << "In gadget open file : " << myfile << "\n";
  in.open(myfile.c_str(),std::ios::in | std::ios::binary);
  if ( ! in.is_open()) {
   in.close();
   in.clear(); // mandatory under win32 !
   // try to open a multiple gadget file
   file0 = myfile+".0";
   std::cerr << "In gadget open file0 : " << file0 << "\n";
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
int GadgetIO::read2(float * pos, float * vel, float * rho, float * rneib, float * temp,const int *index, const int nsel,   const bool load_vel)
{
  bool is_temp=true;
  if (! is_read ) {
    use_gas = false;
    s_gas=1e9;  // starting gas index according to the user
    e_gas=-1;   // ending   gas index according to the user
    // allocate memory
    assert(nsel<=npartTotal);
    is_read=true;
    // allocate memory
    assert(nsel<=npartTotal);
    //t_particle_data_lite * P= new t_particle_data_lite[npartTotal];
    if (! intenerg && header.npartTotal[0]>0) intenerg = new float[header.npartTotal[0]];
    
    int pc_new=0; // particles index
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

      // Read the whole file till a valid block exist
      // or stop = true in case of Gadget1 file      
      while (readBlockName() && !stop) {
        if (version==1) block_name=next_block_name;
        bool ok=false;
        if (block_name=="POS") { // Postions block
          ok=true;
          bytes_counter=0;
          len1 = readFRecord();
          pc_new=pc;
          int ic=0;
          for(int k=0;k<6;k++) {
            for(int n=0;n<header.npart[k];n++) {
              int idx=index[pc_new++];
              if (idx != -1) {
                ioData((char *) &pos[3*ic++], sizeof(float), 3,READ);
              } else {
                //skipData(sizeof(float)*3);
                float tmp3[3];
                ioData((char *) tmp3, sizeof(float), 3,READ);
              }
            }
          }
          len2 = readFRecord();
          assert(in.good() && len1==len2 && len1==bytes_counter);
          if (version==1) next_block_name="VEL";
        }
    
        if (block_name=="VEL") { // Velocities block
          ok=true;
          if (load_vel || version==1) {
            bytes_counter=0;
            len1 = readFRecord();
            pc_new=pc;
            int ic=0;
            for(int k=0;k<6;k++)
              for(int n=0;n<header.npart[k];n++){
                int idx=index[pc_new++];
                if (idx != -1 && load_vel) {
                    ioData((char *) &vel[3*ic++], sizeof(float), 3,READ);                  
                } else {
                  //skipData(sizeof(float)*3);
                  float tmp3[3];
                  ioData((char *) tmp3, sizeof(float), 3,READ);
                }
              }
            len2 = readFRecord();
            assert(in.good() && len1==len2 && len1==bytes_counter);
            //copy3DtoArray(P,vel,index,npartTotal);
          } 
          else {
            skipBlock();
          }
          if (version==1) next_block_name="ID";
        }
  
        if (block_name=="ID") { // IDs block
          ok=true;
          skipBlock();
          if (version==1) next_block_name="MASS";
        }

        if (block_name=="MASS") { // MASS block
          ok=true;
          if (ntot_withmasses>0)
            skipBlock();  
          if (version==1) stop=true; // we stop reading for gadget1
        }

        if (block_name=="U") { // U block (Internal energy)
          assert(header.npart[0]>0);
          ok=true;
          bytes_counter=0;
          len1 = readFRecord();
          ioData((char *) &intenerg[pc], sizeof(float), header.npart[0],READ);
          len2 = readFRecord();
          assert(in.good() && len1==len2 && len1==bytes_counter);
        }

        if (block_name=="NE") { // 
          assert(header.npart[0]>0);
          is_temp = true;
          ok=true;
          bytes_counter=0;
          len1 = readFRecord();
	  // getting NE for gas only
          pc_new=pc;
          int ic=0;
          for(int k=0;k<6;k++)
            for(int n=0;n<header.npart[k];n++){
            int idx=index[pc_new++];
            if (idx != -1) {
              if (k == 0 ) { // gaz component only
                use_gas=true;
                s_gas = std::min(s_gas,ic);
                e_gas = std::max(e_gas,ic+1);
                ioData((char *) &temp[ic++], sizeof(float), 1,READ);
              }
            } else {
              if (k == 0) {
                float tmp;
                ioData((char *) &tmp, sizeof(float), 1,READ);
              }
            }
          }
          //ioData((char *) &temp[pc], sizeof(float),header.npart[0],READ);
          len2 = readFRecord();
          assert(in.good() && len1==len2 && len1==bytes_counter);
        }

        if (block_name=="RHO") { // RHO block (density)
          assert(header.npart[0]>0);
          ok=true;
          bytes_counter=0;
          len1 = readFRecord();
	  // getting RHO for gas only
          pc_new=pc;
          int ic=0;
          for(int k=0;k<6;k++)
            for(int n=0;n<header.npart[k];n++){
            int idx=index[pc_new++];
            if (idx != -1) {
              if (k == 0 ) { // gaz component only
                use_gas=true;
                s_gas = std::min(s_gas,ic);
                e_gas = std::max(e_gas,ic+1);
                ioData((char *) &rho[ic++], sizeof(float), 1,READ);
              }
            } else {
              if (k == 0) {
                float tmp;
                ioData((char *) &tmp, sizeof(float), 1,READ);
              }
            }
          }

          //ioData((char *) &rho[pc], sizeof(float),header.npart[0],READ);
          len2 = readFRecord();
          assert(in.good() && len1==len2 && len1==bytes_counter);
        }

        if (block_name=="HSML") { // HSML block (neighbours size)
          assert(header.npart[0]>0);
          ok=true;
          bytes_counter=0;
          len1 = readFRecord();
	  // getting HSML for gas only
          pc_new=pc;
          int ic=0;
          for(int k=0;k<6;k++)
            for(int n=0;n<header.npart[k];n++){
              int idx=index[pc_new++];
              if (idx != -1) {
                assert(ic<nsel);
                if (k == 0 ) { // gaz component only
                  ioData((char *) &rneib[ic++], sizeof(float), 1,READ);
                }
              } else {
                if (k == 0) {
                  float tmp;
                  ioData((char *) &tmp, sizeof(float), 1,READ);
                }
              }
          }
          //ioData((char *) &rho[pc], sizeof(float),header.npart[0],READ);
          len2 = readFRecord();
          assert(in.good() && len1==len2 && len1==bytes_counter);
        }

        if (block_name=="Z") { // Z block (Metalicity)
          ok=true;
          skipBlock();
        }

        if (block_name=="AGE") { // AGE block
          ok=true;
          skipBlock();
        }

        if (!ok) {
          skipBlock();
        }
      }
    } // end of loop on numfiles

    std::cerr << "Use gas = " << use_gas << " start=" << s_gas << " end=" << e_gas << "\n";
    if (version == 2 && is_temp && use_gas && ! header.flag_cooling) {
      for(int n=0;n<nsel;n++) {
        temp[n] = 1.0;
      }
    }
    // sort particles according to their indexes
    //qsort(P,npartTotal,sizeof(t_particle_data_lite),gadget::compare);
    
    // convert to temperature units only if user has selected gas particles
    if (use_gas && is_temp && header.npartTotal[0] > 0) {
      rhop  = rho;
      tempp = temp;
      unitConversion();
    }
    if (! is_temp) {
      delete [] temp;
      temp=NULL;
    }
    // garbage collecting
  }
  return 1;
}
// ============================================================================
int GadgetIO::read(float * pos, float * vel, float * rho, float * rneib, const int *index, const int nsel,
                   const bool load_vel)
{
  if (! is_read ) {
    is_read=true;
    // allocate memory
    assert(nsel<=npartTotal);
    //t_particle_data_lite * P= new t_particle_data_lite[npartTotal];
    float * P = new float[npartTotal*3];
  
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
      int len1,len2;
      
      // Postions block
      readBlockName();
      bytes_counter=0;
      len1 = readFRecord();
      pc_new=pc;
      for(int k=0;k<6;k++)
        for(int n=0;n<header.npart[k];n++) {
          ioData((char *) &P[3*pc_new++], sizeof(float), 3,READ);
        }
      len2 = readFRecord();
      assert(in.good() && len1==len2 && len1==bytes_counter);
      //copy3DtoArray(P,pos,index,npartTotal);
      copy3DtoArray(P,pos,index,npartTotal);
      
      // Velocities block
      readBlockName();
      if (load_vel) {
        bytes_counter=0;
        len1 = readFRecord();
        pc_new=pc;
        for(int k=0;k<6;k++)
          for(int n=0;n<header.npart[k];n++)
            ioData((char *) &P[3*pc_new++], sizeof(float), 3,READ);
        len2 = readFRecord();
        assert(in.good() && len1==len2 && len1==bytes_counter);
        copy3DtoArray(P,vel,index,npartTotal);
      }
      else skipBlock();
      
      // IDs block
      readBlockName();
      bytes_counter=0;
      len1 = readFRecord();
      pc_new=pc;
      for(int k=0;k<6;k++)
        for(int n=0;n<header.npart[k];n++)
          ioData((char *) &P[pc_new++], sizeof(int), 1,READ);
      len2 = readFRecord();
      assert(in.good() && len1==len2 && len1==bytes_counter);
  
      // MASS block
      
      bytes_counter=0;
      if (ntot_withmasses>0) {
        readBlockName();
        skipBlock();
      }

      // U block
      if (header.npart[0] > 0) { // proceed only if there are gas particles
	readBlockName();
        skipBlock();
      }
      
      // Rho block
      if (header.npart[0] > 0 && version==2) { // proceed only if there are gas particles
	if (version == 2) {
	  while (block_name != "RHO") {         // find RHO blockname
	    readBlockName();
	    if (block_name != "RHO") {
              skipBlock();
	    }
	  }
	}
	if (block_name == "RHO") { //||version==1) {
	  bytes_counter=0;
	  len1 = readFRecord();
	  pc_new=pc;
	  for(int k=0;k<6;k++)
	    for(int n=0;n<header.npart[k];n++) {
	      if (k==0) { // read only gas particles
	       ioData((char *) &P[pc_new++], sizeof(float), 1,READ);
	      } else {
		P[pc_new++] = -1;
	      }
	    }
	  len2 = readFRecord();
	  assert(in.good() && len1==len2 && len1==bytes_counter);
          copy1DtoArray(P,rho,index,npartTotal);
	}
      }
      if (version == 1 && header.flag_sfr) {
          skipBlock();
          skipBlock();       
      }
      // HMSL block
      if (header.npart[0] > 0 && version==2) { // proceed only if there are gas particles
	if (version == 2) {
	  while (block_name != "HSML") {         // find RHO blockname
	    readBlockName();
	    if (block_name != "HSML") {
              skipBlock();
	    }
	  }
	}
	if (block_name == "HSML"||(version==1 && header.flag_cooling)) {
	  bytes_counter=0;
	  len1 = readFRecord();
	  pc_new=pc;
	  for(int k=0;k<6;k++)
	    for(int n=0;n<header.npart[k];n++) {
	      if (k==0 ) { // read only gas particles
	        ioData((char *) &P[pc_new++], sizeof(float), 1,READ);
	      }
              else {
	        P[pc_new++] = -1;
	      }
	    }
	    len2 = readFRecord();
	    assert(in.good() && len1==len2 && len1==bytes_counter);
            copy1DtoArray(P,rneib,index,npartTotal);
	} else {
	    if ((version==1 && !header.flag_cooling)) {
              pc_new=pc;
	      for(int k=0;k<6;k++)
		for(int n=0;n<header.npart[k];n++) {
		  if (k==0 ) { // read only gas particles
		      P[pc_new++] = 1;
		  }
                  else {
		      P[pc_new++] = -1;
		  }
	        }
              copy1DtoArray(P,rneib,index,npartTotal);
            }
	  }
      }
    } // end of loop on numfiles
  
    // sort particles according to their indexes
    //qsort(P,npartTotal,sizeof(t_particle_data_lite),gadget::compare);
    
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
  std::cerr << "len1="<<len1<< " len2="<<len2<< " bytes_counter="<<bytes_counter<<"\n";
  if (in.bad() || len1!=len2 || len1!=bytes_counter)
    return 2;
  if (id==0) {                         // first time
    tframe = header.time;
    npartTotal = 0;
    for(int k=0; k<6; k++)  npartTotal += header.npartTotal[k];  // count particles
    for(int k=0; k<6; k++) {
      if(header.mass[k]==0) {
        ntot_withmasses+= header.npartTotal[k]; // count particles with mass on file
      }
    }
    std::cerr << "Nto with masses : " << ntot_withmasses << "\n";
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
    int dummy;
    char name[9];
    ioData((char *) &dummy, sizeof(int),  1, READ); // read
    ioData((char *) name  , sizeof(char), 8, READ); // read
    ioData((char *) &dummy, sizeof(int),  1, READ); // read
    int i=0; while ((isupper(name[i])||isdigit(name[i]))&& i<4) i++;
    name[i]='\0';
    block_name=name;
    status = in.good();
    std::cerr <<">> Blockname ="<<block_name<<"\n";
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
  std::cerr << "GadgetIO::guessVersion dummy = "<<dummy<<"\n";
  if( dummy != 256  && dummy != 8  ) {           // unknow number
    swap = true;                                 // need bytes swapping
    swapBytes((int *) &dummy, sizeof(int));      // swap
    std::cerr << "GadgetIO::guessVersion dummy swapped= "<<dummy<<"\n";
    if( dummy != 256  && dummy != 8  )           // unknow swapped number
     status = false;                             // not a gadget file    
  }
  if (status) {
    if (dummy==256) version=1; // gadget1
    else            version=2; // gadget2
    in.seekg(0,std::ios::beg); // rewind
    std::cerr << "gadget Version : " <<version << "\n"; 
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
// ============================================================================
// unitConversion()                                                            
void GadgetIO::unitConversion()
{
  double GRAVITY, BOLTZMANN, PROTONMASS;
  double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
  double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;
  double G, Xh, HubbleParam;

  double MeanWeight, u, gamma;
  double RhoUniverse_omegabar;
  assert(use_gas);
  /* physical constants in cgs units */
  GRAVITY   = 6.672e-8;
  BOLTZMANN = 1.3806e-16;
  PROTONMASS = 1.6726e-24;

  /* internal unit system of the code */
  UnitLength_in_cm= 3.085678e21;   /*  code length unit in cm/h */
  UnitMass_in_g= 1.989e43;         /*  code mass unit in g/h */
#ifdef GRAPE_UNIT 
  UnitVelocity_in_cm_per_s= 207.9*1.0e5;
#else
  UnitVelocity_in_cm_per_s= 1.0e5;
#endif

  UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitDensity_in_cgs= UnitMass_in_g/ pow(UnitLength_in_cm,3);
  UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ pow(UnitTime_in_s,2);
  UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);

  G=GRAVITY/ pow(UnitLength_in_cm,3) * UnitMass_in_g * pow(UnitTime_in_s,2);

  Xh= 0.76;  /* mass fraction of hydrogen */
  HubbleParam= 0.65; 
  
  RhoUniverse_omegabar=1.9e-29*0.04;

  assert(intenerg != NULL);

  for(int i=s_gas;i<e_gas;i++) {

    MeanWeight= 4.0/(3*Xh+1+4*Xh*tempp[i]) * PROTONMASS;

    /* convert internal energy to cgs units */

    u  = intenerg[i-s_gas] * UnitEnergy_in_cgs/ UnitMass_in_g;

    gamma= 5.0/3;

    /* get temperature in Kelvin */

    tempp[i] = MeanWeight/BOLTZMANN * (gamma-1) * u;
    if (rhop) {
      rhop[i] *= UnitDensity_in_cgs/RhoUniverse_omegabar;
    }
	  
#if 0
#ifndef NOAGE
    if(P[i].Type==4)
    {
    if(P[i].Age!=0)
    P[i].Age=1.0/P[i].Age-1;
    else
    {
	    //cout<<"Error Age"<<endl;
	    //cout<<"N="<<i<<P[i].Age<< endl;
	    //exit(0);
  }

  }
#endif
#endif
  }

}

}
