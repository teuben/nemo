// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// body.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <body.h>                                  // nbdy::bodies etc          
#include <public/ionl.h>                           // utilities for C++ I/O     
#include <iostream>                                // C++ basic I/O             
#include <fstream>                                 // C++ file I/O              
#include <iomanip>                                 // C++ I/O formating         
#include <cstring>                                 // C++ strings               

#ifdef falcON_NEMO                                 // compiler option           

#include  <public/nmio.h>                          // nbdy NEMO I/O support     
extern "C" {
# include <stdinc.h>                               // NEMO basics               
}

#endif

////////////////////////////////////////////////////////////////////////////////
#define LoopDims for(register int d=0; d!=falcON_NDIM; ++d)
////////////////////////////////////////////////////////////////////////////////
using namespace nbdy;

#define CAST(TYPE,NAME) static_cast<TYPE>(static_cast<void*>(NAME))
#ifdef falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class BodySrce                                                               
//                                                                              
////////////////////////////////////////////////////////////////////////////////
io BodySrce::read_nemo(                            // R: data read              
		       nemo_in const&input,        // I: nemo input             
		       io      const&want,         //[I: what to read]          
		       bool    const&warn,         //[I: warn upon missing data]
		       size_t  const&begin) const  //[I: first body to get]     
{
  if(want & NEMOBITS() == 0) return io::o;
  if(! input.is_open_set(nemo_io::bodies)) 
    falcON_ErrorF("bodies set not open","BodySrce::read_nemo()");
  if(! is_supported(want&NEMOBITS()) )
    falcON_ErrorF("data types not supported","BodySrce::read_nemo()");
  register int end=input.Number()+begin;
  if(end > NB)
    falcON_ErrorF("too many bodies","BodySrce::read_nemo()");
  register io read = io::o;                        // data read                 
  // read the data                                                              
  // 1 positions and velocities are a bit more involved:                        
  // 1.1 both positions and velocities are required                             
  if(want & io::x && want & io::v) {               // IF  x & v  wanted         
    if(input.is_present(nemo_io::pos) ||           //   IF x or                 
       input.is_present(nemo_io::vel)) {           //      v to read            
      if(input.is_present(nemo_io::pos)) {         //     IF x to read          
	input.read(nemo_io::pos,CAST(real*,POS+begin)); //  read x              
	read |= io::x;                             //       add x to read       
      } else if(warn)                              //     ELSE: warning         
	warning("[BodySrce::read_nemo()]: cannot read: x");
      if(input.is_present(nemo_io::vel)) {         //     IF v to read          
	input.read(nemo_io::vel,CAST(real*,VEL+begin)); //  read v              
	read |= io::v;                             //       add v to read       
      } else if(warn)                              //     ELSE: warning         
	warning("BodySrce::read_nemo()]: cannot read: v");
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(input.is_present(nemo_io::posvel)) { //   ELIF phases present     
      input.read(nemo_io::posvel);                 //     read phases           
      for(register int j=0,i=begin; i!=end; ++i,++j)   // LOOP bodies           
	LoopDims {                                 //       LOOP dims           
	  pos(i)[d] = input.bodies_phs(j)[d];      //         copy x            
	  vel(i)[d] = input.bodies_phs(j)[d+Ndim]; //         copy v            
        }                                          //     END LOOPS             
      read |= io::xv;                              //     add phases to read    
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: (x,v)");
  // 1.2 only positions are required                                            
  } else if(want & io::x) {                        // ELIF  x wanted            
    if(input.is_present(nemo_io::pos)) {           //   IF x present            
      input.read(nemo_io::pos,CAST(real*,POS+begin));  // read x                
      read |= io::x;                               //     add x to read         
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(input.is_present(nemo_io::posvel)) { //   ELIF phases present     
      input.read(nemo_io::posvel);                 //     read phases           
      for(register int j=0,i=begin; i!=end; ++i,++j)   // LOOP bodies           
	LoopDims pos(i)[d]=input.bodies_phs(j)[d]; //       LOOP dims: copy x   
      read |= io::x;                               //     add x to read         
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: x");
  // 1.3 only velocities are required                                           
  } else if(want & io::v) {                        // ELIF  v wanted            
    if(input.is_present(nemo_io::vel)) {           //   IF v present            
      input.read(nemo_io::vel,CAST(real*,VEL+begin));  // read v                
      read |= io::v;                               //     add v to read         
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(input.is_present(nemo_io::posvel)) { //   ELIF phases present     
      input.read(nemo_io::posvel);                 //     read phases           
      for(register int j=0,i=begin; i!=end; ++i,++j)   // LOOP bodies           
	LoopDims                                   //       LOOP dims           
	  vel(i)[d]=input.bodies_phs(j)[d+Ndim];   //         copy v            
      read |= io::v;                               //     add v to read         
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: v");
  }                                                // ENDIF                     
  // 2. read remaining source data                                              
  if(want & io::f) {                               // IF f wanted               
    if(input.is_present(nemo_io::flag)) {          //   IF f present            
      input.read(nemo_io::flag,CAST(int*,FLG+begin));  // read f                
      read |= io::f;                               //     add f to read         
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: f");
  }                                                // ENDIF                     
  if(want & io::m) {                               // IF m wanted               
    if(input.is_present(nemo_io::mass)) {          //   IF m present            
      input.read(nemo_io::mass, MAS+begin);        //     read m                
      read |= io::m;                               //     add m to read         
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: m");
  }                                                // ENDIF                     
  if(want & io::e) {                               // IF e wanted               
    if(input.is_present(nemo_io::eps)) {           //   IF e present            
      input.read(nemo_io::eps, EPS+begin);         //     read e                
      read |= io::e;                               //     add e to read         
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: e");
  }                                                // ENDIF                     
  if(want & io::k) {                               // IF k wanted               
    if(input.is_present(nemo_io::key)) {           //   IF k present            
      input.read(nemo_io::key, KEY+begin);         //     read k                
      read |= io::k;                               //     add k to read         
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: k");
  }                                                // ENDIF                     
  return read;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class BodySink                                                               
//                                                                              
////////////////////////////////////////////////////////////////////////////////
io BodySink::read_nemo(                            // R: data read              
		       nemo_in const&input,        // I: nemo input             
		       io      const&want,         //[I: what to read]          
		       bool    const&warn,         //[I: warn upon missing data]
		       size_t  const&begin) const  //[I: first body to get]     
{
  if(want & NEMOBITS() == 0) return io::o;
  if(! input.is_open_set(nemo_io::bodies)) 
    falcON_ErrorF("bodies set not open","BodySink::read_nemo()");
  if(! is_supported(want&NEMOBITS()))
    falcON_ErrorF("data types not supported","BodySink::read_nemo()");
  register int end=input.Number()+begin;
  if(end > NB)
    falcON_ErrorF("too many bodies","BodySink::read_nemo()");
  register io read = io::o;                        // data read                 
  if(want & io::a) {                               // IF a wanted               
    if(input.is_present(nemo_io::acc)) {           //   IF a present            
      input.read(nemo_io::acc,CAST(real*,ACC+begin)); //  read a                
      read |= io::a;                               //     add a to read         
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySink::read_nemo()]: cannot read: a");
  }                                                // ENDIF                     
  if(want & io::p) {                               // IF p wanted               
    if(input.is_present(nemo_io::pot)) {           //   IF p present            
      input.read(nemo_io::pot, POT+begin);         //     read POT              
      read |= io::p;                               //     add p to read         
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySink::read_nemo()]: cannot read: p");
  }                                                // ENDIF                     
  if(want & io::r) {                               // IF r wanted               
    if(input.is_present(nemo_io::rho)) {           //   IF r present            
      input.read(nemo_io::rho, RHO+begin);         //     read r                
      read |= io::r;                               //     add r to read         
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySink::read_nemo()]: cannot read: r");
  }                                                // ENDIF                     
  if(want & io::y) {                               // IF y wanted               
    if(input.is_present(nemo_io::aux)) {           //   IF y present            
      input.read(nemo_io::aux, AUX+begin);         //     read y                
      read |= io::y;                               //     add y to read         
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySink::read_nemo()]: cannot read: y");
  }                                                // ENDIF                     
  if(want & io::l) {                               // IF l wanted               
    if(input.is_present(nemo_io::level)) {         //   IF l present            
      input.read(nemo_io::level,CAST(short*,LEV+begin)); // read l              
      read |= io::l;                               //     add n to read         
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySink::read_nemo()]: cannot read: l");
  }                                                // ENDIF                     
  if(want & io::n) {                               // IF n wanted               
    if(input.is_present(nemo_io::numb)) {          //   IF n present            
      input.read(nemo_io::numb,CAST(int*,NUM+begin)); // read n                 
      read |= io::n;                               //     add n to read         
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySink::read_nemo()]: cannot read: n");
  }                                                // ENDIF                     
  return read;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class BodyPsph                                                               
//                                                                              
////////////////////////////////////////////////////////////////////////////////
io BodyPsph::read_nemo(                            // R: data read              
		       nemo_in const&input,        // I: nemo input             
		       io      const&want,         //[I: what to read]          
		       bool    const&warn,         //[I: warn upon missing data]
		       size_t  const&begin) const  //[I: first body to get]     
{
  if(want & NEMOBITS() == 0) return io::o;
  if(! input.is_open_set(nemo_io::bodies)) 
    falcON_ErrorF("bodies set not open","BodyPsph::read_nemo()");
  if(! is_supported(want&NEMOBITS()))
    falcON_ErrorF("data types not supported","BodyPsph::read_nemo()");
  register int end=input.NumberSPH()+begin;
  if(end > NB)
    falcON_ErrorF("too many bodies","BodyPsph::read_nemo()");
  register io read = io::o;                        // data read                 
  // read SPH data                                                              
  if(want & io::H) {                               // IF H wanted               
    if(input.is_present(nemo_io::h)) {             //   IF H present            
      input.read(nemo_io::h, SIZ+begin);           //     read H                
      read |= io::H;                               //     add H to read         
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: H");
  }                                                // ENDIF                     
#ifdef falcON_SPH
  if(want & io::N) {                               // IF N wanted               
    if(input.is_present(nemo_io::numbSPH)) {       //   IF N present            
      input.read(nemo_io::numbSPH, CAST(int*,NSP+begin)); //     read N         
      read |= io::N;                               //     add N to read         
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: N");
  }                                                // ENDIF                     
  if(want & io::U) {                               // IF U wanted               
    if(input.is_present(nemo_io::uin)) {           //   IF U present            
      input.read(nemo_io::uin, UIN+begin);         //     read U                
      read |= io::U;                               //     add U to read         
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: U");
  }                                                // ENDIF                     
  if(want & io::I) {                               // IF I wanted               
    if(input.is_present(nemo_io::udin)) {          //   IF I present            
      input.read(nemo_io::udin, UDI+begin);        //     read I                
      read |= io::I;                               //     add I to read         
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: I");
  }                                                // ENDIF                     
  if(want & io::E) {                               // IF E wanted               
    if(input.is_present(nemo_io::udex)) {          //   IF E present            
      input.read(nemo_io::udex, UDE+begin);        //     read E                
      read |= io::E;                               //     add E to read         
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: E");
  }                                                // ENDIF                     
  if(want & io::S) {                               // IF S wanted               
    if(input.is_present(nemo_io::entr)) {          //   IF S present            
      input.read(nemo_io::entr, ENT+begin);        //     read S                
      read |= io::S;                               //     add S to read         
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: S");
  }                                                // ENDIF                     
  if(want & io::R) {                               // IF R wanted               
    if(input.is_present(nemo_io::srho)) {          //   IF R present            
      input.read(nemo_io::srho, SRH+begin);        //     read R                
      read |= io::R;                               //     add R to read         
      DATA_CHANGED = 1;                            //     mark data change      
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: R");
  }                                                // ENDIF                     
#endif
  return read;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class sbodies:                                                               
//                                                                              
////////////////////////////////////////////////////////////////////////////////
bool
sbodies::read_nemo_particles(                      // R: was time in range?     
			     nemo_in const&input,  // I: nemo input             
			     io           &read,   // O: what has been read     
			     real         *time,   //[O: time]                  
			     const io      want,   //[I: what to read]          
			     char         *times,  //[I: time range]            
			     const bool    warn)   //[I: warn: missing data]    
{
  register real t;                                 // time                      
  // 1. open parameter set; read N & time; close parameter set;                 
  input.open_set(nemo_io::param);                  // open nemo parameter set   
  input.read_N();                                  //   read N & N_sph          
  if(input.is_present(nemo_io::time))              //   IF time available       
    input.read(nemo_io::time,&t);                  //     read time             
  else {                                           //   ELSE                    
    t = zero;                                      //     default t=0           
    if(time && warn)                               //     IF time wanted: warn  
      warning("[%s]: no time found in snapshot; defaulting to 0",
	      "sbodies::read_nemo_particles()");
  };                                               //   ENDIF                   
  if(time) *time = t;                              //   IF wanted, set time     
  input.close_set(nemo_io::param);                 // close nemo parameter set  
  // 2. check for time range and data support                                   
  if(! time_in_range(t,times)) return false;       // IF(time not in range) DONE
  register io get = want & NEMOBITS();             // bits act'ly NEMO readable 
  if(want != get && warn)
    warning("[%s]: not all data desired are NEMO readable",
	    "sbodies::read_nemo_particles()");
  if(want == io::o) return true;                   // IF(nothing to read) DOME  
  // 3 read data                                                                
  input.open_set(nemo_io::bodies);                 // open nemo particle set    
  reset(input.Number(), get|all_bits(),
	input.NumberSPH());                        // reset data                
  read = BodySrce::read_nemo(input,get,warn);      // read source data          
  read|= BodySink::read_nemo(input,get,warn);      // read sink data            
  read|= BodyPsph::read_nemo(input,get,warn);      // read SPH data             
  input.close_set(nemo_io::bodies);                // close nemo particle set   
  input.reset();                                   // reset nemo input          
  return true;                                     // return                    
}
//------------------------------------------------------------------------------
bool
sbodies::read_nemo_snapshot (                      // R: was time in range?     
			     nemo_in const&input,  // I: nemo input             
			     io           &read,   // O: what has been read     
			     real         *time,   //[O: time]                  
			     const io      want,   //[I: what to read]          
			     char         *times,  //[I: time range]            
			     const bool    warn)   //[I: warn: missing data]    
{
  if(!input.is_present(nemo_io::snap))             // IF no snapshot present    
    falcON_ErrorF("no snapshot found",
		  "sbodies::read_nemo_snapshot()");// ERROR   
  input.open_set(nemo_io::snap);                   // open nemo snapshot        
  bool tt = read_nemo_particles(input,read,time,   // read particles if         
				want,times,warn);  //   time in range           
  input.close_set(nemo_io::snap);                  // close nemo snapshot       
  return tt;                                       // return (time in range)?   
}
//------------------------------------------------------------------------------
bool
sbodies::read_nemo_particles(                      // R: was time in range?     
			     const nemo_in*in,     // I: nemo inputs            
			     int     const&Nin,    // I: # nemo inputs          
			     io           *read,   // O: what has been read     
			     uint         *Nnum,   // O: how many have been read
			     real         *time,   //[O: time]                  
			     const io     *want,   //[I: what to read]          
			     char         *times,  //[I: time range]            
			     const bool    warn)   //[I: warn: missing data]    
{
  register bool time_read=false;                   // time read already?        
  register real t=0,ti;                            // time                      
  register uint NN=0u,NNS=0u;                      // total number to read      
  register io   get=io::mxv;                       // bits act'ly NEMO readable 
  // 1. open parameter sets; read N & time; close parameter sets;               
  for(register int i=0; i!=Nin; ++i) {             // LOOP nemo inputs          
    in[i].open_set(nemo_io::param);                //   open nemo parameter set 
    in[i].read_N();                                //   read N & N_sph          
    NN += Nnum[i] = in[i].Number();                //   accumulate N            
    NNS+= in[i].NumberSPH();                       //   accumulate N_sph        
    if(in[i].is_present(nemo_io::time)) {          //   IF time available       
      in[i].read(nemo_io::time,&ti);               //     read time             
      if(time_read && ti != t)                     //     IF time differs       
	falcON_ErrorF("different times found in snapshots",//  issue ERROR      
		      "sbodies::read_nemo_particles()");
      t = ti;                                      //     set time              
      time_read = true;                            //     set flag: time read   
    }                                              //   ENDIF                   
    in[i].close_set(nemo_io::param);               //   close nemo parameter set
    if(want) get |= want[i];
  }                                                // END LOOP                  
  if(time) {                                       // IF time wanted            
    *time = t;                                     //   assign to output value  
    if(!time_read)                                 //   IF(time not read)       
      warning("[%s]: no time found in snapshots; defaulting to 0",
	      "sbodies::read_nemo_particles()");
  }                                                // ENDIF                     
  // 2. check for time range and data support                                   
  if(! time_in_range(t,times)) return false;       // IF(time not in range) DONE
  if(warn && (get & NEMOBITS()) != get)
    warning("[%s]: not all data desired are NEMO readable",
	    "sbodies::read_nemo_particles()");
  get &= NEMOBITS();                               // bits to allocate          
  if(get == io::o) return true;                    // IF(nothing to read) DOME  
  // 3 read data                                                                
  reset(NN, get|all_bits(),NNS);                   // reset data                
  register size_t begin=0;                         // first body to read in     
  register io wanted;
  for(register int i=0; i!=Nin; ++i) {             // LOOP nemo inputs          
    in[i].open_set(nemo_io::bodies);               //   open nemo particle set  
    if(want) wanted=want[i]&NEMOBITS();            //   data to be read         
    else     wanted=io::mxv;
    read[i] = BodySrce::read_nemo(in[i],wanted,warn,begin); // read source data 
    read[i]|= BodySink::read_nemo(in[i],wanted,warn,begin); // read sink data   
    read[i]|= BodyPsph::read_nemo(in[i],wanted,warn,begin); // read sph data    
    begin  += in[i].Number();                      //   increment begin         
    in[i].close_set(nemo_io::bodies);              //   close nemo particle set 
    in[i].reset();                                 //   reset nemo input        
  }                                                // END LOOP                  
  return true;                                     // return                    
}
//------------------------------------------------------------------------------
bool
sbodies::read_nemo_snapshots(                      // R: was time in range?     
			     const nemo_in*in,     // I: nemo inputs            
			     int     const&Nin,    // I: # nemo inputs          
			     io           *read,   // O: what has been read     
			     uint         *Nnum,   // O: how many have been read
			     real         *time,   //[O: time]                  
			     const io     *want,   //[I: what to read]          
			     char         *times,  //[I: time range]            
			     const bool    warn)   //[I: warn: missing data]    
{
  for(register int i=0; i!=Nin; ++i) {             // LOOP nemo inputs          
    if(!in[i].is_present(nemo_io::snap))           //   IF no snapshot present  
      falcON_ErrorF("no snapshot found",
		    "sbodies::read_nemo_snapshot()"); // ERROR 
    in[i].open_set(nemo_io::snap);                 // open nemo snapshot        
  }                                                // END LOOP                  
  bool tt=read_nemo_particles(in,Nin,read,Nnum,    // read particles IF         
			      time,want,times,warn); // time in range           
  for(register int i=0; i!=Nin; ++i)               // LOOP nemo inputs          
    in[i].close_set(nemo_io::snap);                //   open nemo snapshot      
  return tt;                                       // return (time in range)?   
}
//------------------------------------------------------------------------------
void                                               // write bodies to output    
sbodies::write_nemo_particles(nemo_out const&out,  // I: nemo output            
			      const real    *time, //[I: write time]            
			      io   const    &write,//[I: what to write]         
			      uint const    &K,    //[I: only write first K]    
			      uint const    &begin)//[I: begin with this]       
  const
{
  register uint i,j;                               // index                     
  register uint Nout = N_bodies() > begin?         // # bodies to write out     
    (K? min(begin+K,N_bodies()) - begin : N_bodies()-begin) : 0;
  register uint Sout = N_sph   () > begin?         // # SPH bodies to write out 
    (K? min(begin+K,N_sph   ()) - begin : N_sph   ()-begin) : 0;
  // 1. open parameter set; write out N, NS & time; close parameter set;        
  out.open_set(nemo_io::param);                    // open a parameter set      
  out.write_N(Nout,Sout);                          //   write N, NS             
  if(time)                                         //   IF(time to be written)  
    out.write(nemo_io::time,*time);                //     write out time        
  out.close_set(nemo_io::param);                   // close parameter set       
  // 2. allocate memory for output; open particle set; output data;             
  if(write == io::o) return;                       // nothing to be done        
  if(Nout==0u && Sout==0) return;                  // nothing to be done        
  out.open_set(nemo_io::bodies);                   // open particles set        
  out.write(nemo_io::cart);                        // put out coord system      
  // 2.1 output masses                                                          
  if(write & io::m && has(io::m))                  // IF(masses to be written)  
    out.write(nemo_io::mass,MAS+begin);            //   write out masses        
  // 2.2 output positions                                                       
  if(write & io::x && has(io::x))                  // IF poss for output        
    out.write(nemo_io::pos,CAST(real*,POS+begin)); //   write out poss          
  // 2.3 output velocities                                                      
  if(write & io::v && has(io::v))                  // IF vels for output        
    out.write(nemo_io::vel,CAST(real*,VEL+begin)); //   write out vels          
  // 2.4 output accelerations                                                   
  if(write & io::a && has(io::a))                  // IF accs for output        
    out.write(nemo_io::acc,CAST(real*,ACC+begin)); //   write out accs          
  // 2.5 output potential(s)                                                    
  // 2.5.1 N-body potential + external potential                                
  if(write&io::p && has(io::p) &&
     write&io::q && has(io::q)) {                  // IF N-body + external pot  
    out.allocscalar();                             //   allocate memory for P+p 
    for(j=0,i=begin; j!=Nout; ++i,++j)             //   LOOP bodies             
      out.bodies_scl(j) = pot(i) + pex(i);         //     add up total potential
    out.write(nemo_io::pot);                       //   write out potential     
  // 2.5.2 just N-body potential                                                
  } else if(write & io::p && has(io::p))           // ELIF only N-body pot      
    out.write(nemo_io::pot,POT+begin);             //   write out potential     
  // 2.5.3 just external potential                                              
  else if(write & io::q && has(io::q))             // ELIF only external pot    
    out.write(nemo_io::pot,PEX+begin);             //   write external pot      
  // 2.6 output density                                                         
  if(write & io::r && has(io::r))                  // IF(rhos for output)       
    out.write(nemo_io::rho,RHO+begin);             //     write out scalars     
  // 2.7 output aux data                                                        
  if(write & io::y && has(io::y))                  // IF(auxx for output)       
    out.write(nemo_io::aux,AUX+begin);             //   write out scalars       
  // 2.8 output eps_i                                                           
  if(write & io::e && has(io::e))                  // IF(eps_ifor output)       
    out.write(nemo_io::eps,EPS+begin);             //   write out potential     
  // 2.9 output level data                                                      
  if(write & io::l && has(io::l))                  // IF(level for output)      
    out.write(nemo_io::level,CAST(short*,LEV+begin)); // write out shorts       
  // 2.10 output num data                                                       
  if(write & io::n && has(io::n))                  // IF(num for output)        
    out.write(nemo_io::numb,CAST(int*,NUM+begin)); // write out ints            
  // 2.11 output flag data                                                      
  if(write & io::f && has(io::f))                  // IF(flag for output)       
    out.write(nemo_io::flag,CAST(int*,FLG+begin)); // write out flags           
  // 2.12 output key data                                                       
  if(write & io::k && has(io::k))                  // IF(key for output)        
    out.write(nemo_io::key,KEY+begin);             //   write out keys          
  // 2.13 output h_i data                                                       
  if(write & io::H && has(io::H))                  // IF(key for output)        
    out.write(nemo_io::h,SIZ+begin);               //   write out keys          
#ifdef falcON_SPH
  // 2.14 output N data                                                         
  if(write & io::N && has(io::N))                  // IF(N for output)          
    out.write(nemo_io::numbSPH,CAST(int*,NSP+begin)); //  write out N           
  // 2.15 output U_intern                                                       
  if(write & io::U && has(io::U))                  // IF(internal U for output) 
    out.write(nemo_io::uin,UIN+begin);             //   write out internal U    
  // 2.16 output (dU_intern/dt)_intern                                          
  if(write & io::I && has(io::I))                  // IF(internal I for output) 
    out.write(nemo_io::udin,UDI+begin);            //   write out dU/dt_intern  
  // 2.17 output (dU_intern/dt)_extern                                          
  if(write & io::E && has(io::E))                  // IF(internal I for output) 
    out.write(nemo_io::udex,UDE+begin);            //   write out dU/dt_extern  
  // 2.18 output entropy                                                        
  if(write & io::S && has(io::S))                  // IF(entrpy for output)     
    out.write(nemo_io::entr,ENT+begin);            //   write out entropies     
  // 2.19 output gas-density                                                    
  if(write & io::R && has(io::R))                  // IF(entrpy for output)     
    out.write(nemo_io::srho,SRH+begin);            //   write out entropies     
#endif
  // 3. close particle set & reset nemo output                                  
  out.close_set(nemo_io::bodies);                  // close particles set       
  out.reset();                                     // de-allocate arrays        
}
//------------------------------------------------------------------------------
void                                               // write snapshot            
sbodies::write_nemo_snapshot(nemo_out const&out,   // I: nemo output            
			     const real    *time,  //[I: write time]            
			     io       const&write, //[I: write pot/acc etc?]    
			     uint     const&K,     //[I: only write first K]    
			     uint     const&begin) //[I: begin with this]       
  const
{
  out.open_set(nemo_io::snap);                     // open a new nemo snapshot  
  write_nemo_particles(out,time,write,K,begin);    //   write particle set      
  out.close_set(nemo_io::snap);                    // close nemo snapshot       
}
////////////////////////////////////////////////////////////////////////////////
#endif // ALLLOW_NEMO
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class sbodies: YANC output                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  struct BodyWrite {
    static io            WriteThis;
    static std::ostream *OStream;
  };
  io            BodyWrite::WriteThis;
  std::ostream* BodyWrite::OStream;
  //----------------------------------------------------------------------------
  template<int IO> struct BodyWriteAscii :
    protected BodyWrite, protected nbody_io<IO,sbodies> {
    static void act_on_body(sbodies::body const&B) {
      if(WriteThis & IO)
	(*OStream) << const_access(B) <<' ';
    } };
  //----------------------------------------------------------------------------
  template<int IO> struct BodyWriteBinary :
    protected BodyWrite, protected nbody_io<IO,sbodies> {
    static void act_on_bodies(const sbodies*const&B) {
      if(WriteThis & IO)
	OStream->write(CAST(char*,array_access(B)),number(B)*sizeof(io_type));
    } };
}
//------------------------------------------------------------------------------
void sbodies::write_yanc_ascii(std::ostream &out, const io write) const
{
  if(write == io::o) return;
  BodyWrite::WriteThis = write & all_bits();
  BodyWrite::OStream   = &out;
  out<<BodyWrite::WriteThis<<std::endl;
  LoopBodies(sbodies,this,Bi) {
    LoopIO<BodyWriteAscii>::loop(Bi);
    out<<'\n';
  }
  out.flush();
}
//------------------------------------------------------------------------------
void sbodies::write_yanc_binary(std::ostream&out, const io write) const
{
  if(write == io::o) return;
  BodyWrite::WriteThis = write & all_bits();
  BodyWrite::OStream   = &out;
  out.write(CAST(char*,&(BodyWrite::WriteThis)),sizeof(io));
  LoopIO<BodyWriteBinary>::loop(this);
  out.flush();
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class sbodies: YANC input                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  struct BodyRead {
    static io            ReadThis, Present;
    static std::istream *IStream;
  };
  io            BodyRead::ReadThis, BodyRead::Present;
  std::istream* BodyRead::IStream;
  //----------------------------------------------------------------------------
  template<int IO> struct BodyReadAscii :
    protected BodyRead, protected nbody_io<IO,sbodies> {
    static void act_on_body(sbodies::body &B) {
      if(Present & IO)
	if(ReadThis & IO)
	  (*IStream) >> access(B);
	else {
	  typename nbody_io<IO,sbodies>::io_type dummy;
	  (*IStream) >> dummy;
	}
    } };
  //----------------------------------------------------------------------------
  template<int IO> struct BodyReadBinary :
    protected BodyRead, protected nbody_io<IO,sbodies> {
    static void act_on_bodies(const sbodies*const&B) {
      if(Present&IO) {
	if(ReadThis&IO)
	  IStream->read(CAST(char*,array_access(B)),number(B)*sizeof(io_type));
	else
	  IStream->seekg(number(B)*sizeof(io_type));
      }
    } };
}
//------------------------------------------------------------------------------
io sbodies::read_yanc_ascii(std::istream& in, const io want)
  // NOTE: reset(N,data) has to be called BEFORE.                               
  //       We assume the correct N and try to read all data wanted              
  //       Unlike previous versions, we cannot "enable" holding aux data        
{
  if(!all_bits().contains(want)) {
    char miss[30];
    all_bits().missing(want).make_word(miss);
    warning("[sbodies::read_yanc_ascii()]: "
	    "cannot read %s: no memory allocated",miss);
  }
  in >> BodyRead::Present;
  if(!BodyRead::Present.contains(want)) {
    char miss[30];
    BodyRead::Present.missing(want).make_word(miss);
    warning("[sbodies::read_yanc_ascii()]: "
	    "cannot read %s: not present in file",miss);
  }
  BodyRead::ReadThis = want & BodyRead::Present & all_bits();
  BodyRead::IStream  = &in;
  reset(N_bodies(), BodyRead::ReadThis | all_bits(), N_sph());
  if(BodyRead::ReadThis & io::sphmax) {
    if(N_sph())
      LoopSPHBodies(sbodies,this,Bi) {
        LoopIO<BodyReadAscii>::loop(Bi);
	SwallowRestofLine(in);
      }  
    else
      warning("[sbodies::read_yanc_ascii()]: cannot read SPH data");
    LoopNonSPHBodies(sbodies,this,Bi) {
      LoopIO<BodyReadAscii,0,IO_NOSPH>::loop(Bi);
      SwallowRestofLine(in);
    }
  } else
    LoopBodies(sbodies,this,Bi) {
      LoopIO<BodyReadAscii,0,IO_NOSPH>::loop(Bi);
      SwallowRestofLine(in);
    }
  return BodyRead::ReadThis;
}
//------------------------------------------------------------------------------
void sbodies::read_yanc_binary(std::istream&in, const io want)
{
  if(!all_bits().contains(want)) {
    char miss[30];
    all_bits().missing(want).make_word(miss);
    warning("[sbodies::read_yanc_ascii()]: "
	    "cannot read %s: no memory allocated",miss);
  }
  in >> BodyRead::Present;
  if(!BodyRead::Present.contains(want)) {
    char miss[30];
    BodyRead::Present.missing(want).make_word(miss);
    warning("[sbodies::read_yanc_ascii()]: "
	    "cannot read %s: not present in file",miss);
  }
  BodyRead::ReadThis = want & BodyRead::Present & all_bits();
  BodyRead::IStream  = &in;
  reset(N_bodies(), BodyRead::ReadThis | all_bits(), N_sph());
  if((BodyRead::ReadThis & io::sphmax) && N_sph()==0)
    warning("[sbodies::read_yanc_binary()]: cannot read SPH data");
  LoopIO<BodyReadBinary>::loop(this);
}
//------------------------------------------------------------------------------
namespace {
  template<int IO> struct BodyReadSimple :
    protected BodyRead, protected nbody_io<IO,sbodies> {
    static void read(std::istream&in, sbodies::body &B) {
      in >> access(B); }
    };
  //----------------------------------------------------------------------------
  template<int IO>
  void Read(std::istream&in, sbodies::body &B) {
    BodyReadSimple<IO>::read(in,B);
  }
  //----------------------------------------------------------------------------
  void Ignore(std::istream&in, sbodies::body &B) { }
  //----------------------------------------------------------------------------
  typedef void (*p_reader)(std::istream&, sbodies::body &);
}
//------------------------------------------------------------------------------
void sbodies::read_simple_ascii(std::istream   &in,
				const io* const&item,
				uint      const&Nb,
				uint      const&Ns)
{
  int K;
  io  get=io::o;
  p_reader readall[32]={0}, readnon[32]={0};
  for(K=0; K!=32 && item[K]!=io::o; ++K) {
    get |= item[K];
    switch(item[K]) {
    case io::m: readall[K]=&Read<io::m>; readnon[K]=&Read<io::m>; break;
    case io::x: readall[K]=&Read<io::x>; readnon[K]=&Read<io::x>; break;
    case io::v: readall[K]=&Read<io::v>; readnon[K]=&Read<io::v>; break;
    case io::e: readall[K]=&Read<io::e>; readnon[K]=&Read<io::e>; break;
    case io::f: readall[K]=&Read<io::f>; readnon[K]=&Read<io::f>; break;
    case io::k: readall[K]=&Read<io::k>; readnon[K]=&Read<io::k>; break;
    case io::p: readall[K]=&Read<io::p>; readnon[K]=&Read<io::p>; break;
    case io::q: readall[K]=&Read<io::q>; readnon[K]=&Read<io::q>; break;
    case io::a: readall[K]=&Read<io::a>; readnon[K]=&Read<io::a>; break;
    case io::r: readall[K]=&Read<io::r>; readnon[K]=&Read<io::r>; break;
    case io::y: readall[K]=&Read<io::y>; readnon[K]=&Read<io::y>; break;
    case io::l: readall[K]=&Read<io::l>; readnon[K]=&Read<io::l>; break;
    case io::n: readall[K]=&Read<io::n>; readnon[K]=&Read<io::n>; break;
    case io::H: readall[K]=&Read<io::H>; readnon[K]=&Ignore;      break;
#ifdef falcON_SPH
    case io::N: readall[K]=&Read<io::N>; readnon[K]=&Ignore;      break;
    case io::U: readall[K]=&Read<io::U>; readnon[K]=&Ignore;      break;
    case io::Y: readall[K]=&Read<io::Y>; readnon[K]=&Ignore;      break;
    case io::I: readall[K]=&Read<io::I>; readnon[K]=&Ignore;      break;
    case io::E: readall[K]=&Read<io::E>; readnon[K]=&Ignore;      break;
    case io::S: readall[K]=&Read<io::S>; readnon[K]=&Ignore;      break;
    case io::R: readall[K]=&Read<io::R>; readnon[K]=&Ignore;      break;
    case io::D: readall[K]=&Read<io::D>; readnon[K]=&Ignore;      break;
    case io::V: readall[K]=&Read<io::V>; readnon[K]=&Ignore;      break;
    case io::Q: readall[K]=&Read<io::Q>; readnon[K]=&Ignore;      break;
    case io::T: readall[K]=&Read<io::T>; readnon[K]=&Ignore;      break;
    case io::J: readall[K]=&Read<io::J>; readnon[K]=&Ignore;      break;
#endif
    }
  }
  reset(Nb,get|all_bits(),Ns);                     // reset data                
  LoopSPHBodies(sbodies,this,Bi) {
    if(!in)
      falcON_ErrorF("end of input before data have been read",
		    "sbodies::read_simple_ascii()");
    if(! eat_line(in,'#') ) {
      for(register p_reader* R=readall; *R; ++R) (*R)(in,Bi);
      SwallowRestofLine(in);
    }
  }
  LoopNonSPHBodies(sbodies,this,Bi) {
    if(!in)
      falcON_ErrorF("end of input before data have been read",
		    "sbodies::read_simple_ascii()");
    if(! eat_line(in,'#') ) {
      for(register p_reader* R=readnon; *R; ++R) (*R)(in,Bi);
      SwallowRestofLine(in);
    }
  }
}
