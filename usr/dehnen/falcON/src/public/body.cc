// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// body.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::io                                                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
const char* io::word() const {
  static char w[20];
  make_word(w);
  return w;
}
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
    } else if(input.is_present(nemo_io::posvel)) { //   ELIF phases present     
      input.read(nemo_io::posvel);                 //     read phases           
      for(register int j=0,i=begin; i!=end; ++i,++j)   // LOOP bodies           
	LoopDims {                                 //       LOOP dims           
	  pos(i)[d] = input.bodies_phs(j)[d];      //         copy x            
	  vel(i)[d] = input.bodies_phs(j)[d+Ndim]; //         copy v            
        }                                          //     END LOOPS             
      read |= io::xv;                              //     add phases to read    
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: (x,v)");
  // 1.2 only positions are required                                            
  } else if(want & io::x) {                        // ELIF  x wanted            
    if(input.is_present(nemo_io::pos)) {           //   IF x present            
      input.read(nemo_io::pos,CAST(real*,POS+begin));  // read x                
      read |= io::x;                               //     add x to read         
    } else if(input.is_present(nemo_io::posvel)) { //   ELIF phases present     
      input.read(nemo_io::posvel);                 //     read phases           
      for(register int j=0,i=begin; i!=end; ++i,++j)   // LOOP bodies           
	LoopDims pos(i)[d]=input.bodies_phs(j)[d]; //       LOOP dims: copy x   
      read |= io::x;                               //     add x to read         
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: x");
  // 1.3 only velocities are required                                           
  } else if(want & io::v) {                        // ELIF  v wanted            
    if(input.is_present(nemo_io::vel)) {           //   IF v present            
      input.read(nemo_io::vel,CAST(real*,VEL+begin));  // read v                
      read |= io::v;                               //     add v to read         
    } else if(input.is_present(nemo_io::posvel)) { //   ELIF phases present     
      input.read(nemo_io::posvel);                 //     read phases           
      for(register int j=0,i=begin; i!=end; ++i,++j)   // LOOP bodies           
	LoopDims                                   //       LOOP dims           
	  vel(i)[d]=input.bodies_phs(j)[d+Ndim];   //         copy v            
      read |= io::v;                               //     add v to read         
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: v");
  }                                                // ENDIF                     
  // 2. read remaining source data                                              
  if(want & io::f) {                               // IF f wanted               
    if(input.is_present(nemo_io::flag)) {          //   IF f present            
      input.read(nemo_io::flag,CAST(int*,FLG+begin));  // read f                
      read |= io::f;                               //     add f to read         
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: f");
  }                                                // ENDIF                     
  if(want & io::m) {                               // IF m wanted               
    if(input.is_present(nemo_io::mass)) {          //   IF m present            
      input.read(nemo_io::mass, MAS+begin);        //     read m                
      read |= io::m;                               //     add m to read         
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: m");
  }                                                // ENDIF                     
  if(want & io::e) {                               // IF e wanted               
    if(input.is_present(nemo_io::eps)) {           //   IF e present            
      input.read(nemo_io::eps, EPS+begin);         //     read e                
      read |= io::e;                               //     add e to read         
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySrce::read_nemo()]: cannot read: e");
  }                                                // ENDIF                     
  if(want & io::k) {                               // IF k wanted               
    if(input.is_present(nemo_io::key)) {           //   IF k present            
      input.read(nemo_io::key, KEY+begin);         //     read k                
      read |= io::k;                               //     add k to read         
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
      read |= io::l;                               //     add l to read         
    } else if(warn)                                //   ELSE: warning           
      warning("[BodySink::read_nemo()]: cannot read: l");
  }                                                // ENDIF                     
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
  register uint NN = input.read_N();               //   read N                  
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
  reset(NN,get | all_bits());                      // reset data                
  read = BodySrce::read_nemo(input,get,warn);      // read source data          
  read|= BodySink::read_nemo(input,get,warn);      // read sink data            
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
  register uint NN=0u;                             // total number to read      
  register io   get=io::mxv;                       // bits act'ly NEMO readable 
  // 1. open parameter sets; read N & time; close parameter sets;               
  for(register int i=0; i!=Nin; ++i) {             // LOOP nemo inputs          
    in[i].open_set(nemo_io::param);                //   open nemo parameter set 
    NN += Nnum[i] = in[i].read_N();                //   read & accumulate N     
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
  reset(NN,get | all_bits());                      // reset data                
  register size_t begin=0;                         // first body to read in     
  register io wanted;
  for(register int i=0; i!=Nin; ++i) {             // LOOP nemo inputs          
    in[i].open_set(nemo_io::bodies);               //   open nemo particle set  
    if(want) wanted=want[i]&NEMOBITS();            //   data to be read         
    else     wanted=io::mxv;
    read[i] = BodySrce::read_nemo(in[i],wanted,warn,begin); // read source data 
    read[i]|= BodySink::read_nemo(in[i],wanted,warn,begin); // read sink data   
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
sbodies::write_nemo_particles(nemo_out  &output,   // I: nemo output            
			      const real*time,     //[I: write time]            
			      const io   write,    //[I: what to write]         
			      const uint K,        //[I: only write first K]    
			      const uint begin)    //[I: begin with this]       
  const
{
  register uint i,j;                               // index                     
  register uint Nout = K? K : N_bodies();          // # bodies to write out     
  // 1. open parameter set; write out N & time; close parameter set;            
  output.open_set(nemo_io::param);                 // open a parameter set      
  output.write_N(Nout);                            //   write N                 
  if(time)                                         //   IF(time to be written)  
    output.write(nemo_io::time,*time);             //     write out time        
  output.close_set(nemo_io::param);                // close parameter set       
  // 2. allocate memory for output; open particle set; output data;             
  if(write == io::o) return;                       // nothing to be done        
  if(N_bodies() ==0) return;                       // nothing to be done        
  output.open_set(nemo_io::bodies);                // open particles set        
  output.write(nemo_io::cart);                     // put out coord system      
  // 2.1 output masses                                                          
  if(write & io::m && has(io::m))                  // IF(masses to be written)  
    output.write(nemo_io::mass,MAS+begin);         //   write out masses        
  // 2.2 - 2.3 output of positions & velocities                                 
#ifdef  falcON_OUTPUT_PHASES  // as phase-space coordinates [x,v]_i             
  // 2.2.1 both positions and velocities to be written out                      
  if(write & io::x && has(io::x) &&
     write & io::v && has(io::v)) {                // IF [x,v] to be written    
    output.allocarrays();                          //   allocate memory for x,v 
    for(j=0,i=begin; j!=Nout; ++i,++j) LoopDims {  //   LOOP bodies & dims      
      output.bodies_phs(j)[d]      = pos(i)[d];    //     copy x                
      output.bodies_phs(j)[d+Ndim] = vel(i)[d];    //     copy v                
    }                                              //   END LOOP                
    output.write(nemo_io::posvel);                 //   write out phases        
  // 2.2.2 only positions to be written out                                     
  } else if(write & io::x && has(io::x)) {         // ELIF x to be written      
    output.allocarrays();                          //   allocate memory for x,v 
    for(j=0,i=begin; j!=Nout; ++i,++j) LoopDims    //   LOOP bodies             
      output.bodies_phs(j)[d]      = pos(i)[d];    //     copy x                
    output.write(nemo_io::posvel);                 //   write out phases        
  // 2.2.3 only velocities to be written out                                    
  } else if(write & io::v && has(io::v)) {         // ELIF v to be written      
    output.allocarrays();                          //   allocate memory for x,v 
    for(j=0,i=begin; j!=Nout; ++i,++j) LoopDims    //   LOOP bodies             
      output.bodies_phs(j)[d+Ndim] = vel(i)[d];    //     copy v                
    output.write(nemo_io::posvel);                 //   write out phases        
  }                                                // ENDIF                     
#else                         // separately as positions x_i and velocities v_i 
  // 2.2 output positions                                                       
  if(write & io::x && has(io::x))                  // IF poss for output        
    output.write(nemo_io::pos,CAST(real*,POS+begin));  // write out poss        
  // 2.3 output velocities                                                      
  if(write & io::v && has(io::v))                  // IF vels for output        
    output.write(nemo_io::vel,CAST(real*,VEL+begin));  // write out vels        
#endif
  // 2.4 output accelerations                                                   
  if(write & io::a && has(io::a))                  // IF accs for output        
    output.write(nemo_io::acc,CAST(real*,ACC+begin));  // write out accs        
  // 2.5 output potential(s)                                                    
  // 2.5.1 N-body potential + external potential                                
  if(write&io::p && has(io::p) &&
     write&io::P && has(io::P)) {                  // IF N-body + external pot  
    output.allocscalar();                          //   allocate memory for P+p 
    for(j=0,i=begin; j!=Nout; ++i,++j)             //   LOOP bodies             
      output.bodies_scl(j) = pot(i) + pex(i);      //     add up total potential
    output.write(nemo_io::pot);                    //   write out potential     
  // 2.5.2 just N-body potential                                                
  } else if(write & io::p && has(io::p))           // ELIF only N-body pot      
    output.write(nemo_io::pot,POT+begin);          //   write out potential     
  // 2.5.3 just external potential                                              
  else if(write & io::P && has(io::P))             // ELIF only external pot    
    output.write(nemo_io::pot,PEX+begin);          //   write external pot      
  // 2.6 output density                                                         
  if(write & io::r && has(io::r))                  // IF(rhos for output)       
    output.write(nemo_io::rho,RHO+begin);          //     write out scalars     
  // 2.7 output aux data                                                        
  if(write & io::y && has(io::y))                  // IF(auxx for output)       
    output.write(nemo_io::aux,AUX+begin);          //   write out scalars       
  // 2.8 output eps_i                                                           
  if(write & io::e && has(io::e))                  // IF(eps_ifor output)       
    output.write(nemo_io::eps,EPS+begin);          //   write out potential     
  // 2.9 output level data                                                      
  if(write & io::l && has(io::l))                  // IF(level for output)      
    output.write(nemo_io::level,CAST(short*,LEV+begin)); // write out shorts    
  // 2.10 output flag data                                                      
  if(write & io::f && has(io::f))                  // IF(flag for output)       
    output.write(nemo_io::flag,CAST(int*,FLG+begin)); // write out flags        
  // 2.11 output key data                                                       
  if(write & io::k && has(io::k))                  // IF(key for output)        
    output.write(nemo_io::key,KEY+begin);          //   write out keys          
  // 3. close particle set & reset nemo output                                  
  output.close_set(nemo_io::bodies);               // close particles set       
  output.reset();                                  // de-allocate arrays        
}
//------------------------------------------------------------------------------
void                                               // write snapshot            
sbodies::write_nemo_snapshot(nemo_out  &output,    // I: nemo output            
			     const real*time,      //[I: write time]            
			     const io   write,     //[I: write pot/acc etc?]    
			     const uint K,         //[I: only write first K]    
			     const uint begin)     //[I: begin with this]       
  const
{
  output.open_set(nemo_io::snap);                  // open a new nemo snapshot  
  write_nemo_particles(output,time,write,K,begin); //   write particle set      
  output.close_set(nemo_io::snap);                 // close nemo snapshot       
}
////////////////////////////////////////////////////////////////////////////////
#endif // ALLLOW_NEMO
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class sbodies: YANC I/O                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  inline io writer(io const&write, const sbodies*const&B)
  {
    register io written = io::o;
    for(register int i=0,bit=1; i!=io::N_quant(); ++i,bit<<=1)
      if(write & bit && B->all_bits_contain(bit)) written |= bit;
    return written;
  }
}
//------------------------------------------------------------------------------
void sbodies::write_yanc_ascii(std::ostream &out, const io write) const
{
  if(write == io::o) return;
  register io   written = writer(write,this);
  out<<written<<std::endl;
  for(register uint i=0; i!=N_bodies(); ++i) {
    if(written & io::m) out << mass (i) <<" ";
    if(written & io::x) out << pos  (i) <<" ";
    if(written & io::v) out << vel  (i) <<" ";
    if(written & io::e) out << eps  (i) <<" ";
    if(written & io::p) out << pot  (i) <<" ";
    if(written & io::P) out << pex  (i) <<" ";
    if(written & io::a) out << acc  (i) <<" ";
    if(written & io::r) out << rho  (i) <<" ";
    if(written & io::y) out << aux  (i) <<" ";
    if(written & io::l) out << level(i) <<" ";
    if(written & io::k) out << key  (i) <<" ";
    if(written & io::f) out << flag (i) <<" ";
    if(written & io::n) out << num  (i) <<" ";
    if(i < N_sph()) {
      if(written & io::s) out << size (i) <<" ";
#ifdef falcON_SPH
      if(written & io::T) out << temp (i) <<" ";
      if(written & io::u) out << ein  (i) <<" ";
      if(written & io::S) out << ent  (i) <<" ";
      if(written & io::R) out << srho (i) <<" ";
      if(written & io::D) out << divv (i) <<" ";
      if(written & io::V) out << rotv (i) <<" ";
#endif
    }
    out<<"\n";
  }
  out.flush();
}
//------------------------------------------------------------------------------
void sbodies::write_yanc_binary(std::ostream&out, const io write) const
{
  if(write == io::o) return;
  register io written = writer(write,this);
  int w = written;
  out.write(CAST(char*,&w),sizeof(int));
  if(written&io::m) out.write(CAST(char*,MAS),N_bodies()*sizeof(real));
  if(written&io::x) out.write(CAST(char*,POS),N_bodies()*sizeof(vect));
  if(written&io::v) out.write(CAST(char*,VEL),N_bodies()*sizeof(vect));
  if(written&io::e) out.write(CAST(char*,EPS),N_bodies()*sizeof(real));
  if(written&io::p) out.write(CAST(char*,POT),N_bodies()*sizeof(real));
  if(written&io::P) out.write(CAST(char*,PEX),N_bodies()*sizeof(real));
  if(written&io::a) out.write(CAST(char*,ACC),N_bodies()*sizeof(vect));
  if(written&io::r) out.write(CAST(char*,RHO),N_bodies()*sizeof(real));
  if(written&io::y) out.write(CAST(char*,AUX),N_bodies()*sizeof(real));
  if(written&io::l) out.write(CAST(char*,LEV),N_bodies()*sizeof(indx));
  if(written&io::k) out.write(CAST(char*,KEY),N_bodies()*sizeof(int ));
  if(written&io::f) out.write(CAST(char*,FLG),N_bodies()*sizeof(int ));
  if(written&io::n) out.write(CAST(char*,NUM),N_bodies()*sizeof(uint));
  if(N_sph()) {
    if(written&io::s) out.write(CAST(char*,SIZ),N_sph   ()*sizeof(real));
#ifdef falcON_SPH
    if(written&io::T) out.write(CAST(char*,TEM),N_sph   ()*sizeof(real));
    if(written&io::u) out.write(CAST(char*,EIN),N_sph   ()*sizeof(real));
    if(written&io::S) out.write(CAST(char*,ENT),N_sph   ()*sizeof(real));
    if(written&io::R) out.write(CAST(char*,SRH),N_sph   ()*sizeof(real));
    if(written&io::D) out.write(CAST(char*,DVV),N_sph   ()*sizeof(real));
    if(written&io::V) out.write(CAST(char*,RTV),N_sph   ()*sizeof(vect));
#endif
  }
  out.flush();
}
//------------------------------------------------------------------------------
namespace nbdy {
  inline io reader(io            const&want,
		   io            const&present, 
		   const sbodies*const&B,
		   const char*   const&func)
  {
    register io read = io::o;
    for(register int i=0,bit=1; i!=io::N_quant(); ++i,bit<<=1)
      if(want & bit) {
	if     (! present & bit)
	  warning("[%s]: cannot read %s: not present in file",
		  func,io::S_quant(i));
	else if(!B->all_bits_contain(bit))
	  warning("[%s]: cannot read %s: no memory allocated",
		  func,io::S_quant(i));
	else read |= bit;
      }
    return read;
  }
}
//------------------------------------------------------------------------------
io sbodies::read_yanc_ascii(std::istream& in, const io want)
  // NOTE: reset(N,data) has to be called BEFORE.                               
  //       We assume the correct N and try to read all data wanted              
  //       Unlike previous versions, we cannot "enable" holding aux data        
{
  register int p; in>>p;
  register io pr=p, read=reader(want,pr,this,"sbodies::read_yanc_ascii()");
  reset(N_bodies(), read | all_bits(), N_sph());
  if((read & io::sphmax) && N_sph()==0)
    warning("[sbodies::read_yanc_ascii()]: cannot read SPH data");
  register real rdum;
  register vect vdum;
  for(register uint i=0; i!=N_bodies(); ++i) {
    if(pr & io::m) { if(read & io::m) in>>mass (i); else in>>rdum; }
    if(pr & io::x) { if(read & io::x) in>>pos  (i); else in>>vdum; }
    if(pr & io::v) { if(read & io::v) in>>vel  (i); else in>>vdum; }
    if(pr & io::e) { if(read & io::e) in>>eps  (i); else in>>rdum; }
    if(pr & io::p) { if(read & io::p) in>>pot  (i); else in>>rdum; }
    if(pr & io::P) { if(read & io::P) in>>pex  (i); else in>>rdum; }
    if(pr & io::a) { if(read & io::a) in>>acc  (i); else in>>vdum; }
    if(pr & io::r) { if(read & io::r) in>>rho  (i); else in>>rdum; }
    if(pr & io::y) { if(read & io::y) in>>aux  (i); else in>>rdum; }
    if(pr & io::l) { if(read & io::l) in>>level(i); else in>>rdum; }
    if(pr & io::k) { if(read & io::k) in>>key  (i); else in>>rdum; }
    if(pr & io::f) { if(read & io::f) in>>flg  (i); else in>>rdum; }
    if(pr & io::n) { if(read & io::n) in>>num  (i); else in>>rdum; }
    if(i<N_sph() && read & io::sphmax) {
      if(pr & io::s) { if(read & io::s) in>>size (i); else in>>rdum; }
#ifdef falcON_SPH
      if(pr & io::T) { if(read & io::T) in>>temp (i); else in>>rdum; }
      if(pr & io::u) { if(read & io::u) in>>ein  (i); else in>>rdum; }
      if(pr & io::S) { if(read & io::S) in>>ent  (i); else in>>rdum; }
      if(pr & io::R) { if(read & io::R) in>>srho (i); else in>>rdum; }
      if(pr & io::D) { if(read & io::D) in>>divv (i); else in>>rdum; }
      if(pr & io::V) { if(read & io::V) in>>rotv (i); else in>>rdum; }
#endif
    }
    SwallowRestofLine(in);
  }
  return read;
}
//------------------------------------------------------------------------------
#define READIT(IO,NAME,TYPE,NUM)				\
  if(pr&IO) {							\
    if(read&IO) in.read (CAST(char*,NAME),NUM*sizeof(TYPE));	\
    else        in.seekg(NUM*sizeof(TYPE),std::ios::cur);	\
  }
void sbodies::read_yanc_binary(std::istream&in, const io want)
{
  register int p; in>>p;
  register io pr=p, read=reader(want,pr,this,"sbodies::read_yanc_ascii()");
  reset(N_bodies(), read | all_bits(), N_sph());
  READIT(io::m,MAS,real,N_bodies())
  READIT(io::x,POS,vect,N_bodies())
  READIT(io::v,VEL,vect,N_bodies())
  READIT(io::e,EPS,real,N_bodies())
  READIT(io::p,POT,real,N_bodies())
  READIT(io::P,PEX,real,N_bodies())
  READIT(io::a,ACC,vect,N_bodies())
  READIT(io::r,RHO,real,N_bodies())
  READIT(io::y,AUX,real,N_bodies())
  READIT(io::l,LEV,indx,N_bodies())
  READIT(io::k,KEY,int ,N_bodies())
  READIT(io::f,FLG,int ,N_bodies())
  READIT(io::n,NUM,uint,N_bodies())
  if(N_sph()) {
    READIT(io::s,SIZ,real,N_sph())
#ifdef falcON_SPH
    READIT(io::T,TEM,real,N_sph())
    READIT(io::u,EIN,real,N_sph())
    READIT(io::S,ENT,real,N_sph())
    READIT(io::R,SRH,real,N_sph())
    READIT(io::D,DVV,real,N_sph())
    READIT(io::V,RTV,vect,N_sph())
#endif
  }
#undef READIT
}
//------------------------------------------------------------------------------
void sbodies::read_simple_ascii(std::istream&in, uint const&L, const io*item)
{
  register int K  =0;
  register io  get=io::o;
  while(item[K]!=io::o) get |= item[K++];
  reset(L,get|all_bits());                         // reset data                
  for(register int l=0; l!=L; ) {
    if(!in)
      falcON_ErrorF("end of input before data have been read",
		    "sbodies::read_simple_ascii()");
    if(! eat_line(in,'#') ) {
      for(register int k=0; k!=K; ++k)
	switch(item[k]) {
	case io::m: in >> mass (l); break;
	case io::x: in >> pos  (l); break;
	case io::v: in >> vel  (l); break;
	case io::e: in >> eps  (l); break;
	case io::p: in >> pot  (l); break;
	case io::P: in >> pex  (l); break;
	case io::a: in >> acc  (l); break;
	case io::r: in >> rho  (l); break;
	case io::y: in >> aux  (l); break;
	case io::l: in >> level(l); break;
	case io::k: in >> key  (l); break;
	case io::f: in >> flg  (l); break;
	case io::n: in >> num  (l); break;
	case io::s: in >> size (l); break;
	}
      ++l;
      SwallowRestofLine(in);
    }
  }
}
