// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// body.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2002                                          |
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

#ifdef ALLOW_NEMO                                  // compiler option           

#include  <public/nmio.h>                          // nbdy NEMO I/O support     
extern "C" {
# include <stdinc.h>                               // NEMO basics               
}

#endif

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
#ifdef ALLOW_NEMO
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
    NbdyErrorF("bodies set not open","BodySrce::read_nemo()")
  if(! is_supported(want&NEMOBITS()) )
    NbdyErrorF("data types not supported","BodySrce::read_nemo()")
  register int end=input.Number()+begin;
  if(end > NB)
    NbdyErrorF("too many bodies","BodySrce::read_nemo()")
  register io read = io::o;                        // data read                 
  // read the data                                                              
  // 1 positions and velocities are a bit more involved:                        
  // 1.1 both positions and velocities are required                             
  if(want & io::x && want & io::v) {               // IF([x,v] wanted)     >    
    if(input.is_present(nemo_io::pos) ||           //   IF(positions or         
       input.is_present(nemo_io::vel)) {           //      velocities to read) >
      if(input.is_present(nemo_io::pos)) {         //     IF(positions to read) 
	input.read(nemo_io::pos,CAST(real*,POS+begin)); //  read positions      
	read |= io::x;                             //       add x to read       
      } else if(warn)                              //     ELSE: warning         
	warning("[%s]: cannot read: x","BodySrce::read_nemo()");
      if(input.is_present(nemo_io::vel)) {         //     IF(velocities to read)
	input.read(nemo_io::vel,CAST(real*,VEL+begin)); //  read velocities     
	read |= io::v;                             //       add v to read       
      } else if(warn)                              //     ELSE: warning         
	warning("[%s]: cannot read: v","BodySrce::read_nemo()");
    } else if(input.is_present(nemo_io::posvel)) { //   < ELSE IF(x,v present) >
      input.read(nemo_io::posvel);                 //     read phases           
      for(register int j=0,i=begin; i!=end; ++i,++j) { // loop bodies         > 
	pos(i)[0] = input.bodies_phs(j)[0];
	vel(i)[0] = input.bodies_phs(j)[0+NDIM];
	pos(i)[1] = input.bodies_phs(j)[1];
	vel(i)[1] = input.bodies_phs(j)[1+NDIM];
#if NDIM > 2
	pos(i)[2] = input.bodies_phs(j)[2];
	vel(i)[2] = input.bodies_phs(j)[2+NDIM];
#endif
      }                                            //     <                     
      read |= io::xv;                              //     add phases to read    
    } else if(warn)                                //   < ELSE: warning         
      warning("[%s]: cannot read: (x,v)","BodySrce::read_nemo()");
  // 1.2 only positions are required                                            
  } else if(want & io::x) {                        // < ELSE IF( x wanted)     >
    if(input.is_present(nemo_io::pos)) {           //   IF(x present)          >
      input.read(nemo_io::pos,CAST(real*,POS+begin));  // read positions        
      read |= io::x;                               //     add poss to read      
    } else if(input.is_present(nemo_io::posvel)) { //   < ELSE IF(x,v present) >
      input.read(nemo_io::posvel);                 //     read phases           
      for(register int j=0,i=begin; i!=end; ++i,++j) { // loop bodies         > 
	pos(i)[0] = input.bodies_phs(j)[0];
	pos(i)[1] = input.bodies_phs(j)[1];
#if NDIM > 2
	pos(i)[2] = input.bodies_phs(j)[2];
#endif
      }                                            //     <                     
      read |= io::x;                               //     add poss to read      
    } else if(warn)                                //   < ELSE: warning         
      warning("[%s]: cannot read: x","BodySrce::read_nemo()");
  // 1.3 only velocities are required                                           
  } else if(want & io::v) {                        // < ELSE IF( v wanted)     >
    if(input.is_present(nemo_io::vel)) {           //   IF( v present)         >
      input.read(nemo_io::vel,CAST(real*,VEL+begin));  // read velocities       
      read |= io::v;                               //     add vels to read      
    } else if(input.is_present(nemo_io::posvel)) { //   < ELSE IF(x,v present) >
      input.read(nemo_io::posvel);                 //     read phases           
      for(register int j=0,i=begin; i!=end; ++i,++j) { // loop bodies         > 
	vel(i)[0] = input.bodies_phs(j)[0+NDIM];
	vel(i)[1] = input.bodies_phs(j)[1+NDIM];
#if NDIM > 2
	vel(i)[2] = input.bodies_phs(j)[2+NDIM];
#endif
      }                                            //     <                     
      read |= io::v;                               //     add vels to read      
    } else if(warn)                                //   < ELSE: warning         
      warning("[%s]: cannot read: v","BodySrce::read_nemo()");
  }                                                // <                         
  // 2. read remaining source data                                              
  if(want & io::f) {                               // IF(f wanted)             >
    if(input.is_present(nemo_io::flag)) {          //   IF(f present)          >
      input.read(nemo_io::flag,CAST(int*,FLG+begin));  // read FLG              
      read |= io::f;                               //     add f to read         
    } else if(warn)                                //   < ELSE: warning        <
      warning("[%s]: cannot read: f","BodySrce::read_nemo()");
  }
  if(want & io::m) {                               // IF(m wanted)             >
    if(input.is_present(nemo_io::mass)) {          //   IF(m present)          >
      input.read(nemo_io::mass, MAS+begin);        //     read MAS              
      read |= io::m;                               //     add m to read         
    } else if(warn)                                //   < ELSE: warning        <
      warning("[%s]: cannot read: m","BodySrce::read_nemo()");
  }
  if(want & io::e) {                               // IF(e wanted)             >
    if(input.is_present(nemo_io::eps)) {           //   IF(e present)          >
      input.read(nemo_io::eps, EPS+begin);         //     read EPS              
      read |= io::e;                               //     add e to read         
    } else if(warn)                                //   < ELSE: warning        <
      warning("[%s]: cannot read: e","BodySrce::read_nemo()");
  }
  if(want & io::k) {                               // IF(k wanted)             >
    if(input.is_present(nemo_io::key)) {           //   IF(k present)          >
      input.read(nemo_io::key, KEY+begin);         //     read KEY              
      read |= io::k;                               //     add k to read         
    } else if(warn)                                //   < ELSE: warning        <
      warning("[%s]: cannot read: k","BodySrce::read_nemo()");
  }
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
    NbdyErrorF("bodies set not open","BodySink::read_nemo()")
  if(! is_supported(want&NEMOBITS()))
    NbdyErrorF("data types not supported","BodySink::read_nemo()")
  register int end=input.Number()+begin;
  if(end > NB)
    NbdyErrorF("too many bodies","BodySink::read_nemo()")
  register io read = io::o;                        // data read                 
  if(want & io::a) {                               // IF(a wanted)             >
    if(input.is_present(nemo_io::acc)) {           //   IF(a present)          >
      input.read(nemo_io::acc,CAST(real*,ACC+begin)); //  read ACC              
      read |= io::a;                               //     add a to read         
    } else if(warn)                                //   < ELSE: warning        <
      warning("[%s]: cannot read: a","BodySink::read_nemo()");
  }
  if(want & io::p) {                               // IF(p wanted)             >
    if(input.is_present(nemo_io::pot)) {           //   IF(p present)          >
      input.read(nemo_io::pot, POT+begin);         //     read POT              
      read |= io::p;                               //     add p to read         
    } else if(warn)                                //   < ELSE: warning        <
      warning("[%s]: cannot read: p","BodySink::read_nemo()");
  }
  if(want & io::r) {                               // IF(r wanted)             >
    if(input.is_present(nemo_io::rho)) {           //   IF(r present)          >
      input.read(nemo_io::rho, RHO+begin);         //     read RHO              
      read |= io::r;                               //     add r to read         
    } else if(warn)                                //   < ELSE: warning        <
      warning("[%s]: cannot read: r","BodySink::read_nemo()");
  }
  if(want & io::y) {                               // IF(y wanted)             >
    if(input.is_present(nemo_io::aux)) {           //   IF(y present)          >
      input.read(nemo_io::aux, AUX+begin);         //     read AUX              
      read |= io::y;                               //     add y to read         
    } else if(warn)                                //   < ELSE: warning        <
      warning("[%s]: cannot read: y","BodySink::read_nemo()");
  }
  if(want & io::l) {                               // IF(l wanted)             >
    if(input.is_present(nemo_io::level)) {         //   IF(l present)          >
      input.read(nemo_io::level,CAST(short*,LEV+begin)); // read LEV            
      read |= io::l;                               //     add l to read         
    } else if(warn)                                //   < ELSE: warning        <
      warning("[%s]: cannot read: l","BodySink::read_nemo()");
  }
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
  if(input.is_present(nemo_io::time))              //   IF time available      >
    input.read(nemo_io::time,&t);                  //     read time             
  else {                                           //   < ELSE >                
    t = zero;                                      //     default t=0           
    if(time && warn)                               //     IF time wanted: warn  
      warning("[%s]: no time found in snapshot; defaulting to 0",
	      "sbodies::read_nemo_particles()");
  };                                               //   <                       
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
    NbdyErrorF("no snapshot found","sbodies::read_nemo_snapshot()") // ERROR    
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
  for(register int i=0; i!=Nin; ++i) {             // loop nemo inputs         >
    in[i].open_set(nemo_io::param);                //   open nemo parameter set 
    NN += Nnum[i] = in[i].read_N();                //   read & accumulate N     
    if(in[i].is_present(nemo_io::time)) {          //   IF time available      >
      in[i].read(nemo_io::time,&ti);               //     read time             
      if(time_read && ti != t)                     //     IF(time differs)      
	NbdyErrorF("different times found in snapshots",//  issue ERROR         
		   "sbodies::read_nemo_particles()")
      t = ti;                                      //     set time              
      time_read = true;                            //     set flag: time read   
    }                                              //   <                       
    in[i].close_set(nemo_io::param);               //   close nemo parameter set
    if(want) get |= want[i];
  }                                                // <                         
  if(time) {                                       // IF(time wanted)          >
    *time = t;                                     //   assign to output value  
    if(!time_read)                                 //   IF(time not read)       
      warning("[%s]: no time found in snapshots; defaulting to 0",
	      "sbodies::read_nemo_particles()");
  }                                                // <                         
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
  for(register int i=0; i!=Nin; ++i) {             // loop nemo inputs         >
    in[i].open_set(nemo_io::bodies);               //   open nemo particle set  
    if(want) wanted=want[i]&NEMOBITS();            //   data to be read         
    else     wanted=io::mxv;
    read[i] = BodySrce::read_nemo(in[i],wanted,warn,begin); // read source data 
    read[i]|= BodySink::read_nemo(in[i],wanted,warn,begin); // read sink data   
    begin  += in[i].Number();                      //   increment begin         
    in[i].close_set(nemo_io::bodies);              //   close nemo particle set 
    in[i].reset();                                 //   reset nemo input        
  }                                                // <                         
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
  for(register int i=0; i!=Nin; ++i) {             // loop nemo inputs         >
    if(!in[i].is_present(nemo_io::snap))           //   IF no snapshot present  
      NbdyErrorF("no snapshot found","sbodies::read_nemo_snapshot()") // ERROR  
    in[i].open_set(nemo_io::snap);                 // open nemo snapshot        
  }
  bool tt=read_nemo_particles(in,Nin,read,Nnum,time, // read particles if       
			      want,times,warn);    //   time in range           
  for(register int i=0; i!=Nin; ++i)               // loop nemo inputs         >
    in[i].close_set(nemo_io::snap);                //   open nemo snapshot     <
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
  if(write & io::m)                                // IF(masses to be written) >
    output.write(nemo_io::mass,MAS+begin);         //   write out masses       <
  // 2.2 - 2.3 output of positions & velocities                                 
#if(1)                  // separately as positions x_i and velocities v_i       
  // 2.2 output positions                                                       
  if(write & io::x)                                // IF poss for output     >  
    output.write(nemo_io::pos,CAST(real*,POS+begin));  // write out poss     <  
  // 2.3 output velocities                                                      
  if(write & io::v)                                // IF vels for output     >  
    output.write(nemo_io::vel,CAST(real*,VEL+begin));  // write out vels     <  
#else                   // as phase-space coordinates [x,v]_i                   
  // 2.2.1 both positions and velocities to be written out                      
  if(write & io::x && write & io::v) {             // IF([x,v] to be written)  >
    output.allocarrays();                          //   allocate memory for x,v 

    for(j=0,i=begin; j!=Nout; ++i,++j) {           //   loop bodies           > 
      output.bodies_phs(j)[0]      = pos(i)[0];
      output.bodies_phs(j)[0+NDIM] = vel(i)[0];
      output.bodies_phs(j)[1]      = pos(i)[1];
      output.bodies_phs(j)[1+NDIM] = vel(i)[1];
#if NDIM > 2
      output.bodies_phs(j)[2]      = pos(i)[2];
      output.bodies_phs(j)[2+NDIM] = vel(i)[2];
#endif
    }                                              //   <                       
    output.write(nemo_io::posvel);                 //   write out phases        
  // 2.2.2 only positions to be written out                                     
  } else if(write & io::x) {                       // <ELSE IF(x to be written)>
    output.allocarrays();                          //   allocate memory for x,v 
    for(j=0,i=begin; j!=Nout; ++i,++j) {           //   loop bodies           > 
      output.bodies_phs(j)[0]      = pos(i)[0];
      output.bodies_phs(j)[1]      = pos(i)[1];
#if NDIM > 2
      output.bodies_phs(j)[2]      = pos(i)[2];
#endif
    }                                              //   <                       
    output.write(nemo_io::posvel);                 //   write out phases        
  // 2.2.3 only velocities to be written out                                    
  } else if(write & io::v) {                       // <ELSE IF(v to be written)>
    output.allocarrays();                          //   allocate memory for x,v 
    for(j=0,i=begin; j!=Nout; ++i,++j) {           //   loop bodies           > 
      output.bodies_phs(j)[0+NDIM] = vel(i)[0];
      output.bodies_phs(j)[1+NDIM] = vel(i)[1];
#if NDIM > 2
      output.bodies_phs(j)[2+NDIM] = vel(i)[2];
#endif
    }                                              //   <                       
    output.write(nemo_io::posvel);                 //   write out phases        
  }                                                // <                         
#endif
  // 2.4 output accelerations                                                   
  if(write & io::a)                                // IF accs for output     >  
    output.write(nemo_io::acc,CAST(real*,ACC+begin));  // write out accs     <  
  // 2.5 output potential(s)                                                    
  if(write & io::p)                                // IF pot for output   >     
    output.write(nemo_io::pot,POT+begin);          //   write out potential     
  // 2.6 output density                                                         
  if(write & io::r && has_rho())                   // IF(rhos for output)    >  
    output.write(nemo_io::rho,RHO+begin);          //     write out scalars  <  
  // 2.7 output aux data                                                        
  if(write & io::y && has_aux())                   // IF(auxx for output)    >  
    output.write(nemo_io::aux,AUX+begin);          //   write out scalars    <  
  // 2.8 output eps_i                                                           
  if(write & io::e && has_eps())                   // IF(eps_ifor output)    >  
    output.write(nemo_io::eps,EPS+begin);          //   write out potential  <  
  // 2.9 output level data                                                      
  if(write & io::l && has_lev())                   // IF(level for output)   >  
    output.write(nemo_io::level,CAST(short*,LEV+begin)); // write out shorts <  
  // 2.10 output flag data                                                      
  if(write & io::f)                                // IF(flag for output)    >  
    output.write(nemo_io::flag,CAST(int*,FLG+begin)); // write out flags     <  
  // 2.11 output key data                                                       
  if(write & io::k && has_key())                   // IF(key for output)     >  
    output.write(nemo_io::key,KEY+begin);          //   write out keys       <  
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
// class sbodies: YANC I/O                                                      
////////////////////////////////////////////////////////////////////////////////
void sbodies::write_yanc_ascii(std::ostream &out, const io write) const
{
  if(write == io::o) return;
  register io   written = io::o;
  if(write & io::m)                written |= io::m;
  if(write & io::x)                written |= io::x;
  if(write & io::v)                written |= io::v;
  if(write & io::e && has_eps())   written |= io::e;
  if(write & io::p)                written |= io::p;
  if(write & io::a)                written |= io::a;
  if(write & io::r && has_rho())   written |= io::r;
  if(write & io::y && has_aux())   written |= io::y;
  if(write & io::l && has_level()) written |= io::l;
  if(write & io::k && has_key())   written |= io::k;
  if(write & io::f && has_flg())   written |= io::f;
  if(write & io::n && has_num())   written |= io::n;
  if(write & io::s && has_size())  written |= io::s;
  out<<written<<std::endl;
  for(register uint i=0; i!=N_bodies(); ++i) {
    if(written & io::m) out << mass (i) <<" ";
    if(written & io::x) out << pos  (i) <<" ";
    if(written & io::v) out << vel  (i) <<" ";
    if(written & io::e) out << eps  (i) <<" ";
    if(written & io::p) out << pot  (i) <<" ";
    if(written & io::a) out << acc  (i) <<" ";
    if(written & io::r) out << rho  (i) <<" ";
    if(written & io::y) out << aux  (i) <<" ";
    if(written & io::l) out << level(i) <<" ";
    if(written & io::k) out << key  (i) <<" ";
    if(written & io::f) out << flag (i) <<" ";
    if(written & io::n) out << num  (i) <<" ";
    if(written & io::s) out << size (i) <<" ";
    out<<"\n";
  }
  out.flush();
}
//------------------------------------------------------------------------------
io sbodies::read_yanc_ascii(std::istream& in, const io want)
  // NOTE: reset(N,data) has to be called BEFORE.                               
  //       We assume the correct N and try to read all data wanted              
  //       Unlike previous versions, we cannot "enable" holding aux data        
{
  register io pr,read=io::o;
  { register int p; in >> p; pr=p; }
  if(want&io::m) { if(pr&io::m) read|=io::m;
                   else warning("[sbodies::read_yanc_ascii()]: cannot read m");}
  if(want&io::x) { if(pr&io::x) read|=io::x;
                   else warning("[sbodies::read_yanc_ascii()]: cannot read x");}
  if(want&io::v) { if(pr&io::v) read|=io::v;
                   else warning("[sbodies::read_yanc_ascii()]: cannot read v");}
  if(want&io::e) { if(pr&io::e) read|=io::e;
                   else warning("[sbodies::read_yanc_ascii()]: cannot read e");}
  if(want&io::p) { if(pr&io::p) read|=io::p;
                   else warning("[sbodies::read_yanc_ascii()]: cannot read p");}
  if(want&io::a) { if(pr&io::a) read|=io::a;
                   else warning("[sbodies::read_yanc_ascii()]: cannot read a");}
  if(want&io::r) { if(pr&io::r) read|=io::r;
                   else warning("[sbodies::read_yanc_ascii()]: cannot read r");}
  if(want&io::y) { if(pr&io::y) read|=io::y;
                   else warning("[sbodies::read_yanc_ascii()]: cannot read y");}
  if(want&io::l) { if(pr&io::l) read|=io::l;
                   else warning("[sbodies::read_yanc_ascii()]: cannot read l");}
  if(want&io::k) { if(pr&io::k) read|=io::k;
                   else warning("[sbodies::read_yanc_ascii()]: cannot read k");}
  if(want&io::f) { if(pr&io::f) read|=io::f;
                   else warning("[sbodies::read_yanc_ascii()]: cannot read f");}
  if(want&io::n) { if(pr&io::n) read|=io::n;
                   else warning("[sbodies::read_yanc_ascii()]: cannot read n");}
  if(want&io::s) { if(pr&io::s) read|=io::s;
                   else warning("[sbodies::read_yanc_ascii()]: cannot read s");}
  reset(N_bodies(), read | all_bits());
  register real rdum;
  register vect vdum;
  for(register uint i=0; i!=N_bodies(); ++i) {
    if(pr & io::m) { if(read & io::m) in>>mas(i); else in>>rdum; }
    if(pr & io::x) { if(read & io::x) in>>pos(i); else in>>vdum; }
    if(pr & io::v) { if(read & io::v) in>>vel(i); else in>>vdum; }
    if(pr & io::e) { if(read & io::e) in>>eps(i); else in>>rdum; }
    if(pr & io::e) { if(read & io::e) in>>pot(i); else in>>rdum; }
    if(pr & io::a) { if(read & io::a) in>>acc(i); else in>>vdum; }
    if(pr & io::r) { if(read & io::r) in>>rho(i); else in>>rdum; }
    if(pr & io::y) { if(read & io::y) in>>aux(i); else in>>rdum; }
    if(pr & io::l) { if(read & io::l) in>>lev(i); else in>>rdum; }
    if(pr & io::k) { if(read & io::k) in>>key(i); else in>>rdum; }
    if(pr & io::f) { if(read & io::f) in>>flg(i); else in>>rdum; }
    if(pr & io::n) { if(read & io::n) in>>num(i); else in>>rdum; }
    if(pr & io::s) { if(read & io::s) in>>siz(i); else in>>rdum; }
    SwallowRestofLine(in);
  }
  return read;
}
//------------------------------------------------------------------------------
void sbodies::write_yanc_binary(std::ostream&out, const io write) const
{
  if(write == io::o) return;
  register io written = io::o;
  if(write & io::m)                written |= io::m;
  if(write & io::x)                written |= io::x;
  if(write & io::v)                written |= io::v;
  if(write & io::e && has_eps())   written |= io::e;
  if(write & io::p)                written |= io::p;
  if(write & io::a)                written |= io::a;
  if(write & io::r && has_rho())   written |= io::r;
  if(write & io::y && has_aux())   written |= io::y;
  if(write & io::l && has_level()) written |= io::l;
  if(write & io::k && has_key())   written |= io::k;
  if(write & io::f)                written |= io::f;
  if(write & io::n && has_num())   written |= io::n;
  if(write & io::s && has_size())  written |= io::s;
  int w = written;
  out.write(CAST(char*,&w),sizeof(int));
  if(written&io::m) out.write(CAST(char*,MAS),N_bodies()*sizeof(real));
  if(written&io::x) out.write(CAST(char*,POS),N_bodies()*sizeof(vect));
  if(written&io::v) out.write(CAST(char*,VEL),N_bodies()*sizeof(vect));
  if(written&io::e) out.write(CAST(char*,EPS),N_bodies()*sizeof(real));
  if(written&io::p) out.write(CAST(char*,POT),N_bodies()*sizeof(real));
  if(written&io::a) out.write(CAST(char*,ACC),N_bodies()*sizeof(vect));
  if(written&io::r) out.write(CAST(char*,RHO),N_bodies()*sizeof(real));
  if(written&io::y) out.write(CAST(char*,AUX),N_bodies()*sizeof(real));
  if(written&io::l) out.write(CAST(char*,LEV),N_bodies()*sizeof(indx));
  if(written&io::k) out.write(CAST(char*,KEY),N_bodies()*sizeof(int ));
  if(written&io::f) out.write(CAST(char*,FLG),N_bodies()*sizeof(int ));
  if(written&io::n) out.write(CAST(char*,NUM),N_bodies()*sizeof(uint));
  if(written&io::s) out.write(CAST(char*,SIZ),N_bodies()*sizeof(real));
  out.flush();
}
//------------------------------------------------------------------------------
#define READIT(IO,NAME,TYPE)						\
  if(pr&IO) { 								\
    if(read&IO) in.read (CAST(char*,NAME),N_bodies()*sizeof(TYPE));	\
    else        in.seekg(N_bodies()*sizeof(TYPE),std::ios::cur);	\
  }
void sbodies::read_yanc_binary(std::istream&in, const io want)
{
  register io pr,read=io::o;
  { register int w; in.read(CAST(char*,&w),sizeof(int)); pr=w; }
  if(want&io::m) { if(pr&io::m) read|=io::m;
                   else warning("[%s]: cannot read: m",
				"sbodies::read_yanc_binary()");}
  if(want&io::x) { if(pr&io::x) read|=io::x;
                   else warning("[%s]: cannot read: x",
				"sbodies::read_yanc_binary()");}
  if(want&io::v) { if(pr&io::v) read|=io::v;
                   else warning("[%s]: cannot read: v",
				"sbodies::read_yanc_binary()");}
  if(want&io::e) { if(pr&io::e) read|=io::e;
                   else warning("[%s]: cannot read: e",
				"sbodies::read_yanc_binary()");}
  if(want&io::p) { if(pr&io::p) read|=io::p;
                   else warning("[%s]: cannot read: p",
				"sbodies::read_yanc_binary()");}
  if(want&io::a) { if(pr&io::a) read|=io::a;
                   else warning("[%s]: cannot read: a",
				"sbodies::read_yanc_binary()");}
  if(want&io::r) { if(pr&io::r) read|=io::r;
                   else warning("[%s]: cannot read: r",
				"sbodies::read_yanc_binary()");}
  if(want&io::y) { if(pr&io::y) read|=io::y;
                   else warning("[%s]: cannot read: y",
				"sbodies::read_yanc_binary()");}
  if(want&io::l) { if(pr&io::l) read|=io::l;
                   else warning("[%s]: cannot read: l",
				"sbodies::read_yanc_binary()");}
  if(want&io::k) { if(pr&io::k) read|=io::k;
                   else warning("[%s]: cannot read: k",
				"sbodies::read_yanc_binary()");}
  if(want&io::f) { if(pr&io::f) read|=io::f;
                   else warning("[%s]: cannot read: f",
				"sbodies::read_yanc_binary()");}
  if(want&io::n) { if(pr&io::n) read|=io::n;
                   else warning("[%s]: cannot read: n",
				"sbodies::read_yanc_binary()");}
  if(want&io::s) { if(pr&io::s) read|=io::s;
                   else warning("[%s]: cannot read: s",
				"sbodies::read_yanc_binary()");}
  reset(N_bodies(), read | all_bits());
  READIT(io::m,MAS,real)
  READIT(io::x,POS,vect)
  READIT(io::v,VEL,vect)
  READIT(io::e,EPS,real)
  READIT(io::p,POT,real)
  READIT(io::a,ACC,vect)
  READIT(io::r,RHO,real)
  READIT(io::y,AUX,real)
  READIT(io::l,LEV,indx)
  READIT(io::k,KEY,int )
  READIT(io::f,FLG,int )
  READIT(io::n,NUM,uint)
  READIT(io::s,SIZ,real)
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
      NbdyErrorF("end of input before data have been read",
		 "sbodies::read_simple_ascii()")
    if(! eat_line(in,'#') )
      for(register int k=0; k!=K; ++k)
	switch(item[k]) {
	case io::m: in >> mas(l++); break;
	case io::x: in >> pos(l++); break;
	case io::v: in >> vel(l++); break;
	case io::e: in >> eps(l++); break;
	case io::p: in >> pot(l++); break;
	case io::a: in >> acc(l++); break;
	case io::r: in >> rho(l++); break;
	case io::y: in >> aux(l++); break;
	case io::l: in >> lev(l++); break;
	case io::k: in >> key(l++); break;
	case io::f: in >> flg(l++); break;
	case io::n: in >> num(l++); break;
	case io::s: in >> siz(l++); break;
	}
  }
}
