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
using namespace nbdy;
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class nbdy::bodies                                                           
//                                                                              
////////////////////////////////////////////////////////////////////////////////
//                                                                              
//   1. Data Management                                                         
//                                                                              
//------------------------------------------------------------------------------
// construction: set bits, N, allocate DATA[] if N!=0                           
bodies::bodies(uint nb,                            //[I: # bodies]              
	       io   bits,                          //[I: data bits]             
	       uint ns,                            //[I: # SPH particles]       
	       bool cu) :                          //[I: CUSE flag]             
  bodies_data<0> ( max(nb,ns), ns, cu )            // set base         	        
{
  for(int i=1,b=0; b!=IO_NQUANT; ++b,i<<=1) if(bits&i)
    add_data_void(b, Number(b)? new char[Number(b) * IO_ZQUANT[b]] : 0);
}
//------------------------------------------------------------------------------
// construction: make a partial copy, only copying data specified by 2nd arg    
bodies::bodies(bodies const&B,                     // I: bodies                 
	       io           copy) :                //[I: which data to copy]    
  bodies_data<0> ( B.N_bodies(), 
		   B.N_sph(),
		   B.changes_in_tree_usage_flags() )
{
  if(! B.srce_data_changed() ) mark_srce_data_read();
  if(! B.sph_data_changed () ) mark_sph_data_read ();
  const io bits = B.all_bits() & copy;
  if(bits == io::o || N_bodies() == 0u) return;
  for(int i=1,b=0; b!=IO_NQUANT; ++b,i<<=1) if(bits&i) {
    add_data_void(b, Number(b)? new char[Number(b) * IO_ZQUANT[b]] : 0);
    if(Number(b))
      memcpy(data_void(b), B.data_void(b), Number(b)*IO_ZQUANT[b]);
  }
}
//------------------------------------------------------------------------------
// copy body data indicated. If NBOD, NSPH deviating, adjust first              
void bodies::copy(bodies const&B,                  // I: bodies                 
		  io           copy)               //[I: which data to copy]    
{
  if(B.changes_in_tree_usage_flags()) mark_tree_usage_change();
  else                                after_tree_growth();
  const io bits = B.all_bits() & copy;
  resetN(B.N_bodies(),B.N_sph());
  add_fields(bits);
  if(bits == io::o || N_bodies() == 0u) return;
  for(int i=1,b=0; b!=IO_NQUANT; ++b,i<<=1)
    if(bits&i && Number(b))
      memcpy(data_void(b), B.data_void(b), Number(b)*IO_ZQUANT[b]);
  if(bits & io::source) mark_srce_data_changed();
  if(bits & io::sphmax) mark_sph_data_changed ();
}
//------------------------------------------------------------------------------
// destruction: delete all supported data fields                                
bodies::~bodies()
{
  for(int i=1,b=0; b!=IO_NQUANT; ++b,i<<=1)
    if(has(i) && data_void(b))
      delete[] static_cast<char*>(data_void(b));
}
//------------------------------------------------------------------------------
// adds unsupported DATA[] fields as indicated by arg                           
// supported DATA[] fields are not affected                                     
void bodies::add_fields(io bits)                   // I: which fields to add    
{
  bits &= ~all_bits();
  for(int i=1,b=0; b!=IO_NQUANT; ++b,i<<=1) if(bits&i)
    add_data_void(b, Number(b)? new char[Number(b) * IO_ZQUANT[b]] : 0);
}
//------------------------------------------------------------------------------
// adds unsupported DATA[] field as indicated by arg                            
// supported DATA[] fields are not affected                                     
void bodies::add_field(int b)                      // I: which field to add     
{
  if(! has(1<<b) ) 
    add_data_void(b, Number(b)? new char[Number(b) * IO_ZQUANT[b]] : 0);
}
//------------------------------------------------------------------------------
// removes supported DATA[] fields as indicated by arg                          
void bodies::del_fields(io bits)                   // I: which fields to delete 
{
  bits &= all_bits();
  for(int i=1,b=0; b!=IO_NQUANT; ++b,i<<=1) if(bits&i) {
    if(data_void(b)) delete[] static_cast<char*>(data_void(b));
    del_data_void(b);
  }
}
//------------------------------------------------------------------------------
// removes supported DATA[] field as indicated by arg                           
void bodies::del_field(int b)                      // I: which fields to delete 
{
  if(has(1<<b)) del_data_void(b);
}
//------------------------------------------------------------------------------
// resets DATA[] fields to arg                                                  
// DATA[] fields already supported AND contained in arg are not affected.       
void bodies::reset_fields(io bits)                 // I: which fields to support
{
  io newbits =  bits & ~all_bits();
  io delbits = ~bits &  all_bits();
  for(int i=1,b=0; b!=IO_NQUANT; ++b,i<<=1)
    if     (newbits&i)
      add_data_void(b, Number(b)? new char[Number(b) * IO_ZQUANT[b]] : 0);
    else if(delbits&i) {
      if(data_void(b)) delete[] static_cast<char*>(data_void(b));
      del_data_void(b);
    }
}
//------------------------------------------------------------------------------
// reset NBODIES & NSPH. If different from old, previous data are deleted       
void bodies::resetN(uint nb,                       // I: # bodies               
		    uint ns)                       //[I: # SPH particles]       
{
  for(int i=1,b=0; b!=IO_NQUANT; ++b,i<<=1) {
    uint N = b < IO_NOSPH? nb : ns;
    if(has(i) && N!=Number(b)) {
      if(data_void(b)) delete[] static_cast<char*>(data_void(b));
      add_data_void(b, N? new char[N * IO_ZQUANT[b]] : 0);
    }
  }
  if(has(io::source) && nb!=N_bodies()) mark_srce_data_changed();
  if(has(io::sphmax) && ns!=N_sph()   ) mark_sph_data_changed();
  setN(nb,ns);
}
//------------------------------------------------------------------------------
// change NBODIES & NSPH. Old data are preserved for n <= min(N_new,N_old)      
void bodies::changeN(uint nb,                      // I: # bodies               
		     uint ns)                      //[I: # SPH particles]       
{
  for(int i=1,b=0; b!=IO_NQUANT; ++b,i<<=1) {
    uint N = b < IO_NOSPH? nb : ns;
    if(has(i) && N!=Number(b)) {
      char *NEW = N? new char[N * IO_ZQUANT[b]] : 0;
      if(data_void(b)) {
	uint copy = min(N, Number(b));
	if(copy) memcpy(NEW,data_void(b),copy*IO_ZQUANT[b]);
	delete[] static_cast<char*>(data_void(b));
      }
      add_data_void(b,NEW);
    }
  }
  if(has(io::source) && nb!=N_bodies()) mark_srce_data_changed();
  if(has(io::sphmax) && ns!=N_sph()   ) mark_sph_data_changed();
  setN(nb,ns);
}
//------------------------------------------------------------------------------
//                                                                              
//   2. NEMO Data I/O                                                           
//                                                                              
//------------------------------------------------------------------------------
#ifdef falcON_NEMO
namespace { using namespace nbdy; using nbdy::uint;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // struct nemo_body_io<>                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<int IO> struct nemo_body_io : public field_io<IO,0> {
    static void read_nemo (nemo_in  const&I, bodies*const&B, io&R,
			   size_t F=0)
    {
      if(is_nemo && is_present(I)) {
	B->add_field(bit);
	read_nemo_array (I,B->template data_io<IO>() + F ,R);
      }
    }
    static void write_nemo(nemo_out const&O, const bodies*const&B, 
			   size_t F=0) {
      if(is_nemo && B->has(IO))
	write_nemo_array(O,B->template data_io<IO>() + F);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // struct nemo_read_write                                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  struct nemo_read_write {
    static const nemo_in *input;
    static const nemo_out*output;
    static io             getput;
    static size_t         first;
    static io            *read;
    static void set_read (const nemo_in&i, io g, io&r, size_t f=0) {
      input  = &i;
      getput =  g;
      first  =  f; 
      read   = &r;
    }
    static void set_write(const nemo_out&o, io p, size_t f=0) {
      output = &o;
      getput =  p;
      first  =  f;
    }
  };
  const nemo_in *  nemo_read_write::input;
  const nemo_out*  nemo_read_write::output;
  io               nemo_read_write::getput;
  size_t           nemo_read_write::first;
  io            *  nemo_read_write::read;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nemo_reader<>                                                      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<int IO> struct nemo_reader : 
    protected nemo_body_io<IO>,
    protected nemo_read_write {
    static void act_on_bodies(bodies*const&B) {
      if(is_nemo && getput & IO)
	read_nemo(*input,B,*read,first);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nemo_writer<>                                                      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<int IO> struct nemo_writer : 
    protected nemo_body_io<IO>,
    protected nemo_read_write {
    static void act_on_bodies(const bodies*const&B) { 
      if(is_nemo && getput & IO) write_nemo(*output,B,first);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
}
//------------------------------------------------------------------------------
// read a single nemo particle set into bodies begin to begin+N                 
io bodies::read_nemo(                              // R: data read              
		     nemo_in const&input,          // I: nemo input             
		     io            get,            // I: what to read           
		     size_t        begin)          //[I: first body to get]     
{
  io read = io::o;
  if(get&io::xv &&                                 // IF x,v wanted AND         
     input.is_present(nemo_io::posvel)) {          //    x,v given as phases    
    input.read(nemo_io::posvel);                   //   read phases             
    if(get & io::x) {                              //   IF x wanted             
      add_field(field_io<io::x,0>::bit);           //     allocate positions    
      for(int i=0; i!=N_bodies(); ++i)             //     LOOP bodies           
	pos(i+begin).copy(input.bodies_phs(i));    //       copy phases -> x    
      read |= io::x;                               //     update read           
    }                                              //   ENDIF                   
    if(get & io::v) {                              //   IF v wanted             
      add_field(field_io<io::v,0>::bit);           //     allocate velocities   
      for(int i=0; i!=N_bodies(); ++i)             //     LOOP bodies           
	vel(i+begin).copy(input.bodies_phs(i)+Ndim); //     copy phases -> v    
      read |= io::v;                               //     update read           
    }                                              //   ENDIF                   
  }                                                // ENDIF                     
  nemo_read_write::set_read(input,get& ~read,read,begin);
  LoopIO<nemo_reader>::loop(this);                 // read all data             
  return read;                                     // return read               
}
//------------------------------------------------------------------------------
// read a single nemo snapshot, supposed to be open                             
bool                                               // R: was time in range?     
bodies::read_nemo_particles(nemo_in const&input,   // I: nemo input             
			    io           &read,    // O: what has been read     
			    double       *time,    // O: time                   
			    io            want,    // I: what to read           
			    char         *times,   //[I: time range]            
			    bool          warn)    //[I: warn: missing data]    
{
  double t;                                        // time                      
  char   miss[IO_NQUANT+1];                        // missing data              
  // 1. read parameter set: read time, check for time in range, read NBOD, NSPH 
  input.open_set(nemo_io::param);                  // open nemo parameter set   
  if(input.is_present(nemo_io::time))              //   IF time available       
    input.read(nemo_io::time,&t);                  //     read time             
  else {                                           //   ELSE                    
    t = 0.;                                        //     default t=0           
    if(time && warn)                               //     IF time wanted: warn  
      warning("bodies::read_nemo: no time found in snapshot; defaulting to 0");
  };                                               //   ENDIF                   
  if(time) *time = t;                              //   IF wanted, set time     
  if(! time_in_range(t,times)) {                   //   IF time !in range       
    input.reset();                                 //     reset nemo input      
    input.close_set(nemo_io::param);               //     close parameter set   
    return false;                                  //     RETURN                
  }                                                //   ENDIF                   
  input.read_N();                                  //   read NBOD, NSPH         
  input.close_set(nemo_io::param);                 // close nemo parameter set  
  // 2. check for data support & prepare for data reading                       
  io get = want & NIOBits;                         // bits act'ly NEMO readable 
  if(warn && want != get)                          // IF wanted data not readabl
    warning("bodies::read_nemo: %s desired but not NEMO readable",
	    want.missing(get).make_word(miss));    //   WARNING                 
  if(get == io::o) return true;                    // IF nothing wanted: DOME   
  resetN(input.Number(),input.NumberSPH());        // adjust N if necessary     
  // 3. open bodies set & read data                                             
  input.open_set(nemo_io::bodies);                 // open nemo particle set    
  read = read_nemo(input,get);                     //   read data               
  input.close_set(nemo_io::bodies);                // close nemo particle set   
  // 4. clean up & close                                                        
  input.reset();                                   // reset nemo input          
  if(read & io::source) mark_srce_data_changed();
  if(read & io::sphmax) mark_sph_data_changed ();
  if(warn && get != read) warning("bodies::read_nemo: couldn't read %s",
				  get.missing(read).make_word(miss));
  return true;                                     // return                    
}
//------------------------------------------------------------------------------
// read a single particle set combining several nemo input streams              
bool
bodies::read_nemo_particles(                       // R: was time in range?     
			    const nemo_in*in,      // I: nemo inputs            
			    int     const&Nin,     // I: # nemo inputs          
			    io           *read,    // O: what has been read     
			    uint         *Nnum,    // O: how many have been read
			    double       *time,    //[O: time]                  
			    const io     *want,    //[I: what to read]          
			    char         *times,   //[I: time range]            
			    bool          warn)    //[I: warn: missing data]    
{
  bool   time_read=false;                          // time read already?        
  double t=0,ti;                                   // time                      
  uint   NN=0u,NNS=0u;                             // total number to read      
  char   miss[IO_NQUANT+1];                        // missing data              
  io     wanted = want? io::o : io::mxv;           // what do we want?          
  // 1. open parameter sets; read N & time; close parameter sets;               
  for(int i=0; i!=Nin; ++i) {                      // LOOP nemo inputs          
    in[i].open_set(nemo_io::param);                //   open nemo parameter set 
    in[i].read_N();                                //   read N & N_sph          
    NN += Nnum[i] = in[i].Number();                //   accumulate N            
    NNS+= in[i].NumberSPH();                       //   accumulate N_sph        
    if(in[i].is_present(nemo_io::time)) {          //   IF time available       
      in[i].read(nemo_io::time,&ti);               //     read time             
      if(time_read && ti != t)                     //     IF time differs       
	falcON_ErrorF("different times found in snapshots",//  issue ERROR      
		      "bodies::read_nemo_particles()");
      t = ti;                                      //     set time              
      time_read = true;                            //     set flag: time read   
    }                                              //   ENDIF                   
    in[i].close_set(nemo_io::param);               //   close nemo parameter set
    if(want) wanted |= want[i];                    //   accumulate wish list    
  }                                                // END LOOP                  
  if(time) {                                       // IF time wanted            
    *time = t;                                     //   assign to output value  
    if(!time_read && warn)                         //   IF(time not read)       
      warning("[%s]: no time found in snapshots; defaulting to 0",
	      "bodies::read_nemo_particles()");
  }                                                // ENDIF                     
  if(! time_in_range(t,times)) return false;       // IF(time not in range) DONE
  // 2. check for data support & prepare for data reading                       
  io get = wanted & NIOBits;                       // bits act'ly NEMO readable 
  if(warn && wanted != get)                        // IF wanted data not readabl
    warning("bodies::read_nemo: %s desired but not NEMO readable",
	    wanted.missing(get).make_word(miss));  //   WARNING                 
  if(get == io::o) return true;                    // IF nothing wanted: DOME   
  resetN(NN,NNS);                                  // adjust N if necessary     
  add_fields(wanted);                              // add wanted fields         
  // 3 read data                                                                
  size_t begin=0;                                  // first body to read in     
  for(int i=0; i!=Nin; ++i) {                      // LOOP nemo inputs          
    in[i].open_set(nemo_io::bodies);               //   open nemo particle set  
    if(want) get = want[i] & NIOBits;              //   data to be read         
    else     get = io::mxv;
    read[i] = read_nemo(in[i],get,begin);          //   read data               
    begin  += in[i].Number();                      //   increment begin         
    in[i].close_set(nemo_io::bodies);              //   close nemo particle set 
    in[i].reset();                                 //   reset nemo input        
    if(warn && read[i] != get) warning("bodies::read_nemo: couldn't read %s",
				       get.missing(read[i]).make_word(miss));
    if(read[i] & io::source) mark_srce_data_changed();
    if(read[i] & io::sphmax) mark_sph_data_changed ();
  }                                                // END LOOP                  
  return true;                                     // return                    
}
//------------------------------------------------------------------------------
// read a single nemo snapshop ignoring diagnose etc.                           
bool                                               // R: was time in range?     
bodies::read_nemo_snapshot(nemo_in const&input,    // I: nemo input             
			   io           &read,     // O: what has been read     
			   double       *time,     //[O: time]                  
			   io            want,     //[I: what to read]          
			   char         *times,    //[I: time range]            
			   bool          warn)     //[I: warn: missing data]    
{
  if(!input.is_present(nemo_io::snap))             // IF no snapshot present    
    falcON_ErrorF("no snapshot found","bodies::read_nemo_snapshot()");
  input.open_set(nemo_io::snap);                   // open nemo snapshot        
  bool tt = read_nemo_particles(input,read,time,   // read particles if         
				want,times,warn);  //   time in range           
  input.close_set(nemo_io::snap);                  // close nemo snapshot       
  return tt;                                       // return (time in range)?   
}
//------------------------------------------------------------------------------
// read a single snapshot combining several nemo input streams                  
bool
bodies::read_nemo_snapshots(                       // R: was time in range?     
			    const nemo_in*in,      // I: nemo inputs            
			    int           Nin,     // I: # nemo inputs          
			    io           *read,    // O: what has been read     
			    uint         *Nnum,    // O: how many have been read
			    double       *time,    //[O: time]                  
			    const io     *want,    //[I: what to read]          
			    char         *times,   //[I: time range]            
			    bool          warn)    //[I: warn: missing data]    
{
  for(int i=0; i!=Nin; ++i) {                      // LOOP nemo inputs          
    if(!in[i].is_present(nemo_io::snap))           //   IF no snapshot present  
      falcON_ErrorF("no snapshot found",
		    "bodies::read_nemo_snapshots()"); // ERROR 
    in[i].open_set(nemo_io::snap);                 // open nemo snapshot        
  }                                                // END LOOP                  
  bool tt=read_nemo_particles(in,Nin,read,Nnum,    // read particles IF         
			      time,want,times,warn); // time in range           
  for(int i=0; i!=Nin; ++i)                        // LOOP nemo inputs          
    in[i].close_set(nemo_io::snap);                //   open nemo snapshot      
  return tt;                                       // return (time in range)?   
}
//------------------------------------------------------------------------------
// write a single particle set                                                  
void                                               // write bodies to output    
bodies::write_nemo_particles(nemo_out const&out,   // I: nemo output            
			     const double  *time,  //[I: write time]            
			     io             write, //[I: what to write]         
			     uint           K,     //[I: only write first K]    
			     uint           begin) //[I: begin with this]       
  const
{
  io put = write & all_bits() & NIOBits;           // data to be written out    
  if(put == io::o) return;                         // nothing to be done        
  uint Nout = N_bodies() > begin?                  // # bodies to write out     
    (K? min(begin+K,N_bodies()) - begin : N_bodies()-begin) : 0;
  uint Sout = N_sph   () > begin?                  // # SPH bodies to write out 
    (K? min(begin+K,N_sph   ()) - begin : N_sph   ()-begin) : 0;
  // 1. open parameter set; write out N, NS & time; close parameter set;        
  out.open_set(nemo_io::param);                    // open a parameter set      
  out.write_N(Nout,Sout);                          //   write N, NS             
  if(time)                                         //   IF(time to be written)  
    out.write(nemo_io::time,*time);                //     write out time        
  out.close_set(nemo_io::param);                   // close parameter set       
  // 2. allocate memory for output; open particle set; output data;             
  if(Nout==0u && Sout==0) return;                  // nothing to be done        
  out.open_set(nemo_io::bodies);                   // open particles set        
  out.write(nemo_io::cart);                        // put out coord system      
  //    output of data is trivial, except for potential                         
  io written = io::o;
  if(put&io::p && put&io::q) {
    out.allocscalar();
    for(int i=0; i!=Nout; ++i)
      out.bodies_scl(i) = pot(i+begin) + pex(i+begin);
    out.write(nemo_io::pot);
    written |= io::pq;
  }
  nemo_read_write::set_write(out,put & ~written,begin);
  LoopIO<nemo_writer>::const_loop(this);
  // 3. close particle set & reset nemo output                                  
  out.close_set(nemo_io::bodies);
  out.reset();
}
//------------------------------------------------------------------------------
// write a single snapshot                                                      
void                                               // write snapshot            
bodies::write_nemo_snapshot(nemo_out const&out,    // I: nemo output            
			    const double  *time,   //[I: write time]            
			    io             write,  //[I: write pot/acc etc?]    
			    uint           K,      //[I: only write first K]    
			    uint           begin)  //[I: begin with this]       
  const
{
  out.open_set(nemo_io::snap);
  write_nemo_particles(out,time,write,K,begin);
  out.close_set(nemo_io::snap);
}
#endif                                             // falcON_NEMO               
//------------------------------------------------------------------------------
namespace {
  template<int IO> struct BodyReadSimple :
    protected field_io<IO> {
    static void read(std::istream&in, bodies::iterator &B) {
      in >> B. template data_io<IO>(); }
  };
  //----------------------------------------------------------------------------
  template<int IO>
  void Read(std::istream&in, bodies::iterator &B) {
    BodyReadSimple<IO>::read(in,B);
  }
  //----------------------------------------------------------------------------
 void Ignore(std::istream&in, bodies::iterator &B) { }
  //----------------------------------------------------------------------------
  typedef void (*p_reader)(std::istream&, bodies::iterator &);
}
//------------------------------------------------------------------------------
// read simple ascii formatted input                                            
void bodies::read_simple_ascii(std::istream   &in,
			       const io* const&item,
			       uint            Nb,
			       uint            Ns)
{
  // 1. create table of readers                                                 
  io  get=io::o;
  p_reader readall[32]={0}, readnon[32]={0};
  for(int K=0; K!=32 && item[K]!=io::o; ++K) {
    if(get & item[K]) 
      warning("bodies::read_simple_ascii: reading item more than once");
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
  // 2. reset N & add fields                                                    
  resetN(Nb,Ns);
  add_fields(get);
  // 3. loop bodies & read data whereby ignoring lines starting with '#'        
  LoopSPHBodies(bodies,this,Bi) {
    while( in && eat_line(in,'#') );
    if(!in) falcON_ErrorF("end of input before data have been read",
			  "bodies::read_simple_ascii()");
    for(register p_reader* R=readall; *R; ++R) (*R)(in,Bi);
    SwallowRestofLine(in);
  }
  LoopNonSPHBodies(bodies,this,Bi) {
    while( in && eat_line(in,'#') );
    if(!in) falcON_ErrorF("end of input before data have been read",
			  "bodies::read_simple_ascii()");
    for(register p_reader* R=readnon; *R; ++R) (*R)(in,Bi);
    SwallowRestofLine(in);
  }
}
//==============================================================================
