// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// yanc.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
//                             Yet Another N-body Code                         |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines                                                                     |
//                                                                             |
// class yanc                                                                  |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_yanc_h
#define falcON_included_yanc_h

#ifndef falcON_included_iostream
#  include <iostream>
#  define falcON_included_iostream
#endif

#ifndef falcON_included_nbio_h
#  include <public/nbio.h>
#endif
#ifndef falcON_included_deft_h
#  include <public/deft.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {                               // yanc is in namespace nbdy     
  class ndata;                                 // forward declare ndata         
  class sbodies;                               // forward declare bodies        
  class basic_nbody;                           // forward declare nbody         
  class extpot;                                // forward declare extpot        
  //----------------------------------------------------------------------------
  class pot_provider {
    public: virtual extpot* get_pot(const char*) = 0;
  };
  //----------------------------------------------------------------------------
  class yanc {
  private:
    yanc();                                    // not implemented               
    yanc            (const yanc&);             // not implemented               
    yanc& operator= (const yanc&);             // not implemented               
  protected:
    ndata        *MYNDATA;                     // our N-body data               
    basic_nbody  *MYNBODY;                     // our N-body code               
    const sbodies*mybodies() const;            // our N bodies                  
    const double  r_max   () const;            // max r for bodies in tree      
  public:
    //--------------------------------------------------------------------------
    yanc(                                      // construction via input file   
	 const char*,                          // I: yanc-formatted input file  
	 pot_provider* =0);                    //[I: provider for external pot] 
#ifdef falcON_NEMO
    //--------------------------------------------------------------------------
    yanc(                                      // constructor taking parameters 
	 const char*,                          // I: input file (yanc/nemo)     
	 const bool,                           // I: do nemo I/O?               
	 const double,                         // I: theta                      
	 const int,                            // I: hgrow                      
	 const int,                            // I: Ncrit                      
	 const double,                         // I: eps                        
	 const int,                            // I: kern                       
	 const int,                            // I: hmin                       
	 const int,                            // I: Nsteps                     
	 const double,                         // I: f_a                        
	 const double  = 0.,                   //[I: f_p]                       
	 const double  = 0.,                   //[I: f_c]                       
	 const double  = 0.,                   //[I: f_e]                       
#ifdef falcON_INDI
	 const double  = 0.,                   //[I: Nsoft]                     
	 const int     = 32,                   //[I: Nref]                      
	 const double  = 0.,                   //[I: emin]                      
	 const int     = 0,                    //[I: softening:global]          
#endif
	 const double  = 1.,                   //[I: Grav]                      
	 const bool    = false,                //[I: resume old (if nemo I/O)   
	 const extpot* = 0,                    //[I: external potential]        
	 const io      = io::o,                //[I: what else to read?]        
	 const int     = 2,                    //[I: # nemo output streams]     
	 const int[4]  = Default::direct);     //[I: direct sum control]        
    //--------------------------------------------------------------------------
    void  describe_nemo   (                    // describe simulation           
			   std::ostream&,      // I: output stream              
			   const char*);       // I: command line               
    //--------------------------------------------------------------------------
    void  open_nemo       (                    // open NEMO output stream       
			   const int   =0,     //[I: index of nemo stream]      
			   const char* =0,     //[I: file name ]                
			   const bool  =false);//[I" resume old sim?]           
    //--------------------------------------------------------------------------
    void  close_nemo      (                    // close NEMO output stream      
			   const int =0);      //[I: index of nemo stream]      
    //--------------------------------------------------------------------------
    void  write_nemo      (                    // NEMO-format output            
			   const io  =io::mxv, //[I: what to write out]         
			   const int =0,       //[I: index of nemo stream]      
			   const bool=true);   //[I: output diagnostics?]       
    //--------------------------------------------------------------------------
    int   nemo_devices    () const;            // R: # nemo output streams      
    //--------------------------------------------------------------------------
    bool  nemo_is_open    (                    // ready for NEMO output ?       
			   const int =0) const;//[I: index of nemo stream]      
    //--------------------------------------------------------------------------
#endif // falcON_NEMO
    //--------------------------------------------------------------------------
    ~yanc();                                   // destructor                    
    //--------------------------------------------------------------------------
    void  full_step       ();                  // perform one full time step    
    //--------------------------------------------------------------------------
    void  describe        (                    // describe simulation           
			   std::ostream&,      // I: output stream              
			   const char* =0,     //[I: name of main]              
			   const char* =0);    //[I: base name: output files]   
    //--------------------------------------------------------------------------
    static
    void  describe_yanc_format(                // describe yanc I/O format      
			       std::ostream&); // I/O: stream to write to       
    //--------------------------------------------------------------------------
    bool  write_yanc_ascii(                    // YANC-format output to file    
			   const char*,        // I: file name                  
			   const char* =0);    //[I: name of main ]             
    //--------------------------------------------------------------------------
    bool  write_yanc_binary(                   // write yanc form, bodies binary
			    const char*,       // I:   name of output file      
			    const char* =0)const;// I: name of calling program  
    //--------------------------------------------------------------------------
    void   stats          (std::ostream&) const; // write short stats           
    void   stats_head     (std::ostream&) const; // write header for short stats
    bool   okay           ()         const;    // is all fine?                  
    int    Nsteps         ()         const;    // # time steps levels           
    int    Nbodies        ()         const;    // # bodies                      
    const  char* data_file()         const;    // input data file               
    float  eps            ()         const;    // softening length              
    float  tau_min        ()         const;    // minimum time step             
    float  tau_max        ()         const;    // maximum time step             
#ifdef falcON_INDI
    int    softening      ()         const;    // type of softening             
    int    Nref           ()         const;    // # bodies in cell for rho estim
    float  Nsoft          ()         const;    // # bodies in eps sphere        
    float  emin           ()         const;    // # lower limit for eps_i       
#endif
    int    kernel         ()         const;    // type of kernel                
    double Grav           ()         const;    // Newton's G                    
    float  time           ()         const;    // current simulation time       
    float  initial_time   ()         const;    // time at input                 
    io     writing        ()         const;    // what to write out by default  
    float  theta          ()         const;    // opening parameter             
    int    Ncrit          ()         const;    // # max bodies / un-split cell  
    int    hgrow          ()         const;    // # grow fresh tree every 2^h   
                                               //   shortest steps              
    float  facc           ()         const;    // time-stepping parameter       
    float  fpot           ()         const;    // time-stepping parameter       
    float  fcom           ()         const;    // time-stepping parameter       
    float  feps           ()         const;    // time-stepping parameter       
    float  cpu_total      ()         const;    // accumulated cpu time [sec]    
    //--------------------------------------------------------------------------
  };
}
//------------------------------------------------------------------------------
#endif                                         // falcON_included_yanc_h        
