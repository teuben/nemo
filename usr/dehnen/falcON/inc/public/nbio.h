// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// nbio.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_nbio_h
#define falcON_included_nbio_h

#ifndef falcON_included_cstring
#  include <cstring>                               // C style charaters         
#  define falcON_included_cstring
#endif
#ifndef falcON_included_iostream
#  include <iostream>                              // C++ I/O                   
#  define falcON_included_iostream
#endif
#ifndef falcON_included_auxx_h
#  include <public/auxx.h>
#endif
#ifndef falcON_included_flag_h
#  include <public/flag.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  const int         IO_NOSPH   = 13;
#ifdef falcON_SPH
  const int         IO_NQUANT  = 26;
  const char* const IO_SQUANT  ="mxvefkpqarylnHNUYIESRDVQTJ";
#else
  const int         IO_NQUANT  = 14;
  const char* const IO_SQUANT  ="mxvefkpqarylnH";
#endif
  //----------------------------------------------------------------------------
  const char* const IO_QNAME[IO_NQUANT] = 
    { "mass", "position", "velocity", "epsilon", "flag", "key",
      "potential", "potential_external", "acceleration",
      "density",  "auxiliary", "level", "number_of_partners", "size"
#ifdef falcON_SPH
     ,"number_of_SPH_partners", "U_internal", "U_predicted", "(dU/dt)_internal",
      "(dU/dt)_external", "entropy", "gas-density", "dgas-density/dt",
      "velocity_predicted", "sigma-squared","temperature","dsize/dt"
#endif
    };
  //----------------------------------------------------------------------------
  const size_t IO_ZQUANT[IO_NQUANT] = {
    //            source properties: 6
    sizeof(real), //  0: mass
    sizeof(vect), //  1: position
    sizeof(vect), //  2: velocity
    sizeof(real), //  3: softening length
    sizeof(flag), //  4: flag
    sizeof(int),  //  5: key
    //            sink properties: 7
    sizeof(real), //  6: potential
    sizeof(real), //  7: external potential
    sizeof(vect), //  8: acceleration
    sizeof(real), //  9: density (mass)
    sizeof(real), // 10: auxiliary
    sizeof(indx), // 11: level
    sizeof(uint), // 12: number of partners
    //            SPH properties: 1/12
    sizeof(real)  // 13: size
#ifdef falcON_SPH
   ,sizeof(uint), // 14: number of SPH partners
    sizeof(real), // 15: U_internal
    sizeof(real), // 16: U_predicted
    sizeof(real), // 17: (dU/dt)_internal
    sizeof(real), // 18: (dU/dt)_external
    sizeof(real), // 19: entropy
    sizeof(real), // 20: gas-density
    sizeof(real), // 21: dgas-density/dt
    sizeof(vect), // 22: velocity_predicted
    sizeof(real), // 23: sigma-squared
    sizeof(real), // 24: temperature
    sizeof(real)  // 25: dh/dt
#endif
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::io                                                           //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class io {
    //--------------------------------------------------------------------------
    int bits;
  public:
    //--------------------------------------------------------------------------
    static int    const N_quant = IO_NQUANT;
    static char   const&S_quant(int const&i) { return IO_SQUANT[i]; }
    static size_t const&Z_quant(int const&i) { return IO_ZQUANT[i]; }
    static size_t       P_srce ()            { return  0; }
    static size_t       N_srce ()            { return  6; }
    static size_t       P_sink ()            { return  6; }
    static size_t       N_sink ()            { return  7; }
    static size_t       P_psph ()            { return 13; }
#ifdef falcON_SPH
    static size_t       N_psph ()            { return 13; }
#else
    static size_t       N_psph ()            { return  1; }
#endif
    static int          first  (io const&bits) 
    {
      if(bits == 0) return -1;
      for(register int b=1,i=0; i!=N_quant; ++i,b<<=1)
	if(b & bits) return i;
      return -1;
    }
    static const char* first_name  (io const&bits)
    { return bits? IO_QNAME[first(bits)] : "nothing"; }
    //--------------------------------------------------------------------------
    enum { o       = 0,                            // nothing written/to be read
	   // data in BodySrce                                                  
	   m       = 1 << 0,                       // masses                    
	   x       = 1 << 1,                       // positions                 
	   v       = 1 << 2,                       // velocities                
	   e       = 1 << 3,                       // indiv eps_i               
	   f       = 1 << 4,                       // integer flag              
	   k       = 1 << 5,                       // integer key               
	   // data in BodySink                                                  
	   p       = 1 << 6,                       // N-body potentials         
	   q       = 1 << 7,                       // external potentials       
	   a       = 1 << 8,                       // accelarations             
	   r       = 1 << 9,                       // mass-density              
	   y       = 1 << 10,                      // auxiliary data            
	   l       = 1 << 11,                      // short level               
	   n       = 1 << 12,                      // number of neighbours      
	   // useful combinations                                               
	   ap      = a    | p,                     // just gravity              
	   xv      = x    | v,                     // phases = [x,v]            
	   mx      = m    | x,                     // masses & positions        
	   mxv     = m    | xv,                    // masses & phases           
	   mxvk    = mxv  | k,                     // masses, phases & keys     
	   mxvp    = mxv  | p,                     // masses, phases & pots     
	   mxve    = mxv  | e,                     // masses, phases & eps_i    
	   mxvke   = mxvk | e,                     // masses, phases, keys & eps
	   mxvf    = mxv  | f,                     // masses, phases & flags    
           mxvpaf  = mxvf | a | p,                 // elementary body data      
	   gravity = mxvpaf,
	   // data in BodyPsph are all uppercase                                
	   H       = 1 << 13,                      // size h of SPH particle    
#ifdef falcON_SPH                                  // these are only for SPH:   
	   N       = 1 << 14,                      //   # SPH partners          
	   U       = 1 << 15,                      //   U = inner energy        
	   Y       = 1 << 16,                      //   predicted U             
	   I       = 1 << 17,                      //   internal dU/dt          
	   E       = 1 << 18,                      //   external dU/dt          
	   S       = 1 << 19,                      //   entropy                 
	   R       = 1 << 20,                      //   rho := gas-density      
	   D       = 1 << 21,                      //   drho/dt                 
	   V       = 1 << 22,                      //   predicted V             
	   Q       = 1 << 23,                      //   sigma^2                 
	   T       = 1 << 24,                      //   temperature             
	   J       = 1 << 25,                      //   dh/dt                   
	   // more useful combinations                                          
	   NHRVD   = N|H|R|V|D,                    // elementary SPH data       
	   NHRVDJ  = NHRVD|J,                      // elementary SPH data       
	   Ud      = U|Y|I,                        // internal energy data      
	   NHRVDU  = NHRVD|Ud,                     // SPH data with internal U  
	   SPH     = mxvpaf | NHRVDJ,              // SPH particle data         
	   SPHU    = SPH | Ud,                     // SPH particle data         
	   sphmin  = N|H|R|V,                      // min sph quantities        
	   sphdef  = N|H|R|V|D|J,                  // def sph quantities        
	   sphmax  = N|H|T|U|Y|I|E|S|R|D|V|Q|J,    // all sph quantities        
	   sphnemo = N|H|U|I|E|S|R,                // sph quants with nemo I/O  
#else
	   sphmin  = H,                            // min sph quantities        
	   sphdef  = H,                            // def sph quantities        
	   sphmax  = H,                            // max sph quantities        
	   sphnemo = H,                            // sph quants with nemo I/O  
#endif
	   NEMO    = mxv|e|f|k|p|q|a|r|y|l|n|sphnemo
    };
    //==========================================================================
    // construction                                                             
    //--------------------------------------------------------------------------
    io()             : bits(io::o) {}
    io(const int  i) : bits(i) {}
    io(const io&  i) : bits(i.bits) {}
    //--------------------------------------------------------------------------
    io(char const&c) : bits(io::o) {
      for(register int i=0; i!=N_quant; ++i)
	if(c == S_quant(i)) bits |= 1<<i;
    }
    //--------------------------------------------------------------------------
    io(const char*c) : bits(io::o) {
      for(register int i=0; i!=N_quant; ++i)
	if(std::strchr(c,S_quant(i))) bits |= 1<<i;
    }
    //==========================================================================
    void make_word(char* w) const {
      char *letter = w;
      for(register int i=0; i!=N_quant; ++i)
	if(bits & 1<<i) *(letter++) = S_quant(i);
      *letter = 0;
    }
    //--------------------------------------------------------------------------
    inline size_t  bytes       ()            const {
      register size_t n = 0;
      for(register int i=0; i!=N_quant; ++i)
	if(bits & 1<<i) n += Z_quant(i);
      return n;
    }
    //==========================================================================
    io&  operator=      (const io &i)       { bits =i.bits; return *this; }
    io&  operator|=     (const io &i)       { bits|=i.bits; return *this; }
    io&  operator|=     (const int&i)       { bits|=i; return *this; }
    io&  operator&=     (const io &i)       { bits&=i.bits; return *this; }
    io&  operator&=     (const int&i)       { bits&=i; return *this; }
    //==========================================================================
    operator int        ()            const { return bits; }
    //==========================================================================
    friend std::ostream& operator<< (std::ostream&s, const io&i)
    { return s<<i.bits; }
    friend std::istream& operator>> (std::istream&s,       io&i)
    { return s>>i.bits; }
    //==========================================================================
    io   operator|      (const io &i) const { return io(bits | i.bits); }
    io   operator|      (const int&i) const { return io(bits | i); }
    io   operator&      (const io &i) const { return io(bits & i.bits); }
    io   operator&      (const int&i) const { return io(bits & i); }
    bool operator==     (const io &i) const { return bits == i.bits; }
    bool operator==     (const int&i) const { return bits == i; }
    bool operator!=     (const io &i) const { return bits != i.bits; }
    bool operator!=     (const int&i) const { return bits != i; }
    bool contains       (const io &i) const { return (bits & i.bits)==i.bits; }
    io   missing        (const io &i) const { return io((i^bits) & i); }
    io   operator~      ()            const { return io(~bits); }
  };
}
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_nbio_h
