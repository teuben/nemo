// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// nbio.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
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
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
#ifdef falcON_SPH
  const int         IO_NQUANT = 20;
  const char* const IO_SQUANT ="mxvefkpParylnsTuSRDV";
#else
  const int         IO_NQUANT = 14;
  const char* const IO_SQUANT ="mxvefkpParylns";
#endif
  const size_t      IO_ZQUANT[IO_NQUANT] =
    { sizeof(real),   // m
      sizeof(vect),   // x
      sizeof(vect),   // v
      sizeof(real),   // e
      sizeof(int),    // f
      sizeof(int),    // k
      sizeof(real),   // p
      sizeof(real),   // P
      sizeof(vect),   // a
      sizeof(real),   // r
      sizeof(real),   // y
      sizeof(indx),   // l
      sizeof(uint),   // n
      sizeof(real)    // s
#ifdef falcON_SPH
     ,sizeof(real),   // T
      sizeof(real),   // u
      sizeof(real),   // S
      sizeof(real),   // R
      sizeof(real),   // D
      sizeof(vect)    // V
#endif
       };
  const char* const IO_QNAME[IO_NQUANT] = 
    { "mass",
      "pos",
      "vel",
      "eps",
      "flag",
      "key",
      "pot",
      "pex",
      "acc",
      "rho",
      "aux",
      "level",
      "num",
      "size"
#ifdef falcON_SPH
     ,"temp",
      "ein",
      "ent",
      "srho",
      "divv",
      "rotv"
#endif
    };
  class io {
    //--------------------------------------------------------------------------
    int bits;
  public:
    //--------------------------------------------------------------------------
    static int    const&N_quant()            { return IO_NQUANT; }
    static char   const&S_quant(int const&i) { return IO_SQUANT[i]; }
    static size_t const&Z_quant(int const&i) { return IO_ZQUANT[i]; }
    static size_t const P_srce ()            { return  0; }
    static size_t const N_srce ()            { return  6; }
    static size_t const P_sink ()            { return  6; }
    static size_t const N_sink ()            { return  7; }
    static size_t const P_psph ()            { return 13; }
#ifdef falcON_SPH
    static size_t const N_psph ()            { return  7; }
#else
    static size_t const N_psph ()            { return  1; }
#endif
    static int          first  (io const&bits) 
    {
      if(bits == 0) return -1;
      for(register int b=1,i=0; i!=IO_NQUANT; ++i,b<<=1)
	if(b & bits) return i;
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
	   f       = 1 << 4,                       // integer key               
	   k       = 1 << 5,                       // integer flag              
	   // data in BodySink                                                  
	   p       = 1 << 6,                       // N-body potentials         
	   P       = 1 << 7,                       // external potentials       
	   a       = 1 << 8,                       // accelarations             
	   r       = 1 << 9,                       // mass-density              
	   y       = 1 << 10,                      // auxiliary data            
	   l       = 1 << 11,                      // short level               
	   n       = 1 << 12,                      // number of neighbours      
	   // useful combinations                                               
	   xv      = x    | v,                     // phases = [x,v]            
	   mx      = m    | x,                     // masses & positions        
	   mxv     = m    | xv,                    // masses & phases           
	   mxvk    = mxv  | k,                     // masses, phases & keys     
	   mxvp    = mxv  | p,                     // masses, phases & pots     
	   mxve    = mxv  | e,                     // masses, phases & eps_i    
	   mxvf    = mxv  | f,                     // masses, phases & flags    
           mxvpaf  = mxvf | a | p,                 // elementary body data      
	   // data in BodyPsph                                                  
	   s       = 1 << 13,                      // size of SPH particle      
#ifdef falcON_SPH
	   T       = 1 << 14,                      // temperature (SPH only)    
	   u       = 1 << 15,                      // inner energy (SPH only)   
	   S       = 1 << 16,                      // entropy (SPH only)        
	   R       = 1 << 17,                      // gas-density (SPH only)    
	   D       = 1 << 18,                      // div(v)  (SPH only)        
	   V       = 1 << 19,                      // rot(v)  (SPH only)        
	   // more useful combinations                                          
	   nsR     = n | s | R,                    // elementary SPH data       
	   SPH     = mxvpaf | nsR,                 // SPH particle data         
	   sphmin  = s|R|u,                        // min sph quantities        
	   sphmax  = s|T|u|S|R|D|V                 // all sph quantities        
#else
	   sphmin  = s,                            // min sph quantities        
	   sphmax  = s                             // max sph quantities        
#endif
    };
    //==========================================================================
    // construction                                                             
    //--------------------------------------------------------------------------
    io()             : bits(io::o) {}
    io(const int  i) : bits(i) {}
    io(const io&  i) : bits(i.bits) {}
    //--------------------------------------------------------------------------
    io(char const&c) : bits(io::o) {
      for(register int i=0; i!=N_quant(); ++i)
	if(c == S_quant(i)) bits |= 1<<i;
    }
    //--------------------------------------------------------------------------
    io(const char*c) : bits(io::o) {
      for(register int i=0; i!=N_quant(); ++i)
	if(std::strchr(c,S_quant(i))) bits |= 1<<i;
    }
    //==========================================================================
    void make_word(char* w) const {
      char *letter = w;
      for(register int i=0; i!=N_quant(); ++i)
	if(bits & 1<<i) *(letter++) = S_quant(i);
      *letter = 0;
    }
    //--------------------------------------------------------------------------
    inline size_t  bytes       ()            const {
      register size_t n = 0;
      for(register int i=0; i!=N_quant(); ++i)
	if(bits & 1<<i) n += Z_quant(i);
      return n;
    }
    //--------------------------------------------------------------------------
    const char* word() const;
    friend const char* word(const io &i) { return i.word(); }
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
