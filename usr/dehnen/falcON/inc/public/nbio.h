// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// nbio.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2001                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef included_nbio_h
#define included_nbio_h

#ifndef included_cstring
#  include <cstring>                               // C style charaters         
#  define included_cstring
#endif
#ifndef included_iostream
#  include <iostream>                              // C++ I/O                   
#  define included_iostream
#endif
#ifndef included_auxx_h
#  include <public/auxx.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  class io {
    int bits;
  public:
    enum { o      = 0,                             // nothing written/to be read
	   m      = 1 << 0,                        // masses                    
	   x      = 1 << 1,                        // positions                 
	   v      = 1 << 2,                        // velocities                
	   e      = 1 << 3,                        // indiv eps_i               
	   p      = 1 << 4,                        // N-body potentials         
	   a      = 1 << 5,                        // accelarations             
	   r      = 1 << 6,                        // density                   
	   y      = 1 << 7,                        // auxiliary data            
	   l      = 1 << 8,                        // short level               
	   k      = 1 << 9,                        // integer key               
	   f      = 1 << 10,                       // integer flag              
	   n      = 1 << 11,                       // number of neighbours      
	   s      = 1 << 12,                       // size of SPH particle      
	   xv     = x    | v,                      // phases = [x,v]            
	   mx     = m    | x,                      // masses & positions        
	   mxv    = m    | xv,                     // masses & phases           
	   mxvk   = mxv  | k,                      // masses, phases & keys     
	   mxvp   = mxv  | p,                      // masses, phases & pots     
	   mxve   = mxv  | e,                      // masses, phases & eps_i    
	   mxvf   = mxv  | f,                      // masses, phases & flags    
           mxvpaf = mxvf | a | p } ;               // elementary body data      
    //==========================================================================
    // construction                                                             
    //--------------------------------------------------------------------------
    io()             : bits(io::o) {}
    io(const int  i) : bits(i) {}
    io(const io&  i) : bits(i.bits) {}
    //--------------------------------------------------------------------------
    io(char const&c) : bits(io::o) {
      const char* set ="mxveparylkfns";
      for(register int i=0; i!=12; ++i)
	if(c == set[i]) bits |= 1<<i;
    }
    //--------------------------------------------------------------------------
    io(const char*c) : bits(io::o) {
      const char* set ="mxveparylkfns";
      for(register int i=0; i!=12; ++i)
	if(std::strchr(c,set[i])) bits |= 1<<i;
    }
    //==========================================================================
    void make_word(char* w) const {
      const char* set ="mxveparylkfns";
      char *letter = w;
      for(register int i=0; i!=12; ++i)
	if(bits & 1<<i) *(letter++) = set[i];
      *letter = 0;
    }
    //--------------------------------------------------------------------------
    const char* word() const;
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
    io   operator~      ()            const { return io(~bits); }
    size_t  bytes       ()            const {
      register size_t n = 0;
      if(bits & m) n+= sizeof(real);
      if(bits & x) n+= sizeof(vect);
      if(bits & v) n+= sizeof(vect);
      if(bits & e) n+= sizeof(real);
      if(bits & p) n+= sizeof(real);
      if(bits & a) n+= sizeof(vect);
      if(bits & r) n+= sizeof(real);
      if(bits & y) n+= sizeof(real);
      if(bits & l) n+= sizeof(indx);
      if(bits & k) n+= sizeof(int );
      if(bits & f) n+= sizeof(int );
      if(bits & n) n+= sizeof(uint);
      if(bits & s) n+= sizeof(real);
      return n;
    }
  };
}
////////////////////////////////////////////////////////////////////////////////
#endif // included_nbio_h
