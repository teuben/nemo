// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// flag.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines                                                                     |
//                                                                             |
// class flag                     fully inline                                 |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_flag_h
#define falcON_included_flag_h 1

////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  // class nbdy::flag                                                           
  //////////////////////////////////////////////////////////////////////////////
  class flag {
    //--------------------------------------------------------------------------
    // type of class flag                                                       
    //--------------------------------------------------------------------------
  public:
    enum {
      ACTIVE        =  1<<0,                    // tree node is active          
      NOT_IN_TREE   =  1<<1,                    // body ignored by tree         
      SPH           =  1<<2,                    // tree node has sph bodies     
      STICKY        =  1<<3,                    // tree node has sticky bodies  
      AL_ACTIVE     =  1<<4,                    // all children are actives     
      AL_SPH        =  1<<5,                    // all children are sph         
      AL_STICKY     =  1<<6,                    // all children are sticky      
      // compounds of flags                                                     
      LEAF_FLAGS    =  ACTIVE|SPH|STICKY,       // sum of all leaf flags        
      BODY_FLAGS    =  LEAF_FLAGS|NOT_IN_TREE,  // sum of all body flags        
      ACTIVE_FLAGS  =  ACTIVE|AL_ACTIVE,        // sum of all active flags      
      SPH_FLAGS     =  SPH   |AL_SPH,           // sum of all sph  flags        
      STICKY_FLAGS  =  STICKY|AL_STICKY         // sum of all sticky flags      
    };
    //--------------------------------------------------------------------------
    // data of class flag                                                       
    //--------------------------------------------------------------------------
  private:
    int FLAG;                                   // our flag                     
    //--------------------------------------------------------------------------
    // constructors                                                             
    //--------------------------------------------------------------------------
  public:
    flag           ()             : FLAG(0) {}
    flag           (const int &F) : FLAG(F) {}
    flag           (const flag&F) : FLAG(F.FLAG) {}
    flag& operator=(const flag&F) { FLAG = F.FLAG; return *this; }
    //--------------------------------------------------------------------------
    // non const methods                                                        
    //--------------------------------------------------------------------------
    void reset       ()                         { FLAG = 0; }
    void set_to      (const int&F)              { FLAG = F; }
    void set_to_part (const int&F, const int&P) { FLAG = F&P; }
    void add         (const int&F)              { FLAG|= F; }
    void add_part    (const int&F, const int&P) { FLAG|= F&P; }
    void set_part    (const int&F, const int&P) { FLAG = FLAG&~P | F&P; }
    void un_set      (const int&F)              { FLAG&= ~F; }
    //--------------------------------------------------------------------------
    void set_to      (const flag*F)             { FLAG = F->FLAG; }
    void set_to_part (const flag*F,const int&P) { FLAG = F->FLAG&P;}
    void add         (const flag*F)             { FLAG|= F->FLAG; }
    void add_part    (const flag*F,const int&P) { FLAG|= F->FLAG&P;}
    void set_part    (const flag*F,const int&P) { FLAG = FLAG&~P | F->FLAG&P; }
    void un_set      (const flag*F)             { FLAG&= ~F->FLAG; }
    //--------------------------------------------------------------------------
    void set_to      (const flag&F)             { FLAG = F.FLAG; }
    void set_to_part (const flag&F,const int&P) { FLAG = F.FLAG&P;}
    void add         (const flag&F)             { FLAG|= F.FLAG; }
    void add_part    (const flag&F,const int&P) { FLAG|= F.FLAG&P;}
    void set_part    (const flag&F,const int&P) { FLAG = FLAG&~P | F.FLAG&P; }
    void un_set      (const flag&F)             { FLAG&= ~F.FLAG; }
    //--------------------------------------------------------------------------
    // boolean information via members                                          
    //--------------------------------------------------------------------------
    bool is_set(const int &T) const { return FLAG & T; }
    //--------------------------------------------------------------------------
    // boolean information via friends                                          
    //--------------------------------------------------------------------------
    friend bool is_active (const flag*F) { return F->is_set(ACTIVE); }
    friend bool al_active (const flag*F) { return F->is_set(AL_ACTIVE); }
    friend bool is_in_tree(const flag*F) { return!F->is_set(NOT_IN_TREE); }
    friend bool is_sph    (const flag*F) { return F->is_set(SPH); }
    friend bool al_sph    (const flag*F) { return F->is_set(AL_SPH); }
    friend bool is_sticky (const flag*F) { return F->is_set(STICKY); }
    friend bool al_sticky (const flag*F) { return F->is_set(AL_STICKY); }
    friend bool is_set    (const flag*F, int const&T) { return F->is_set(T); }
    //--------------------------------------------------------------------------
    friend bool is_active (flag const&F) { return F.is_set(ACTIVE); }
    friend bool al_active (flag const&F) { return F.is_set(AL_ACTIVE); }
    friend bool is_in_tree(flag const&F) { return!F.is_set(NOT_IN_TREE); }
    friend bool is_sph    (flag const&F) { return F.is_set(SPH); }
    friend bool al_sph    (flag const&F) { return F.is_set(AL_SPH); }
    friend bool is_sticky (flag const&F) { return F.is_set(STICKY); }
    friend bool al_sticky (flag const&F) { return F.is_set(AL_STICKY); }
    friend bool is_set    (flag const&F, int const&T) { return F.is_set(T); }
    //--------------------------------------------------------------------------
    friend bool is_active (int const&F) { return F & ACTIVE; }
    friend bool al_active (int const&F) { return F & AL_ACTIVE; }
    friend bool is_in_tree(int const&F) { return!(F & NOT_IN_TREE); }
    friend bool is_sph    (int const&F) { return F & SPH; }
    friend bool al_sph    (int const&F) { return F & AL_SPH; }
    friend bool is_sticky (int const&F) { return F & STICKY; }
    friend bool al_sticky (int const&F) { return F & AL_STICKY; }
    friend bool is_set    (int const&F, int const&T) { return F & T; }
    //--------------------------------------------------------------------------
    // conversion to int                                                        
    //--------------------------------------------------------------------------
    operator const int& () const { return FLAG; }
    operator       int& ()       { return FLAG; }
    //--------------------------------------------------------------------------
    // I/O                                                                      
    //--------------------------------------------------------------------------
    friend std::ostream& operator<< (std::ostream& o, const flag& F)
    { return o<<F.FLAG; }
    friend std::istream& operator>> (std::istream& i, flag& F)
    { return i>>F.FLAG; }
  };
}
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_flag_h
