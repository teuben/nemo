// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// stic.h                                                                      |
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
// defines                                                                     |
//                                                                             |
// class sticky_soul                                                           |
// class sticky_cell                                                           |
// class sticky_tree                                                           |
// class basic_finder                                                          |
// class sticky_finder                                                         |
// class neighbour_finder                                                      |
// class neighbour_counter                                                     |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_stic_h
#define falcON_included_stic_h

#ifndef falcON_included_grat_h
#  include <public/grat.h>
#endif
#ifndef falcON_included_ionl_h
#  include <public/ionl.h>
#endif
////////////////////////////////////////////////////////////////////////////////

namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::sticky_soul                                                  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class sticky_soul : public basic_soul {
    sticky_soul           (const sticky_soul&);    // not implemented           
    sticky_soul& operator=(const sticky_soul&);    // not implemented           
    //--------------------------------------------------------------------------
    // data of class sticky_soul                                                
    //--------------------------------------------------------------------------
  private:
    real SIZE;                                     // physical size             
    real SIZEQ;                                    // (physical size)^ 2        
    vect VEL;                                      // velocity                  
    uint NUM;                                      // # neighbours              
  public:
    //--------------------------------------------------------------------------
    // non-const methods                                                        
    //--------------------------------------------------------------------------
    void inc() { NUM++; }
    //--------------------------------------------------------------------------
    // const data access via friends                                            
    //--------------------------------------------------------------------------
    friend uint const&mybody(const sticky_soul*const&O) { return O->mybody(); } 
    friend uint const&num   (const sticky_soul*const&O) { return O->NUM; } 
    friend vect const&pos   (const sticky_soul*const&O) { return O->pos(); } 
    friend vect const&vel   (const sticky_soul*const&O) { return O->VEL; } 
    friend real const&size  (const sticky_soul*const&O) { return O->SIZE; } 
    friend real const&sizeq (const sticky_soul*const&O) { return O->SIZEQ; } 
    //--------------------------------------------------------------------------
    // simple manipulations for use with bodies                                 
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void set_sticky(const bodies_type*const&B) {
      SIZE=B->size(mybody());
      VEL =B->vel (mybody());
      NUM =0;
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void set_sph   (const bodies_type*const&B) {
      SIZE =B->size(mybody());
      SIZEQ=SIZE*SIZE;
      NUM  =0;
    }
    //--------------------------------------------------------------------------
    // simple manipulations for use with barrays                                
    //--------------------------------------------------------------------------
    void set_sticky(const barrays*const&B) {
      SIZE  =B->size (mybody());
      VEL[0]=B->vel_x(mybody());
      VEL[1]=B->vel_y(mybody());
#if falcON_NDIM > 2
      VEL[2]=B->vel_z(mybody());
#endif
      NUM=0;
    }
    //--------------------------------------------------------------------------
    // dump data                                                                
    //--------------------------------------------------------------------------
    static void dump_head(std::ostream& o) {
      o<<"     #";
      basic_soul::dump_head(o);
      o<<"           size [velocity]";
    }
    //--------------------------------------------------------------------------
    static void dump(std::ostream            &o, 
		     const sticky_soul* const&S,
		     int                const&index) {
      o<<' '<<setw(5)<<index;
      basic_soul::dump(o,S);
      o<<' '<<setw(5)<<setprecision(4)<<nbdy::size(S);
      if(S->is_set(flag::STICKY))
	for(register indx d=0; d!=Ndim; ++d)
	  o<<' '<<setw(7)<<setprecision(4)<<nbdy::vel(S)[d];
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::sticky_cell                                                  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class sticky_cell : public basic_cell {
    sticky_cell           (const sticky_cell&);    // not implemented           
    sticky_cell& operator=(const sticky_cell&);    // not implemented           
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
  public:
    typedef sticky_soul soul_type;                 // type of associated souls  
    //--------------------------------------------------------------------------
    // data of class cell                                                       
    //--------------------------------------------------------------------------
  private:
    vect    POS;                                   // position center           
    vect    VEL;                                   // velocity center           
    real    SIZE;                                  // size                      
    union {
      real  VRAD;                                  // radius of velocity sphere 
      real  RMAX;                                  // radius of position sphere 
    };
  public:
    //--------------------------------------------------------------------------
    // non-const data access via members                                        
    //--------------------------------------------------------------------------
    vect&vel () { return VEL; }
    vect&pos () { return POS; }
    real&rmax() { return RMAX; }
    real&size() { return SIZE; }
    real&vrad() { return VRAD; }
    //--------------------------------------------------------------------------
    // const data access via friends                                            
    //--------------------------------------------------------------------------
    friend real const&radius(const sticky_cell*const&O) { return O->radius(); } 
    friend vect const&center(const sticky_cell*const&O) { return O->center(); } 
    friend int  const&number(const sticky_cell*const&O) { return O->number(); } 
    friend vect const&vel   (const sticky_cell*const&O) { return O->VEL; } 
    friend vect const&pos   (const sticky_cell*const&O) { return O->POS; } 
    friend real const&rmax  (const sticky_cell*const&O) { return O->RMAX; } 
    friend real const&size  (const sticky_cell*const&O) { return O->SIZE; } 
    friend real const&vrad  (const sticky_cell*const&O) { return O->VRAD; } 
    //--------------------------------------------------------------------------
    // boolean information via friends                                          
    //--------------------------------------------------------------------------
    friend bool has_cell_kids(const sticky_cell*const&C) {
      return C->has_cell_kids(); }
    friend bool has_soul_kids(const sticky_cell*const&C) {
      return C->has_soul_kids(); }
    friend bool is_twig          (const sticky_cell*const&C) {
      return C->is_twig(); }
    //--------------------------------------------------------------------------
    // flag manipulation                                                        
    //--------------------------------------------------------------------------
    void add_active_flag(const sticky_soul* const&S)
    {
      add_active_flag_from_soul(S);
    }
    //--------------------------------------------------------------------------
    void add_active_flag(const sticky_cell* const&C)
    {
      add_active_flag_from_cell(C);
    }
    //--------------------------------------------------------------------------
    // dump  data                                                               
    //--------------------------------------------------------------------------
    static void dump_head(std::ostream&o) {
      o<<"     #";
      basic_cell::dump_head(o);
      o<<"           size     rmax/ velocity           [vrad]";
    }
    //--------------------------------------------------------------------------
    static void dump(std::ostream            &o, 
		     const sticky_cell* const&C,
		     int                const&index) {
      o<<' '<<setw(5)<<index;
      basic_cell::dump(o,C);
      o<<' '<<setw(6)<<setprecision(4)<<nbdy::size(C);
      if(C->is_set(flag::STICKY)) {
	for(register indx d=0; d!=Ndim; ++d)
	  o<<' '<<setw(9)<<setprecision(4)<<nbdy::vel(C)[d];
	o<<' '<<setw(5)<<setprecision(4)<<nbdy::vrad(C);
      }
      else if (C->is_set(flag::SPH))
	o<<' '<<setw(5)<<setprecision(4)<<nbdy::rmax(C);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::sticky_tree                                                    
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class sticky_tree : public basic_tree<sticky_tree,sticky_cell> {
    sticky_tree           (const sticky_tree&);    // not implemented           
    sticky_tree& operator=(const sticky_tree&);    // not implemented           
    //--------------------------------------------------------------------------
    // static data                                                              
    //--------------------------------------------------------------------------
    static const int SPEC=flag::SPH|flag::STICKY;  // default flag for sub-tree 
  public:
    typedef uint elem_pair[2];                     // element: interaction list 
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
    sticky_tree(const grav_tree*const&parent,      // I: parent tree            
		int             const&spec=SPEC) : //[I: flag for sub-tree]     
      base_tree(parent, spec) {}
    //--------------------------------------------------------------------------
    void count_neighbours() const;
    //--------------------------------------------------------------------------
    void make_iaction_list(elem_pair *,
			   uint const&,
			   uint      &,
			   real const&) const;
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // const access to stick_cell data from cell_iterator                       //
  //                                                                          //
  // these are not really necessary, since they are already defined for taking//
  // pointer to cell and cell_iterator has a type-conversion operator defined,//
  // but some older compilers get confused if we omit them (eg. gcc 2.95.3    //
  // needs them, gcc 3.2 not)                                                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
#define CSCI sticky_tree::cell_iterator const&I
  inline vect const&vel (CSCI) { return vel (I.c_pter()); }
  inline vect const&pos (CSCI) { return pos (I.c_pter()); }
  inline real const&rmax(CSCI) { return rmax(I.c_pter()); }
  inline real const&size(CSCI) { return size(I.c_pter()); }
  inline real const&vrad(CSCI) { return vrad(I.c_pter()); }
#undef CSCI
}
//------------------------------------------------------------------------------
#endif // falcON_included_stic_h
