// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// stic.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2002                                          |
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
#ifndef included_stic_h
#define included_stic_h

#ifndef included_grat_h
#  include <public/grat.h>
#endif
#ifndef included_iact_h
#  include <public/iact.h>
#endif
#ifndef included_ionl_h
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
    //--------------------------------------------------------------------------
    // constructor                                                              
    //--------------------------------------------------------------------------
  public:
    sticky_soul() {}
    //--------------------------------------------------------------------------
    // non-const methods                                                        
    //--------------------------------------------------------------------------
    void inc() { NUM++; }
    //--------------------------------------------------------------------------
    // const data access via friends                                            
    //--------------------------------------------------------------------------
    friend vect const&cofm  (const sticky_soul*const&O) { return O->cofm(); } 
    friend uint const&mybody(const sticky_soul*const&O) { return O->mybody(); } 
    friend uint const&num   (const sticky_soul*const&O) { return O->NUM; } 
    friend vect const&pos   (const sticky_soul*const&O) { return O->cofm(); } 
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
    // simple manipulations for use with arrays                                 
    //--------------------------------------------------------------------------
    void set_sticky(const areal*const&s, const areal*v[NDIM]) {
      SIZE=s[mybody()];
      for(register int d=0; d<NDIM; ++d) VEL[d]=v[d][mybody()];
      NUM=0;
    }
    //--------------------------------------------------------------------------
    void set_sph   (const areal*const&s) {
      SIZE =s[mybody()];
      SIZEQ=SIZE*SIZE;
      NUM  =0;
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
	for(register indx d=0; d<NDIM; d++)
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
    //--------------------------------------------------------------------------
    // constructor                                                              
    //--------------------------------------------------------------------------
  public:
    sticky_cell() {}
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
    void add_sink_flag(const sticky_soul* const&S) {add_sink_flag_from_soul(S);}
    void add_sink_flag(const sticky_cell* const&C) {add_sink_flag_from_cell(C);}
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
	for(register indx d=0; d<NDIM; ++d)
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
    // data of class sticky_tree                                                
    //--------------------------------------------------------------------------
  private:
    const areal  *S,**V;                           // arrays with sizes and vels
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    sticky_tree(grav_tree *parent) :
      base_tree(parent), S(0), V(0)
    { if(!use_sbodies()) NbdyError("bodies/array mismatch") }
    //--------------------------------------------------------------------------
    sticky_tree(grav_tree *parent,
		const areal *s,
		const areal *v[NDIM]) :
      base_tree(parent), S(s), V(v)
    { if(!use_arrays()) NbdyError("bodies/array mismatch") }
    //--------------------------------------------------------------------------
    void prepare_sticky();
    void prepare_sph();
    void set_num(grav_soul*) const;
    void set_num(int *) const;
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
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::basic_finder                                                   
  //                                                                            
  // for finding soul pairs in a sticky_tree;                                   
  // derived from basic_iactor of iact.h, which satisfies the requirements for  
  // an INTERACTOR template parameter to class MutualInteractor<>;              
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class basic_finder : public basic_iactor<sticky_tree> {
    //--------------------------------------------------------------------------
    // types of class basic_finder                                              
    //--------------------------------------------------------------------------
  public:
    typedef uint elem_pair[2];                     // pair of indices           
    //--------------------------------------------------------------------------
    // data of class basic_finder                                               
    //--------------------------------------------------------------------------
  private:
    const uint    MAX;                             // maximal size of list      
    elem_pair    *BL;                              // list of interaction pairs 
    mutable uint  N;                               // actual size of list       
    //--------------------------------------------------------------------------
    // other methods                                                            
    //--------------------------------------------------------------------------
  protected:
    void add_pair(soul_iter const&, soul_iter const&) const;
  public:
    uint const &actual_size_of_list() const { return N; }
    //--------------------------------------------------------------------------
    // constructors                                                             
    //--------------------------------------------------------------------------
    basic_finder(const uint n, elem_pair*l) : MAX(n),BL(l),N(0){}
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::sticky_finder                                                  
  //                                                                            
  // for finding sticky_soul pairs which satisfy:                               
  // - at least one is flagged as sink                                          
  // - there is a time t in [0,tau] such that:                                  
  //   | (x_i+v_i*t) - (x_j+v+j*t) | < size_i+size_j                            
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class sticky_finder : public basic_finder {
  private:
    const real    TAU;                             // time period               
  protected:
    void many       (bool const&, soul_iter const&,
		     soul_iter const&, soul_iter const&) const;
    void single     (soul_iter const&, soul_iter const&) const;
    bool discard    (cell_iter const&, soul_iter const&) const;
    bool discard    (cell_iter const&, cell_iter const&) const;
  public:
    sticky_finder(const real t, const uint n, elem_pair*l) :
      basic_finder(n,l), TAU(t) {}
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return is_twig(B) || size(A)+TAU*vrad(A) > size(B)+TAU*vrad(B);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::neighbour_finder                                               
  //                                                                            
  // for finding sticky_soul pairs which satisfy:                               
  // - at least one is flagged as sink                                          
  // - | x_i - x_j | < max(size_i,size_j)                                       
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class neighbour_finder : public basic_finder {
  protected:
    void many       (bool const&, soul_iter const&,
		     soul_iter const&, soul_iter const&) const;
    void single     (soul_iter const&, soul_iter const&) const;
    bool discard    (cell_iter const&, soul_iter const&) const;
    bool discard    (cell_iter const&, cell_iter const&) const;
  public:
    neighbour_finder(const uint n, elem_pair*l) : basic_finder(n,l) {}
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return is_twig(B) || size(A) > size(B);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::neighbour_counter                                            //
  //                                                                          //
  // for each sink soul: count # souls with |R| < size(A);                    //
  // for individual sizes(soft_type=individual): requires:                    //
  //     - each sink soul to hold its size and size^2.                        //
  //     - each sink cell to hold its size                                    //
  // for global size, it is given as argument to the constructor              //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename TREE, soft_type SOFT>  class neighbour_counter;
  //----------------------------------------------------------------------------
  template<typename TREE>
  class neighbour_counter<TREE,individual> : public basic_iactor<TREE> {
    typedef typename TREE::soul_iterator soul_iter;
    typedef typename TREE::cell_iterator cell_iter;
  protected:
    void many       (bool const&, soul_iter const&,
		     soul_iter const&, soul_iter const&) const;
    void single     (soul_iter const&, soul_iter const&) const;
    bool discard    (cell_iter const&, soul_iter const&) const;
    bool discard    (cell_iter const&, cell_iter const&) const;
  public:
    neighbour_counter() {}
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return nbdy::is_twig(B) || size(A) > size(B);
    }
  };
  //----------------------------------------------------------------------------
  template<typename TREE>
  class neighbour_counter<TREE,global>
    : public basic_iactor<TREE> {
    typedef typename TREE::soul_iterator soul_iter;
    typedef typename TREE::cell_iterator cell_iter;
  private:
    const real EPS,EPQ;
  protected:
    void many       (bool const&, soul_iter const&,
		     soul_iter const&, soul_iter const&) const;
    void single     (soul_iter const&, soul_iter const&) const;
    bool discard    (cell_iter const&, soul_iter const&) const;
    bool discard    (cell_iter const&, cell_iter const&) const;
  public:
    neighbour_counter(const real e) : EPS(e), EPQ(e*e) {}
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return is_twig(B) || size(A) > size(B); }
  };
  //////////////////////////////////////////////////////////////////////////////
}
//------------------------------------------------------------------------------
#endif // included_stic_h
