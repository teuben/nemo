// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// stic.h                                                                      |
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
// class stsp_leaf                                                             |
// class stsp_cell                                                             |
// class stsp_estimator                                                        |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_stic_h
#define falcON_included_stic_h

#ifndef falcON_included_grav_h
#  include <public/grav.h>
#endif
#ifndef falcON_included_ionl_h
#  include <public/ionl.h>
#endif
////////////////////////////////////////////////////////////////////////////////

namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::stsp_leaf                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class stsp_leaf : public basic_leaf {
    stsp_leaf           (const stsp_leaf&);        // not implemented           
    stsp_leaf& operator=(const stsp_leaf&);        // not implemented           
    //--------------------------------------------------------------------------
    // data of class stsp_leaf (only static)                                    
    //--------------------------------------------------------------------------
  public:
    struct leaf_data {
      vect VEL;                                    // velocity / size^2         
    };
    //--------------------------------------------------------------------------
    // private data access                                                      
    //--------------------------------------------------------------------------
  private:
#define __data  static_cast<leaf_data*>(PROP)
    //--------------------------------------------------------------------------
    real           &size ()       { return SCAL; }
    uint           &num  ()       { return AUXU; }
    real           &sizeq()       { return __data->VEL[0]; }
    vect           &vel  ()       { return __data->VEL; }
    //--------------------------------------------------------------------------
    real      const&size () const { return SCAL; }
    uint      const&num  () const { return AUXU; }
    real      const&sizeq() const { return __data->VEL[0]; }
    vect      const&vel  () const { return __data->VEL; }
#undef __data
    //--------------------------------------------------------------------------
    // non-const methods                                                        
    //--------------------------------------------------------------------------
  public:
    void inc() { ++(num()); }
    void set_data(leaf_data*const&d) { PROP = static_cast<void*>(d); }
    //--------------------------------------------------------------------------
    // const data access via friends                                            
    //--------------------------------------------------------------------------
    friend uint const&mybody(const stsp_leaf*const&L) { return L->mybody(); } 
    friend uint const&num   (const stsp_leaf*const&L) { return L->num(); }
    friend vect const&pos   (const stsp_leaf*const&L) { return L->pos(); } 
    friend vect const&vel   (const stsp_leaf*const&L) { return L->vel(); } 
    friend real const&size  (const stsp_leaf*const&L) { return L->size(); } 
    friend real const&sizeq (const stsp_leaf*const&L) { return L->sizeq(); } 
    //--------------------------------------------------------------------------
    // copy data from body to leaf                                              
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void set_sticky(const bodies_type*const&B) {
      size() = B->size(mybody());
      vel () = B->vel (mybody());
      num () = 0u;
    };
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void set_sph(const bodies_type*const&B) {
      size () = B->size(mybody());
      sizeq() = square(size());
      num  () = 0u;
    };
    //--------------------------------------------------------------------------
    // copy data from leaf to body                                              
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void copy_to_bodies_num(const bodies_type*const&B) const {
      B->num(mybody()) = num();
    }
    //--------------------------------------------------------------------------
    // dump data                                                                
    //--------------------------------------------------------------------------
    static void dump_head(std::ostream& o) {
      basic_leaf::dump_head(o);
      o<<"           size [velocity]";
    }
    //--------------------------------------------------------------------------
    void dump(std::ostream &o) const
    {
      basic_leaf::dump(o);
      o<<' '<<setw(5)<<setprecision(4)<<size();
      if(this->is_set(flag::STICKY))
	for(register indx d=0; d!=Ndim; ++d)
	  o<<' '<<setw(7)<<setprecision(4)<<vel()[d];
    }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::stsp_cell                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class stsp_cell : public basic_cell {
    stsp_cell           (const stsp_cell&);        // not implemented           
    stsp_cell& operator=(const stsp_cell&);        // not implemented           
    //--------------------------------------------------------------------------
    // friendships                                                              
    //--------------------------------------------------------------------------
    friend class stsp_lister;                      // for alloc of srce_data    
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
  public:
    typedef stsp_leaf leaf_type;                   // type of associated leafs  
    //--------------------------------------------------------------------------
    // data of class cell                                                       
    //--------------------------------------------------------------------------
    struct srce_data {
      vect    POS;                                 // position center           
      vect    VEL;                                 // velocity center           
      real    SIZE;                                // size including leaf sizes 
      union {
	real  VRAD;                                // radius of position sphere 
	real  RMAX;                                // radius of velocity sphere 
      };
    };
#define __srce static_cast<srce_data*>(AUX1.PTER)
    //--------------------------------------------------------------------------
    // private data access                                                      
    //--------------------------------------------------------------------------
  private:
    vect      const&pos () const { return __srce->POS; }
    vect      const&vel () const { return __srce->VEL; }
    real      const&vrad() const { return __srce->VRAD; }
    real      const&rmax() const { return __srce->RMAX; }
    real      const&size() const { return __srce->SIZE; }
    uint      const&numb() const { return AUX2.NUMB; }
  public:
    //--------------------------------------------------------------------------
    void set_srce(srce_data*const&srce)
    {
      AUX1.PTER = static_cast<void*>(srce);
    }
    //--------------------------------------------------------------------------
    // non-const data access via members                                        
    //--------------------------------------------------------------------------
    vect&vel () { return __srce->VEL; }
    vect&pos () { return __srce->POS; }
    real&vrad() { return __srce->VRAD; }
    real&rmax() { return __srce->RMAX; }
    real&size() { return __srce->SIZE; }
    uint&numb() { return AUX2.NUMB; }
#undef __srce
    //--------------------------------------------------------------------------
    // const data access via friends                                            
    //--------------------------------------------------------------------------
    friend vect const&vel   (const stsp_cell*const&C) { return C->vel(); } 
    friend vect const&pos   (const stsp_cell*const&C) { return C->pos(); } 
    friend real const&size  (const stsp_cell*const&C) { return C->size(); } 
    friend real const&vrad  (const stsp_cell*const&C) { return C->vrad(); } 
    friend real const&rmax  (const stsp_cell*const&C) { return C->rmax(); } 
    friend uint const&numb  (const stsp_cell*const&C) { return C->numb(); } 
    //--------------------------------------------------------------------------
    // dump  data                                                               
    //--------------------------------------------------------------------------
    static void dump_head(std::ostream&o) {
      basic_cell::dump_head(o);
      o<<"           size         rmax/velocity            vrad]";
    }
    //--------------------------------------------------------------------------
    void dump(std::ostream &o) const
    {
      basic_cell::dump(o);
      o<<' '<<setw(6)<<setprecision(4)<<size();
      if       (this->is_set(flag::STICKY)) {
	for(register indx d=0; d!=Ndim; ++d)
	  o<<' '<<setw(9)<<setprecision(4)<<vel()[d];
	o<<' '<<setw(5)<<setprecision(4)<<vrad();
      } else if(this->is_set(flag::SPH))
	o<<' '<<setw(5)<<setprecision(4)<<rmax();
    }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::stsp_estimator                                                 
  //                                                                            
  // NOTE on the meaning of the activity and stsp flags                         
  //                                                                            
  // Each cell knows both the total number of leafs (basic_cell::number()) and  
  // the number of stsp cells (stsp_cell::numb()).                              
  // If    numb(stsp_cell*) == 0                   then  is_sph(stsp_cell*)==0  
  // If    numb(stsp_cell*) == number(stsp_cell*)  then  al_sph(stsp_cell*)==1  
  //                                                                            
  // The activity flag only refers to the stsp leafs, not all leafs:            
  // If # active stsp leaf descendants == 0      then  is_active()==0           
  // If # active stsp leaf descendants == numb() then  al_active()==1           
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class stsp_estimator {
    stsp_estimator           (const stsp_estimator&);
    stsp_estimator& operator=(const stsp_estimator&);
    //--------------------------------------------------------------------------
    // data:                                                                    
    //--------------------------------------------------------------------------
  private:
    const oct_tree       *TREE;                    // the tree to be used       
    stsp_leaf::leaf_data *LEAF_DATA;               // memory for leafs          
    stsp_cell::srce_data *CELL_SRCE;               // memory for cell srce      
    mutable bool          ALL_STSP;                // all leafs are stsp        
    mutable bool          ALL_ACTIVE;              // all stsp leafs are active 
    mutable bool          SPH_UPTODATE;            // tree ready for sph search 
    mutable bool          STC_UPTODATE;            // tree ready for sticky --  
    mutable uint          NL,NC;                   // # stsp leafs & cells      
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void copy_to_bodies_num (const bodies_type*const&) const;
    template<typename bodies_type>
    void update_leafs_sph   (const bodies_type*const&);
//     void update_leafs_sph   (const abodies*const&);
    //--------------------------------------------------------------------------
    void update_leafs_sticky();
    void update_leafs_sph   ();
    void prepare_sticky     ();
    void prepare_sph        ();
    //--------------------------------------------------------------------------
    // tree stuff to be superseeded                                             
    //--------------------------------------------------------------------------
  public:
    typedef stsp_cell                       cell_type;
    typedef stsp_leaf                       leaf_type;
    typedef oct_tree::CellIter<stsp_cell>   cell_iterator;
    typedef leaf_type*                      leaf_iterator;
    //--------------------------------------------------------------------------
    const oct_tree*const&my_tree() const { return TREE; }
    cell_iterator root          () const {
      return cell_iterator(TREE,static_cast<stsp_cell*>(TREE->FstCell())); }
    //--------------------------------------------------------------------------
    // dump cell and leaf data                                                  
    //--------------------------------------------------------------------------
    void dump_cells(std::ostream&) const;
    void dump_leafs(std::ostream&) const;
    //--------------------------------------------------------------------------
    // public type                                                              
    //--------------------------------------------------------------------------
    typedef uint elem_pair[2];                     // element: interaction list 
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
    stsp_estimator(const oct_tree*const&T) :       // I: tree to be used        
      TREE         ( T ),
      LEAF_DATA    ( 0 ),
      CELL_SRCE    ( 0 ),
      ALL_STSP     ( 0 ),
      ALL_ACTIVE   ( 0 ),
      SPH_UPTODATE ( 0 ),
      STC_UPTODATE ( 0 ) {}
    //--------------------------------------------------------------------------
    void reset() {
      if(CELL_SRCE) { delete[] CELL_SRCE; CELL_SRCE=0; }
      if(LEAF_DATA) { delete[] LEAF_DATA; LEAF_DATA=0; }
      SPH_UPTODATE = 0;
      STC_UPTODATE = 0;
    }
    //--------------------------------------------------------------------------
    void new_tree(const oct_tree*const&T) {
      TREE = T;
      reset();
    }
    //--------------------------------------------------------------------------
    // destruction                                                              
    //--------------------------------------------------------------------------
    ~stsp_estimator() {
      if(CELL_SRCE) delete[] CELL_SRCE;
      if(LEAF_DATA) delete[] LEAF_DATA;
    }
    //--------------------------------------------------------------------------
    void make_sticky_list (elem_pair *,            // I/O: interaction list     
			   uint const&,            // I: physical size of list  
			   uint      &,            // O: # pairs found          
			   real const&);           // I: tau                    
    //--------------------------------------------------------------------------
    void make_sph_list    (elem_pair *,            // I/O: interaction list     
			   uint const&,            // I: physical size of list  
			   uint      &,            // O: # pairs found          
			   bool const&);           // I: r < max(hi,hj) or hi+hj
    //--------------------------------------------------------------------------
    void count_sph_partners();                     // counter sph iaction ptners
    //--------------------------------------------------------------------------
  };
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_stic_h
