// -*- C++ -*-                                                                  
//------------------------------------------------------------------------------
//                                                                              
// inc/public/partner.h                                                         
//                                                                              
// Copyright (C) 2000-2005  Walter Dehnen                                       
//                                                                              
// This program is free software; you can redistribute it and/or modify         
// it under the terms of the GNU General Public License as published by         
// the Free Software Foundation; either version 2 of the License, or (at        
// your option) any later version.                                              
//                                                                              
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// General Public License for more details.                                     
//                                                                              
// You should have received a copy of the GNU General Public License            
// along with this program; if not, write to the Free Software                  
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                    
//                                                                              
//------------------------------------------------------------------------------
#ifndef falcON_included_partner_h
#define falcON_included_partner_h

#ifndef falcON_included_tree_h
#  include <public/tree.h>
#endif
#ifndef falcON_included_inline_io_h
#  include <public/inline_io.h>
#endif
////////////////////////////////////////////////////////////////////////////////

namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::PartnerLeaf                                                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class PartnerLeaf : public BasicLeaf {
    PartnerLeaf           (const PartnerLeaf&);    // not implemented           
    PartnerLeaf& operator=(const PartnerLeaf&);    // not implemented           
    //--------------------------------------------------------------------------
    // data of class PartnerLeaf (only static)                                  
    //--------------------------------------------------------------------------
  public:
    struct leaf_data {
      vect VEL;                                    // velocity / size^2         
    };
    //--------------------------------------------------------------------------
    // private data access                                                      
    //--------------------------------------------------------------------------
  private:
    //--------------------------------------------------------------------------
    real      &size ()       { return SCAL; }
    unsigned  &num  ()       { return AUXU; }
    real      &sizeq()       { return static_cast<leaf_data*>(PROP)->VEL[0]; }
    vect      &vel  ()       { return static_cast<leaf_data*>(PROP)->VEL; }
    //--------------------------------------------------------------------------
    real const&size () const { return SCAL; }
    unsigned const&num  () const { return AUXU; }
    real const&sizeq() const { return static_cast<leaf_data*>(PROP)->VEL[0]; }
    vect const&vel  () const { return static_cast<leaf_data*>(PROP)->VEL; }
    //--------------------------------------------------------------------------
    // non-const methods                                                        
    //--------------------------------------------------------------------------
  public:
    void inc() { ++(num()); }
    void set_data(leaf_data*const&d) { PROP = static_cast<void*>(d); }
    //--------------------------------------------------------------------------
    // const data access via friends                                            
    //--------------------------------------------------------------------------
    friend bodies::index const&mybody(const PartnerLeaf*L) {
      return L->mybody(); } 
    friend unsigned const&num   (const PartnerLeaf*L) { return L->num(); }
    friend vect     const&pos   (const PartnerLeaf*L) { return L->pos(); } 
    friend vect     const&vel   (const PartnerLeaf*L) { return L->vel(); } 
    friend real     const&size  (const PartnerLeaf*L) { return L->size(); } 
    friend real     const&sizeq (const PartnerLeaf*L) { return L->sizeq(); } 
    //--------------------------------------------------------------------------
    // copy data from body to leaf                                              
    //--------------------------------------------------------------------------
    void set_sticky(const bodies*const&B) {
      size() = B->size(mybody());
      vel () = B->vel (mybody());
      num () = 0u;
    }
    //--------------------------------------------------------------------------
    void set_sph(const bodies*const&B) {
      size () = B->size(mybody());
      sizeq() = square(size());
      num  () = 0u;
    }
    //--------------------------------------------------------------------------
    // copy data from leaf to body                                              
    //--------------------------------------------------------------------------
    void copy_to_bodies_num(const bodies*const&B) const {
      B->num(mybody()) = num();
    }
    //--------------------------------------------------------------------------
    // dump data                                                                
    //--------------------------------------------------------------------------
    static void dump_head(std::ostream& o) {
      BasicLeaf::dump_head(o);
      o<<"           size [velocity]";
    }
    //--------------------------------------------------------------------------
    void dump(std::ostream &o) const
    {
      BasicLeaf::dump(o);
      o<<' '<<setw(5)<<setprecision(4)<<size();
      if(flg().is_set(flag::STICKY))
	for(register indx d=0; d!=Ndim; ++d)
	  o<<' '<<setw(7)<<setprecision(4)<<vel()[d];
    }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::PartnerCell                                                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class PartnerCell : public BasicCell {
    PartnerCell           (const PartnerCell&);    // not implemented           
    PartnerCell& operator=(const PartnerCell&);    // not implemented           
    //--------------------------------------------------------------------------
    // friendships                                                              
    //--------------------------------------------------------------------------
    friend class stsp_lister;                      // for alloc of srce_data    
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
  public:
    typedef PartnerLeaf leaf_type;                 // type of associated leafs  
    //--------------------------------------------------------------------------
    // data of class cell                                                       
    //--------------------------------------------------------------------------
    struct srce_data {
      vect    VEL;                                 // velocity center           
    };
    //--------------------------------------------------------------------------
    // private data access                                                      
    //--------------------------------------------------------------------------
  private:
    vect     const&pos () const { return POS; }
    vect     const&vel () const { return 
	static_cast<srce_data*>(AUX1.PTER)->VEL; }
    real     const&vrad() const { return RAD; }
    real     const&rmax() const { return RAD; }
    real     const&size() const { return AUX3.SCAL; }
    unsigned const&numb() const { return AUX2.NUMB; }
  public:
    //--------------------------------------------------------------------------
    void set_srce(srce_data*const&srce)
    {
      AUX1.PTER = static_cast<void*>(srce);
    }
    //--------------------------------------------------------------------------
    // non-const data access via members                                        
    //--------------------------------------------------------------------------
    vect    &vel () { return static_cast<srce_data*>(AUX1.PTER)->VEL; }
    vect    &pos () { return POS; }
    real    &vrad() { return RAD; }
    real    &rmax() { return RAD; }
    real    &size() { return AUX3.SCAL; }
    unsigned&numb() { return AUX2.NUMB; }
    //--------------------------------------------------------------------------
    // const data access via friends                                            
    //--------------------------------------------------------------------------
    friend vect     const&vel   (const PartnerCell*C) { return C->vel(); } 
    friend vect     const&pos   (const PartnerCell*C) { return C->pos(); } 
    friend real     const&size  (const PartnerCell*C) { return C->size(); } 
    friend real     const&vrad  (const PartnerCell*C) { return C->vrad(); } 
    friend real     const&rmax  (const PartnerCell*C) { return C->rmax(); } 
    friend unsigned const&numb  (const PartnerCell*C) { return C->numb(); } 
    //--------------------------------------------------------------------------
    // dump  data                                                               
    //--------------------------------------------------------------------------
    static void dump_head(std::ostream&o) {
      BasicCell::dump_head(o);
      o<<"           size         rmax/velocity            vrad]";
    }
    //--------------------------------------------------------------------------
    void dump(std::ostream &o) const
    {
      BasicCell::dump(o);
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
  // class falcON::PartnerEstimator                                             
  //                                                                            
  // NOTE on the meaning of the activity and stsp flags                         
  //                                                                            
  // Each cell knows both the total number of leafs (BasicCell::number()) and   
  // the number of SPH/sticky leafs (PartnerCell::numb()).                      
  // If      numb(PartnerCell*) == 0                                            
  //   then  is_sph(PartnerCell*)==0                                            
  // If      numb(PartnerCell*) == number(PartnerCell*)                         
  //   then  al_sph(PartnerCell*)==1                                            
  //                                                                            
  // The activity flag only refers to the stsp leafs, not all leafs:            
  // If # active stsp leaf descendants == 0      then  is_active()==0           
  // If # active stsp leaf descendants == numb() then  al_active()==1           
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class PartnerEstimator {
    PartnerEstimator           (const PartnerEstimator&);
    PartnerEstimator& operator=(const PartnerEstimator&);
    //--------------------------------------------------------------------------
    // public type                                                              
    //--------------------------------------------------------------------------
  public:
    typedef bodies::index indx_pair[2];            // element: interaction list 
    //--------------------------------------------------------------------------
    // data:                                                                    
    //--------------------------------------------------------------------------
  private:
    const OctTree            *TREE;                // the tree to be used       
    PartnerLeaf::leaf_data   *LEAF_DATA;           // memory for leafs          
    PartnerCell::srce_data   *CELL_SRCE;           // memory for cell srce      
    mutable bool              ALL_STSP;            // all leafs are stsp        
    mutable bool              ALL_ACTIVE;          // all stsp leafs are active 
    mutable bool              SPH_UPTODATE;        // tree ready for sph search 
    mutable bool              STC_UPTODATE;        // tree ready for sticky --  
    mutable unsigned          NL,NC;               // # stsp leafs & cells      
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
    void copy_to_bodies_num (bool) const;
    void update_leafs_sph   (const bodies*const&);
    //--------------------------------------------------------------------------
    void update_leafs_sticky();
    void update_leafs_sph   ();
    void prepare_sticky     ();
    void prepare_sph        ();
    //--------------------------------------------------------------------------
    template<bool> void make_st_list (indx_pair*, unsigned, unsigned&, real);
    template<bool> void make_sp_list (indx_pair*, unsigned, unsigned&, bool);
    //--------------------------------------------------------------------------
    // tree stuff to be superseeded                                             
    //--------------------------------------------------------------------------
  public:
    typedef PartnerCell                         cell_type;
    typedef PartnerLeaf                         leaf_type;
    typedef OctTree::CellIter<PartnerCell>      cell_iterator;
    typedef leaf_type*                          leaf_iterator;
    //--------------------------------------------------------------------------
    const OctTree*const&my_tree() const { return TREE; }
    cell_iterator root         () const {
      return cell_iterator(TREE,static_cast<PartnerCell*>(TREE->FstCell())); }
    //--------------------------------------------------------------------------
    // dump cell and leaf data                                                  
    //--------------------------------------------------------------------------
    void dump_cells(std::ostream&) const;
    void dump_leafs(std::ostream&) const;
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
    PartnerEstimator(const OctTree*const&T) :      // I: tree to be used        
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
    void new_tree(const OctTree*const&T) {
      TREE = T;
      reset();
    }
    //--------------------------------------------------------------------------
    // destruction                                                              
    //--------------------------------------------------------------------------
    ~PartnerEstimator() {
      if(CELL_SRCE) delete[] CELL_SRCE;
      if(LEAF_DATA) delete[] LEAF_DATA;
    }
    //--------------------------------------------------------------------------
    void make_sticky_list (indx_pair*,             // I/O: interaction list     
			   unsigned  ,             // I: physical size of list  
			   unsigned &,             // O: # pairs found          
			   real      ,             // I: tau                    
			   bool      )             // I: count as well?         
      falcON_THROWING;
    //--------------------------------------------------------------------------
    void make_sph_list    (indx_pair*,             // I/O: interaction list     
			   unsigned  ,             // I: physical size of list  
			   unsigned &,             // O: # pairs found          
			   bool      ,             // I: r < max(hi,hj) or hi+hj
			   bool      )             // I: count as well?         
      falcON_THROWING;
    //--------------------------------------------------------------------------
    void count_sph_partners(                       // count sph iaction partners
			    bool) falcON_THROWING; // I: r < max(hi,hj) or hi+hj
    //--------------------------------------------------------------------------
  };
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_partner_h  
