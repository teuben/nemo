// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tool.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_tool_h
#define falcON_included_tool_h 1

#ifndef falcON_included_body_h
#  include <body.h>
#endif
#ifndef falcON_included_tree_h
#  include <public/tree.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //----------------------------------------------------------------------------
  // compute centre of mass of bodies                                           
  //----------------------------------------------------------------------------
  vect centre_of_mass(const bodies*const&);        // I: bodies                 
  //----------------------------------------------------------------------------
  inline                                           // american spelling         
  vect center_of_mass(const bodies*const&B)
  {
    return centre_of_mass(B);
  }
  //----------------------------------------------------------------------------
  // compute iteratively the |Phi|^alpha weighted centre                        
  //                                                                            
  //   The centre is found iteratively as the m*|phi|^alpha weighted position   
  //   of all bodies that are within r_cut from the last centre position. At    
  //   each iteration r_cut is reduce by the factor f. The iteration is stopped 
  //   if the number of bodies in the final sphere is the smallest greater than 
  //   Nmin.                                                                    
  //   If the bool argument is true, we use the input values for centre position
  //   and radius as starting point for the iteration. Otherwise, we use the    
  //   mean position and mean-square radius of the whole distribution.          
  //----------------------------------------------------------------------------
  void find_centre(const bodies*const&,            // I  : bodies               
		   real         const&,            // I  : reduction factor f   
		   uint         const&,            // I  : alpha                
		   uint         const&,            // I  : Nmin                 
		   vect              &,            // I/O: centre position      
		   real              &,            // I/O: centre radius        
		   bool         const& = 0,        //[I  : use input as initial 
		   vect              * = 0,        //[O  : centre velocity]     
		   real              * = 0,        //[O  : centre density]      
		   uint         const& = 0,        //[I  : begin of bodies]     
		   uint         const& = 0);       //[I  : end   of bodies]     
  //----------------------------------------------------------------------------
  inline                                           // american spelling         
  void find_center(const bodies*const&b,           // I  : bodies               
		   real         const&f,           // I  : reduction factor f   
		   uint         const&a,           // I  : alpha                
		   uint         const&n,           // I  : Nmin                 
		   vect              &x,           // I/O: centre position      
		   real              &r,           // I/O: centre radius        
		   bool         const&i = 0,       //[I  : use input as initial 
		   vect              *v = 0,       //[O  : centre velocity]     
		   real              *d = 0,       //[O  : centre density]      
		   uint         const&b0= 0,       //[I  : begin of bodies]     
		   uint         const&bn= 0)       //[I  : end   of bodies]     
  {
    find_centre(b,f,a,n,x,r,i,v,d,b0,bn);
  }
  //----------------------------------------------------------------------------
  // compute iteratively the density centre                                     
  //                                                                            
  //----------------------------------------------------------------------------
  void find_centre(const bodies*const&,            // I  : bodies               
		   uint         const&,            // I  : Nmin                 
		   vect              &,            // I/O: centre position      
		   real              &,            // I/O: centre radius        
		   vect              * = 0,        //[O  : centre velocity]     
		   real              * = 0,        //[O  : centre density]      
		   uint         const& = 0,        //[I  : begin of bodies]     
		   uint         const& = 0);       //[I  : end   of bodies]     
  //----------------------------------------------------------------------------
  inline                                           // american spelling         
  void find_center(const bodies*const&b,           // I  : bodies               
		   uint         const&n,           // I  : Nmin                 
		   vect              &x,           // I/O: centre position      
		   real              &r,           // I/O: centre radius        
		   vect              *v = 0,       //[O  : centre velocity]     
		   real              *d = 0,       //[O  : centre density]      
		   uint         const&b0= 0,       //[I  : begin of bodies]     
		   uint         const&bn= 0)       //[I  : end   of bodies]     
  {
    find_centre(b,n,x,r,v,d,b0,bn);
  }
  //----------------------------------------------------------------------------
  // estimate position of density maximum weighted by w = mass * |pot|^alpha    
  //     given an oct_tree which holds the bodies we                            
  //     - find the cell with the maximum weight density                        
  //     - compute its weighted center                                          
  //     - estimate h (center radius) to yield Nmin bodies                      
  //----------------------------------------------------------------------------
  void estimate_density_peak(oct_tree*const&,      // I  : oct_tree             
			     uint     const&,      // I  : alpha                
			     uint     const&,      // I  : Nmin                 
			     vect          &,      // O  : center position      
			     real          &);     // O  : center radius        
  //----------------------------------------------------------------------------
  // find the lagrange radii around the origin for a given set of cumulative    
  // masses (relative to total).                                                
  //----------------------------------------------------------------------------
  void find_lagrange_rad(const bodies*const&,      // I: bodies                 
			       int    const&,      // I: size of tables         
			 const double*const&,      // I: table: masses          
			       double*const&);     // O: table: lagrange radii  
  //----------------------------------------------------------------------------
}
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_tool_h
