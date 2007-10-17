// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/manip/symmetrize_pairs.cc                                
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2004-2007                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2004-2007 Walter Dehnen                                        
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
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>

namespace falcON { namespace Manipulate {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class symmetrize_pairs                                                     
  //                                                                            
  /// manipulator: symmetrizes vectors wrt origin for pairs of bodies           
  ///                                                                           
  /// This manipulator symmetrizes all vector properties (position, velocity,   
  /// acceleration, jerk, auxiliary vector, predicted velocity) for pairs of    
  /// adjacent bodies (bodies 0 and 1 are the first pair) of which the          
  /// first is in_subset() (default: all, see set_subset), so that the          
  /// sum of the positions of the two bodies in each pair is zero (and          
  /// the same for velocities etc).                                             
  /// \note we also ensure that the masses of each of a pair of bodies are      
  /// equal.                                                                    
  ///                                                                           
  /// No parameters or file used.                                               
  ///                                                                           
  /// Usage of pointers: none\n                                                 
  /// Usage of flags:    uses in_subset()\n                                     
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class symmetrize_pairs : public manipulator {
  private:
    //--------------------------------------------------------------------------
    static void make_equal(real&m1, real&m2) {
      m1 += m2;
      m1 *= half;
      m2  = m1;
    }
    //--------------------------------------------------------------------------
    static void ensure_equal_masses(const bodies*B) {
      for(body b(B->begin_all_bodies()); b; ++b) 
	if(in_subset(b)) {
	  real &m(b.mass());
	  if(++b) make_equal(m, b.mass());
	} else
	  ++b;
    }
    //--------------------------------------------------------------------------
    static void make_symmetric(vect&x1, vect&x2) {
      x1 -= x2;
      x1 *= half;
      x2  =-x1;
    }
    //--------------------------------------------------------------------------
    template<int BIT>
    static void symmetrize(const bodies*B) {
      for(body b(B->begin_all_bodies()); b; ++b) 
	if(in_subset(b)) {
	  typename field_traits<BIT>::type &x(b.template datum<BIT>());
	  if(++b) make_symmetric(x, b. template datum<BIT>());
	} else
	  ++b;
    }
  public:
    //--------------------------------------------------------------------------
    const char* name    () const { return "symmetrize_pairs"; }
    const char* describe() const {
      return
	"symmetrize vectors w.r.t. origin for pairs of bodies "
	"passing 'filter' (default: all)";
    }
    //--------------------------------------------------------------------------
    fieldset          need    () const { return fieldset::o; }
    fieldset          provide () const { return fieldset::o; }
    fieldset          change  () const { return fieldset::vectors; }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    symmetrize_pairs(const char*pars,
		     const char*file)
    {
      if(pars && debug(1) || debug(2))
	std::cerr<<
	  " Manipulator \""<<name()<<"\":\n"
	  " symmetrize  vectors  w.r.t. origin"
	  " for pairs of bodies passing 'filter' (default: all)\n";
      if(pars  && debug(1))
	warning(" Manipulator \"%s\": skipping all parameters\n",name());
    }
    //--------------------------------------------------------------------------
    ~symmetrize_pairs() {}
  };
  //////////////////////////////////////////////////////////////////////////////
  bool symmetrize_pairs::manipulate(const snapshot*S) const {
    if(S->have(fieldbit::m)) ensure_equal_masses(S);
    if(S->have(fieldbit::x)) symmetrize<fieldbit::x>(S);
    if(S->have(fieldbit::v)) symmetrize<fieldbit::v>(S);
    if(S->have(fieldbit::a)) symmetrize<fieldbit::a>(S);
    if(S->have(fieldbit::j)) symmetrize<fieldbit::j>(S);
    if(S->have(fieldbit::w)) symmetrize<fieldbit::w>(S);
    if(S->have(fieldbit::z)) symmetrize<fieldbit::z>(S);
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} }

__DEF__MAN__ALT(falcON::Manipulate::symmetrize_pairs)
