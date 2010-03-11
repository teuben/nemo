// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   src/public/manip/centre_of_mass.cc
///
/// \author Walter Dehnen
/// \date   2010
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Walter Dehnen
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
//
// history:
//
// v 0.0    09/03/2010  WD created
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>

namespace falcON {
namespace Manipulate {
  //
  // class centre_of_mass
  //
  /// manipulator: sets the centre to centre of mass of subset
  ///                                                                       
  /// Usage of parameters: none\n
  /// Usage of file:       if given, report centre
  /// Usage of pointers:   sets 'xcen', 'vcen', and 'mcen'\n
  /// Usage of flags:      uses in_subset()
  ///
  class centre_of_mass : public manipulator {
  private:
    mutable output OUT;
    mutable vect   XCEN,VCEN;
    mutable real   MCEN;
    mutable bool   FIRST;
    //
    static std::ostream&print_line(std::ostream&to) {
      return to << 
	"#-------------"
	"----------------"
	"----------------"
	"----------------"
	"----------------"
	"----------------"
	"----------------\n";
    }
    static std::ostream&print_head(std::ostream&to) {
      return to << 
	"#        time "
	"              x "
	"              y "
	"              z "
	"            v_x "
	"            v_y "
	"            v_z\n";
    }
    //
  public:
    const char*name    () const
    { return "centre_of_mass"; }
    const char*describe() const
    {
      return 
	"provides 'xcen', 'vcen', and 'mcen' as position, velocity, and mass "
	"of the centre of mass of subset (default: all)";
    }
    //
    fieldset need   () const { return fieldset::basic; }
    fieldset provide() const { return fieldset::empty; }
    fieldset change () const { return fieldset::empty; }
    //
    centre_of_mass(const double*, int npar, const char*file) falcON_THROWING
    : OUT  ( file, true ),
      XCEN ( vect(zero) ),
      VCEN ( vect(zero) ),
      FIRST( true )
    {
      if(debug(2) || npar)
	std::cerr <<
	  " Manipulator \"centre_of_mass\":\n"
	  " find centre of mass of subset and put it in 'xcen' and 'vcen'\n";
      if(npar)
	falcON_WarningN("Manipulator \"%s\": ignoring parameters\n",name());
    }
    //
    bool manipulate(const snapshot*) const;
  };
  //
  bool centre_of_mass::manipulate(const snapshot*S) const
  {
    // 0  first call: write header to output file, if any
    if(FIRST && OUT) {
      FIRST = false;
      print_line(OUT);
      OUT  << "#\n# output from Manipulator \"" << name() << "\"\n#\n";
      RunInfo::header(OUT);
      print_head(OUT);
      print_line(OUT);
    }
    // 1  find centre of mass of subset
    double M(0);
    vect_d X(0.), V(0.);
    if(S->have_flag())
      LoopSubsetBodies(S,b) {
	M += mass(b);
	X += mass(b) * pos(b);
	V += mass(b) * vel(b);
      }
    else
      LoopAllBodies(S,b) {
	M += mass(b);
	X += mass(b) * pos(b);
	V += mass(b) * vel(b);
      }
    MCEN = M;
    X /= M; XCEN=X;
    V /= M; VCEN=V;
    // 2  output
    if(OUT)
      OUT << ' '
	  << print(S->time(),12,8) << ' '
	  << print(XCEN,15,8) << ' '
	  << print(VCEN,15,8) << std::endl;
    // 3  put centre under 'xcen' and 'vcen'
    S->set_pointer(&XCEN,"xcen");
    S->set_pointer(&VCEN,"vcen");
    S->set_pointer(&MCEN,"mcen");
    return false;
  }
  //
} }

__DEF__MAN(falcON::Manipulate::centre_of_mass)
