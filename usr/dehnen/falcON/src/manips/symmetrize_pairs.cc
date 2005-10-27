// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// symmetrize_pairs.cc                                                         |
//                                                                             |
// Copyright (C) 2004, 2005 Walter Dehnen                                      |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/defman.h>

namespace {
  using namespace falcON;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class symmetrize_pairs                                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class symmetrize_pairs : public manipulator {
  private:
    const int K, N;
    //--------------------------------------------------------------------------
    static void make_symmetric(vect&x1, vect&x2) {
      x1 -= x2;
      x1 *= half;
      x2  =-x1;
    }
    //--------------------------------------------------------------------------
    template<int BIT>
    static void symmetrize(body const&bMin, body const&bMax) {
      for(body b(bMin); ; ) {
	typename field_traits<BIT>::type &x(b. template datum<BIT>());
	if(++b == bMax) break;
	make_symmetric(x, b. template datum<BIT>());
	if(++b == bMax) break;
      }
    }
  public:
    //--------------------------------------------------------------------------
    const char* name    () const { return "symmetrize_pairs"; }
    const char* describe() const {
      if(N>0) return message("symmetrize pairs of bodies [%d:%d] w.r.t. origin",
			     K,N);
      if(K>0) return message("symmetrize pairs of bodies [i>=%d] w.r.t. origin",
			     K);
      else    return "symmetrize pairs of bodies w.r.t. origin";
    }
    //--------------------------------------------------------------------------
    fieldset          need    () const { return fieldset::o; }
    fieldset          provide () const { return fieldset::o; }
    fieldset          change  () const { return fieldset::vectors; }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    symmetrize_pairs(const double*pars,
		     int          npar,
		     const char  *file) :
      K   (2* ((npar>0? int (pars[0]) : 0)/2)),
      N   (2* ((npar>1? int (pars[1]) : 0)/2))
    {
      if(K >= N)
	error(" Manipulator \"symmetrize_pairs\":"
	      " par[0]=%d >= par[1]=%d",K,N);
      if(npar<2 && debug(1) || debug(2))
	std::cerr<<
	  " Manipulator \"symmetrize_pairs\":\n"
	  " symmetrize  xvaw  w.r.t. origin"
	  " for pairs of bodies with k in [ par[0],par[1] [\n"
	  " Default: all\n";
      if(npar>2  && debug(1))
	warning(" Manipulator \"symmetrize_pairs\":"
		" skipping parameters beyond 2\n");
    }
    //--------------------------------------------------------------------------
    ~symmetrize_pairs() {}
  };
  //////////////////////////////////////////////////////////////////////////////
  bool symmetrize_pairs::manipulate(const snapshot*S) const {
    if(S->have_not(fieldset::vectors)) return false;
    const body bMin = S->bodyNo(K>0? K : 0);
    const body bMax = S->bodyNo(N<=0?            S->N_bodies() : 
				N<S->N_bodies()? N : S->N_bodies() );
    if(S->have(fieldbit::x)) symmetrize<fieldbit::x>(bMin,bMax);
    if(S->have(fieldbit::v)) symmetrize<fieldbit::v>(bMin,bMax);
    if(S->have(fieldbit::a)) symmetrize<fieldbit::a>(bMin,bMax);
    if(S->have(fieldbit::j)) symmetrize<fieldbit::j>(bMin,bMax);
    if(S->have(fieldbit::w)) symmetrize<fieldbit::w>(bMin,bMax);
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} // namespace{

__DEF__MAN(symmetrize_pairs)
