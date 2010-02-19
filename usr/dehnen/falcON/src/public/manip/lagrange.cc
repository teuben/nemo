// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   src/public/manip/lagrange.cc
///
/// \author Walter Dehnen
/// \date   2004-2008
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004-2008 Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
//
// history:                                                                     
//
// v 1.0    27/06/2006  WD using bodyset in 'subset' 
// v 1.1    27/06/2006  WD using filter in 'filter' instead of 'subset' 
// v 1.2    07/07/2006  WD using flags::ignore instead of filter
// v 1.3    11/09/2008  WD erased direct use of nemo functions
// v 1.4    03/11/2008  WD appending, print(), RunInfo::header()
// v 1.5    07/11/2008  WD using FindPercentile from numerics.h, sanity checks
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <utils/heap.h>
#include <utils/numerics.h>
#include <ctime>
#include <cmath>

namespace falcON {
namespace Manipulate {
  // ///////////////////////////////////////////////////////////////////////////
  //
  // class lagrange
  //
  /// manipulator: computes Lagrange radii for subset, writes them to file
  ///
  /// This manipulator computes the Lagrange radii w.r.t. centre 'xcen' for
  /// all bodies in_subset() (default: all, see set_subset), and writes them
  /// to file.
  ///
  /// Meaning of the parameters:\n
  /// par[0-n]: mass fractions for which Lagrange radii shall be computed\n
  /// file:     file name for output of table with Lagrange radii\n
  ///
  /// Usage of pointers: uses 'xcen'\n
  /// Usage of flags:    uses in_subset()\n
  ///
  // ///////////////////////////////////////////////////////////////////////////
  class lagrange : public manipulator {
  private:
    int              N;    ///< number of mass fractions
    double          *F;    ///< table with mass fractions
    mutable output   OUT;  ///< output
    mutable bool     FST;  ///< true before first manipulation
    //--------------------------------------------------------------------------
    void print_line() const {
      OUT  << "#-----------";
      for(int i=0; i!=N; ++i)
	OUT<< "---------";
      OUT  << '\n';
    }
    //--------------------------------------------------------------------------
    void print_head() const {
      OUT  << "# time      ";
      for(int i=0; i!=N; ++i)
	OUT<< " r[" << std::setw(4) << 100*F[i] << "%]";
      OUT  << '\n';
    }
    //--------------------------------------------------------------------------
  public:
    const char* name    () const { return "lagrange"; }
    const char* describe() const {
      return
	"computes lagrange radii of bodies passing 'filter' (default: all) "
	"w.r.t. \"xcen\" (default: origin)";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::m | fieldset::x; }
    fieldset provide() const { return fieldset::empty; }
    fieldset change () const { return fieldset::empty; }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    lagrange(const double*pars,
	     int          npar,
	     const char  *file) :
      N   ( npar ),
      F   ( npar>0? falcON_NEW(double,npar) : 0 ),
      OUT ( file? file : ".", true),
      FST ( true )
    {
      if(debug(2) || file == 0 || npar == 0) {
	std::cerr<<
	  " Manipulator \""<<name()<<"\":\n"
	  " computes lagrange radii w.r.t. 'xcen' (default: origin)\n"
	  " at given mass fractions (parameters) for bodies in_subset()\n"
	  " (default: all) and writes them given data file.\n";
      }
      // 1 ensure that that output is okay
      if(file == 0)
	falcON_ErrorN("Manipulator \"%s\": no output file given\n",name());
      if(!OUT.is_open())
	falcON_ErrorN("Manipulator \"%s\": couldn't open output\n",name());
      // 2 ensure that mass fractions are okay
      if(npar == 0)
	falcON_ErrorN("Manipulator \"%s\": no mass fractions given\n",name());
      // 2.1 ensure that mass are fractions in ascending order
      bool okay=true;
      for(int i=0; i!=N; ++i) {
	F[i] = pars[i];
	if(i && F[i-1] > F[i])
	  okay=false;
      }
      if(!okay) {
	falcON_WarningN("Manipulator \"%s\": "
			"mass fractions not in ascending order: "
			"will order them.\n",name());
	HeapSortAsc(F,N);
      }
      // 2.2 ensure that mass fractions are distinct and in range [0,1[
      int Nf=0;
      for(int i=0; i!=N; ++i)
	if(F[i]>0 && F[i]<1 && (i==0 || F[i]>F[i-1])) ++Nf;
      if(Nf==0)
	falcON_ErrorN("Manipulator \"%s\": no mass fraction in range ]0,1[\n",
		      name());
      if(Nf <N) {
	falcON_WarningN("Manipulator \"%s\": "
			"only %d of %d mass fractions given are in range "
			"]0,1[ and distinct: will only use these.\n",
			name(),Nf,N);
	double *newF = falcON_NEW(double,Nf);
	for(int i=0,j=0; i!=N; ++i)
	  if(F[i]>0 && F[i]<1 && (i==0 || F[i]>F[i-1]))
	    newF[j++] = F[i];
	falcON_DEL_A(F); F=newF;
	N=Nf;
      }
    }
    //--------------------------------------------------------------------------
    ~lagrange() {
      if(F) { falcON_DEL_A(F); F=0; }
      N=0;
    }
  };
} }
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace falcON;
  vect X;
  body B;
  void SetAll(unsigned, double&q, double&m)
  {
    if(B) {
      q = dist_sq(X,pos(B));
      m = mass(B);
      ++B;
    } else
      falcON_ErrorN("Manipulator \"lagrange\": body number mismatch\n");
  }
  void SetSub(unsigned, double&q, double&m)
  {
    if(B) {
      q = dist_sq(X,pos(B));
      m = mass(B);
      B.next_in_subset();
    } else
      falcON_ErrorN("Manipulator \"lagrange\": body number mismatch\n");
  }
}
//
namespace falcON {
namespace Manipulate {
  bool lagrange::manipulate(const snapshot*S) const
  {
    // 0 check for required data
    if( !CheckMissingBodyData(S,need()) )
      return false;
    // 1 first ever call: print header
    if(FST) {
      FST = false;
      print_line();
      OUT  << "#\n# output from Manipulator \"lagrange\"\n#\n";
      RunInfo::header(OUT);
      print_head();
      print_line();
      OUT.stream().setf(std::ios::left, std::ios::adjustfield);
    }
    DebugInfo(4,"Manipulator \"lagrange\": working out lagrange radii ...\n");
    // 2   find lagrange radii w.r.t. centre and write them to file.
    // 2.1 get centre
    const vect*pX = S->pointer<vect>("xcen");
    if(pX) X = *pX;
    else   X = vect(zero);
    // 2.2 subset or all? : set first body and construct percentile finder
    unsigned Ns = S->N_subset();
    bool     all= Ns==S->N_bodies();
    B = all? S->begin_all_bodies() : S->begin_subset_bodies();
    FindPercentile<double> FP(Ns, all? SetAll : SetSub, N);
    if(B) falcON_ErrorN("Manipulator \"lagrange\": body number mismatch\n");
    double Mtot=FP.TotalWeight();
    // 2.3 loop mass fractions and get lagrange radii
    OUT << print(S->time(),12,8);
    for(int i=0; i!=N; ++i) {
      double M = F[i] * Mtot;
      double R = std::sqrt(FP.PositionOfCumulativeWeight(M));
      OUT << ' ' << print(R,8,2);
    }
    OUT << std::endl;
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} }

__DEF__MAN(falcON::Manipulate::lagrange)
