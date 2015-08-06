// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/mains/snapprop.cc                                               
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2004-2008                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2004-2008 Walter Dehnen                                        
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
// v 0.0    16/07/2004 WD created snapsum                                       
// v 0.1    19/07/2004 WD new bodyfunc                                          
// v 0.2    20/07/2004 WD first tests                                           
// v 1.0    21/07/2004 WD first running version, seems to work as desired       
// v 1.1    20/08/2004 WD much improved & extended version of bodyfunc.h        
// v 1.2    25/08/2004 WD deBUGged error in bodyfunc (Segfault with Min/Max{})  
// v 1.3    08/11/2004 WD parameters to bodiesfunc expressions                  
// v 2.0    24/06/2004 WD new falcON                                            
// v 2.1    13/06/2005 WD changes in fieldset                                   
// v 2.2    04/07/2006 WD made public (along with bodyfunc)                     
// v 2.2.1  10/09/2008 WD happy gcc 4.3.1                                      
////////////////////////////////////////////////////////////////////////////////
#define falcON_VERSION   "2.2.1"
#define falcON_VERSION_D "10-sep-2008 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "snapprop"
#endif
#define falcON_RepAction 0                         // no action reporting       
//-----------------------------------------------------------------------------+
#include <public/bodyfunc.h>                       // body functions            
#include <main.h>                                  // main & NEMO stuff         
using namespace falcON;
//------------------------------------------------------------------------------
const char*defv[] = {
  "in=???\n           input file                                         ",
  "prop=???\n         bodiesfunc expression (see man page) to evaluate   ",
  "pars=\n            parameters, must match requirements from prop      ",
  "times=all\n        times to process                                   ",
  "givetime=t\n       print time as well as property                     ",
  "prec=6\n           significant digits for floating-point numbers      ",
  "zeromissing=f\n    zero missing body properties (or error out)?       ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
const char*usage =
    "snapprop -- evaluates bodies function over snapshot, reports to stdout";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  // set up parameters                                                          
  const bool    Z   (getbparam("zeromissing"));
  const nemo_in IN  (getparam("in"));
  bodiesfunc    BF  (getparam("prop"));
  const int     prec(getiparam("prec"));
  const bool    time(getbparam("givetime"));
  real         _P[10], *P=0;
  snapshot      SHOT;
  if(BF.npar()) {
    int n = getaparam_z("pars",_P,10);
    if(n == 0)
      falcON_THROW("prop=\"%s\" requires %d parameters, none given",
		   getparam("prop"),BF.npar());
    if(n < BF.npar())
      falcON_THROW("prop=\"%s\" requires %d parameters, only %d given",
		   getparam("prop"),BF.npar(),n);
    P = _P;
  }
  // loop snapshots, process them and make outputs                              
  while(IN.has_snapshot()) {
    fieldset got;
    if(!SHOT.read_nemo(IN,got,BF.need(),getparam("times"),0)) continue;
    if(!got.contain(BF.need())) {
      fieldset miss = got.missing(BF.need());
      if(Z) {
	falcON_Warning("data '%s' missing; will zero them",word(miss));
	SHOT.reset_data(miss);
      } else
	falcON_THROW("data on '%s' missing",word(miss));
    }
    if(time)
      std::cout << std::setprecision(prec) << SHOT.time() << ": ";
    switch(BF.type()) {
    case 'b':
      std::cout << BF.func<bool>(SHOT,P);
      break;
    case 'i':
      std::cout << BF.func<int>(SHOT,P);
      break;
    case 'r':
      std::cout << std::setprecision(prec) << BF.func<real>(SHOT,P);
      break;
    case 'v': {
      vect x = BF.func<vect>(SHOT,P);
      std::cout << std::setprecision(prec) << x[0] << ','
		<< std::setprecision(prec) << x[1] << ','
		<< std::setprecision(prec) << x[2];
    } break;
    default: falcON_THROW ("unknown type");
    }
    std::cout<<std::endl;
  }
}
//------------------------------------------------------------------------------
