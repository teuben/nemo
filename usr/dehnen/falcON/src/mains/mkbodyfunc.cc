// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    src/mains/mkbodyfunc.cc                                            
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2006                                                               
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2006  Walter Dehnen                                            
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
// v 0.0    28/07/2006  WD created                                              
////////////////////////////////////////////////////////////////////////////////
#define falcON_VERSION   "0.0"
#define falcON_VERSION_D "28-jul-2006 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "mkbodyfunc"
#endif
#define falcON_RepAction 0                         // no action reporting       
//-----------------------------------------------------------------------------+
#include <iostream>                                // C++ I/O                   
#include <fstream>                                 // C++ file I/O              
#include <iomanip>                                 // C++ I/O formatting        
#include <public/bodyfunc.h>                       // body functions            
#include <main.h>                                  // main & NEMO stuff         
//------------------------------------------------------------------------------
string defv[] = {
  "expr=\n            bodyfunc expression                                ",
  "show=f\n           show existing bodyfunc expressions?                ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string
usage = "create a bodyfunc expression";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  bodyfunc *bF = hasvalue("expr")? new bodyfunc(getparam("expr")) : 0;
  if(getbparam("show"))
    bodyfunc::print_db(std::cout);
  if(bF) falcON_DEL_O(bF);
}
