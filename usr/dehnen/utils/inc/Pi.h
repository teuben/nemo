// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    inc/Pi.h                                                           
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    1994-2005                                                          
///                                                                             
/// \brief   contains the definition of some constants around Pi                
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1994-2005  Walter Dehnen                                       
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
#ifndef WDutils_included_Pi_h
#define WDutils_included_Pi_h
////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  const double Pi   = 3.14159265358979323846264338328;  // Pi                   
  const double Pih  = 0.5  * Pi;                        // Pi/2                 
  const double Piq  = 0.25 * Pi;                        // Pi/4                 
  const double Pi3h = 3.   * Pih;                       // 3*Pi/2               
  const double TPi  = 2.   * Pi;                        // 2*Pi                 
  const double FPi  = 4.   * Pi;                        // 4*Pi                 
  const double iFPi = 0.25 / Pi;                        // 1/(4*Pi)             
  const double FPit = 4.   * Pi/3.;                     // 4*Pi/3               
  const double iFPit= 0.75 / Pi;                        // 3/(4*Pi)             
  const double SPi  = 1.772453850905516027298167483341;	// Sqrt[Pi]             
  const double STPi = 2.506628274631000502415765284811; // Sqrt[2 Pi]           
#if  defined(WDutils_included_complex) || defined(_CPP_COMPLEX)
  const std::complex<double> IMAG(0,1);	                // i                    
  const std::complex<double> ITPi(0,TPi);               // 2 i Pi               
#endif
}
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_Pi_h
