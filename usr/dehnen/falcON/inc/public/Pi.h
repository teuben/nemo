// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// Pi.h                                                                        |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 1994-2003                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_Pi_h
#define falcON_included_Pi_h
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  const double Pi   = 3.14159265358979323846264338328;  // Pi                   
  const double Pih  = 0.5  * Pi;                        // Pi/2                 
  const double Piq  = 0.25 * Pi;                        // Pi/4                 
  const double Pi3h = 3.   * Pih;                       // 3*Pi/2               
  const double TPi  = 2.   * Pi;                        // 2*Pi                 
  const double FPi  = 4.   * Pi;                        // 4*Pi                 
  const double FPit = 4.   * Pi/3.;                     // 4*Pi/3               
  const double SPi  = 1.772453850905516027298167483341;	// Sqrt[Pi]             
  const double STPi = 2.506628274631000502415765284811; // Sqrt[2 Pi]           
#if defined(__COMPLEX__) || defined(_CPP_COMPLEX) || defined(__STD_COMPLEX)
  const std::complex<double> IMAG = std::complex<double>(0,1);	// i            
  const std::complex<double> ITPi = std::complex<double>(0,TPi);// 2 i Pi       
#endif
}
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_Pi_h
