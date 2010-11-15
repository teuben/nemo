// ============================================================================
// Copyright Jean-Charles LAMBERT - 2009 - 2010                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

#include <iostream>
#include <cmath>
#include "cgaussian.h"

using namespace jclut;

// ----------------------------------------------------------------------------
// constructor
template <class T> CGaussian<T>::CGaussian(const int _pixel, const T _g)
{
  pixel=_pixel;
  g=_g;
  // allocate gaussian array
  gaussian =new T[pixel*pixel];
      
  T pi=atan(1.0)*4.;
  T sqrtpi=sqrt(2*pi);
  T halfp=(T )pixel/2.;
  for (int i=0; i<pixel; i++) {
    for (int j=0; j<pixel; j++) {
      T distance=sqrt((i-halfp)*(i-halfp)+(j-halfp)*(j-halfp));
      T gauss=exp(-(distance)*(distance)/g/g/2.0)/sqrtpi/g;      
      gaussian[i*pixel+j]=gauss;
    }
  }   
}
template CGaussian<float >::CGaussian(const int _pixel, const float  _g);
template CGaussian<double>::CGaussian(const int _pixel, const double _g);

// ----------------------------------------------------------------------------
// applyOnArrayXY
template <class T> void CGaussian<T>::applyOnArrayXY(T * tab, const int dimx,const int dimy, const int x, const int y)
{
  T halfp=(float )pixel/2.;
  int nindex=0;
  for (int i=0; i<pixel; i++) {
    for (int j=0; j<pixel; j++) {
      if ((x-halfp+j)>=0 && (x-halfp+j)<dimx && (y-halfp+i)>=0 && (y-halfp+i)<dimy) {
	int index=(int) ((y-halfp+i)*dimx+(x-halfp+j));
	if ( index <0 or index > (dimx*dimx)) {
	  std::cerr << "error index = " << index << "\n";
	  nindex++;
	} else {
	  tab[index]+=gaussian[i*pixel+j];
	}
      }
    }
  }  
}

template void CGaussian<float >::applyOnArrayXY(float  * tab, const int dimx,const int dimy, const int x, const int y);
template void CGaussian<double>::applyOnArrayXY(double * tab, const int dimx,const int dimy, const int x, const int y);



//
