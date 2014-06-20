// ============================================================================
// Copyright Jean-Charles LAMBERT - 2009-2013
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
#include <cstdlib>
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
      
#if 0
  T pi=atan(1.0)*4.;
  T sqrtpi=sqrt(2*pi);
  T ig=1./g;
  T isqrtpi=1./sqrtpi;
  T halfp=(T )pixel/2.;
  for (int i=0; i<pixel; i++) {
    for (int j=0; j<pixel; j++) {
      T distance=sqrt((i-halfp)*(i-halfp)+(j-halfp)*(j-halfp));
      //T gauss=exp(-(distance)*(distance)/g/g/2.0)/sqrtpi/g;      
      T gauss=exp(-(distance)*(distance)*ig*ig*0.5)*isqrtpi*ig;
      gaussian[i*pixel+j]=gauss;
    }
  }
#else
  createGaussianMap(pixel);
#endif
}
template CGaussian<float >::CGaussian(const int _pixel, const float  _g);
template CGaussian<double>::CGaussian(const int _pixel, const double _g);

// ============================================================================
//
template <class T> void CGaussian<T>::createGaussianMap(const int pixel)
{
  int N=pixel;
  T *M = new T[2*N*N];
  //T *B = new T[N*N];
  T X,Y,Y2,Dist;
  T Incr = 2.0f/N;
  int i=0;
  int j = 0;
  Y = -1.0f;
    //float mmax = 0;
  for (int y=0; y<N; y++, Y+=Incr)
  {
    Y2=Y*Y;
    X = -1.0f;
    for (int x=0; x<N; x++, X+=Incr, i+=2, j++)
    {
      Dist = (float)sqrtf(X*X+Y2);
      if (Dist>1) Dist=1;
      //M[i+1] = M[i] = evalHermite(1.0f,0.1,0.1,0.1,Dist);
      M[i+1] = M[i] = evalHermite(1.0f,0,0,0,Dist);
      gaussian[j] =  M[i];
      //std::cerr << "j="<<j<< " g="<<gaussian[j] << "\n";
    }
  }
  delete [] M;
  //return(B);
}
// ============================================================================
//
template <class T> T CGaussian<T>::evalHermite(const T pA, const T pB, const T vA, const T vB, const T u)
{
  T u2=(u*u), u3=u2*u;
  T B0 = 2*u3 - 3*u2 + 1;
  T B1 = -2*u3 + 3*u2;
  T B2 = u3 - 2*u2 + u;
  T B3 = u3 - u;
  return( B0*pA + B1*pB + B2*vA + B3*vB );
}
// ----------------------------------------------------------------------------
// applyOnArrayXY
template <class T> void CGaussian<T>::applyOnArrayXY(T * tab, const int dimx,const int dimy,
                                                     const int x, const int y,
                                                     const T weight, const int psort)
{
  int halfp=pixel/2.;
  int nindex=0;
  for (int i=0; i<pixel; i++) {
    for (int j=0; j<pixel; j++) {
      if ((x-halfp+j)>=0 && (x-halfp+j)<dimx && (y-halfp+i)>=0 && (y-halfp+i)<dimy) {
        int index=(int) ((y-halfp+i)*dimx+(x-halfp+j));
        if ( index <0 or index > (dimx*dimx)) {
          std::cerr << "error index = " << index << "\n";
          nindex++;
        } else {

          switch (psort) {
          case 0:
            tab[index]+=(gaussian[i*pixel+j]*weight); // accumulated
            break;
          case 1:
            tab[index]=std::max(tab[index],gaussian[i*pixel+j]*weight); // max front
            break;
          default:
            std::cerr << "bad psort value ["<<psort<<"], file:"<< __FILE__<< " at line:" << __LINE__ << "\n";
            std::exit(1);
            break;
          }
          //tab[index]+=(gaussian[i*pixel+j]*weight); // accumulated
          //std::cerr << x <<" " <<y <<" "  << " " << gaussian[i*pixel+j] << "\n";
          //tab[index]=std::max(tab[index],gaussian[i*pixel+j]*weight); // max front
        }
      }
    }
  }
#if 0
  nindex=0;
  for (int i=0; i<pixel; i++) {
    for (int j=0; j<pixel; j++) {
      if ((x-halfp+j)>=0 && (x-halfp+j)<dimx && (y-halfp+i)>=0 && (y-halfp+i)<dimy) {
          int index=(int) ((y-halfp+i)*dimx+(x-halfp+j));
          if ( index <0 or index > (dimx*dimx)) {
            std::cerr << "error index = " << index << "\n";
            nindex++;
          } else {
            //tab[index]+=(gaussian[i*pixel+j]*weight);
            tab[index]+=weight;
          }
      }
    }
  }
#endif
}

// ----------------------------------------------------------------------------
// applyOnArrayXY
template <class T> void CGaussian<T>::computeOnArrayXY(T * tab, const int dimx,const int dimy,
                                                     const int x, const int y,
                                                     const T weight, const int pixel)
{
  delete [] gaussian;
  gaussian =new T[pixel*pixel];
  createGaussianMap(pixel);
  int halfp=pixel/2.;
  int nindex=0;
  for (int i=0; i<pixel; i++) {
    for (int j=0; j<pixel; j++) {
      if ((x-halfp+j)>=0 && (x-halfp+j)<dimx && (y-halfp+i)>=0 && (y-halfp+i)<dimy) {
        int index=(int) ((y-halfp+i)*dimx+(x-halfp+j));
        if ( index <0 or index > (dimx*dimx)) {
          std::cerr << "error index = " << index << "\n";
          nindex++;
        } else {
          tab[index]+=(gaussian[i*pixel+j]*weight);
          //std::cerr << x <<" " <<y <<" "  << " " << gaussian[i*pixel+j] << "\n";
        }
      }
    }
  }
}
// ----------------------------------------------------------------------------
// applyOnArrayXYStep
template <class T> void CGaussian<T>::applyOnArrayXYStep(T * tab, const int dimx,const int dimy,
                                                     const int x, const int y,
                                                     const T weight)
{
  int halfp=pixel/2.;
  int nindex=0;
  for (int i=0; i<pixel; i++) {
    for (int j=0; j<pixel; j++) {
      if ((x-halfp+j)>=0 && (x-halfp+j)<dimx && (y-halfp+i)>=0 && (y-halfp+i)<dimy) {
        int index=(int) ((y-halfp+i)*dimx+(x-halfp+j));
        if ( index <0 or index > (dimx*dimx)) {
          std::cerr << "error index = " << index << "\n";
          nindex++;
        } else {
          //tab[index]+=(gaussian[i*pixel+j]*weight);
          tab[index]=std::max(tab[index],gaussian[i*pixel+j]*weight);
        }
      }
    }
  }
}
template void CGaussian<float >::applyOnArrayXY(float  * tab, const int dimx,const int dimy, const int x, const int y, const float, const int);
template void CGaussian<double>::applyOnArrayXY(double * tab, const int dimx,const int dimy, const int x, const int y, const double, const int);

template void CGaussian<float >::computeOnArrayXY(float  * tab, const int dimx,const int dimy, const int x, const int y, const float, const int);
template void CGaussian<double>::computeOnArrayXY(double * tab, const int dimx,const int dimy, const int x, const int y, const double, const int);



//
