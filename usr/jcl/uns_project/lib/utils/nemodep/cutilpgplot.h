#ifndef CUTILPGPLOT_H
#define CUTILPGPLOT_H

#include <vector>
#include <string>
#include <cpgplot.h>
#include <iostream>
namespace jclut {

#define INIT_VECTOR(type, name, ...) \
static const type name##_a[] = __VA_ARGS__; \
std::vector<type> name##_v(name##_a, name##_a + sizeof(name##_a) / sizeof(*name##_a)); \
name = name##_v;


class CRainBow;
class CGray;
class CHeat;

class CPalet {

public:
  CPalet() {}
  void set() {
    cpgctab(&L[0], &R[0], &G[0], &B[0], L.size(), contrast, bright);
  }

protected:
  std::vector<float > L,R,G,B;
  float contrast, bright;

};

class CRainBow : public CPalet {
public:
  CRainBow() {
    INIT_VECTOR(float, L , {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7 } );
    //INIT_VECTOR(float, L , {0.0, 0.125, 0.25, 0.375, 0.50, 0.625, 0.75, 0.875, 1.0 } );
    INIT_VECTOR(float, R , { 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0 } );
    INIT_VECTOR(float, G , { 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0 } );
    INIT_VECTOR(float, B , { 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0 } );
    contrast = 1.0;
    bright   = 0.5;
  }
};

class CHeat : public CPalet {
public:
  CHeat() {
    INIT_VECTOR(float, L , { 0.0, 0.2, 0.4, 0.6, 1.0 } );
    INIT_VECTOR(float, R , { 0.0, 0.5, 1.0, 1.0, 1.0} );
    INIT_VECTOR(float, G , { 0.0, 0.0, 0.5, 1.0, 1.0} );
    INIT_VECTOR(float, B , { 0.0, 0.0, 0.0, 0.3, 1.0} );
    contrast = 1.0;
    bright   = 0.5;
  }
};

class CGray : public CPalet {
public:
  CGray() {
    INIT_VECTOR(float, L , { 0.0, 1.0 } );
    INIT_VECTOR(float, R , { 0.0, 1.0 } );
    INIT_VECTOR(float, G , { 0.0, 1.0 } );
    INIT_VECTOR(float, B , { 0.0, 1.0 } );
    contrast = 1.0;
    bright   = 0.5;
  }
};

class CIraf : public CPalet {
public:
  CIraf() {
    INIT_VECTOR(float, L , { 0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0 } );
    INIT_VECTOR(float, R , { 0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0 } );
    INIT_VECTOR(float, G , { 0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0 } );
    INIT_VECTOR(float, B , { 0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0 } );
    contrast = 1.0;
    bright   = 0.5;
  }
};

class CAips : public CPalet {
public:
  CAips() {
    INIT_VECTOR(float, L , { 0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
                             0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0 } );
    INIT_VECTOR(float, R , { 0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 } );
    INIT_VECTOR(float, G , { 0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
                             0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0 } );
    INIT_VECTOR(float, B , { 0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} );
    contrast = 1.0;
    bright   = 0.5;
  }
};

class CUtilPgplot {
public:
  CUtilPgplot();

  void selectCMap(const int cmap) {
    switch (cmap) {
    case 0: {
      CRainBow rb;
      rb.set();
      break;
    }
    case 1: {
      CHeat heat;
      heat.set();
      break;
    }
    case 2: {
      CGray gray;
      gray.set();
      break;
    }
    }
  }

};

}
#endif // CUTILPGPLOT_H
