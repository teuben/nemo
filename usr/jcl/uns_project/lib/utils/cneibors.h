#ifndef CNEIBORS_H
#define CNEIBORS_H

#include <vector>
#include <ctree.h>
#include <cmath>

namespace jcltree {
//------------------------------------------------------------------------------
// CRadiusId
// class to store radius and id of the particles
class CDistanceId {
public:
  CDistanceId(double _r, int _id) {
    distance2 = _r;
    id        = _id;
  }
  int getId()           const { return id;             }
  double getDistance2() const { return distance2;      }
  double getDistance()  const { return sqrt(distance2);}

  static bool sortD(const CDistanceId& a, const CDistanceId& b) {
    return (a.distance2<b.distance2);
  }

private:
  int id;
  double distance2;
};

//------------------------------------------------------------------------------
// CNeibors
// class to find neighbors
class CNeibors {

public:
    CNeibors(const CTree * tree,const double _rneib0=0.00001);
    template <class T>  void process(const T * _pos0, const int _nneib, std::vector<CDistanceId> * _neib);
                        void process(const int i, const int _nneib, std::vector<CDistanceId> * _neib);
    void setRneib0(const double r) {
      rneib = r;
    }
    void setMaxRadius(const double _mr) {
      max_radius = _mr;
    }
    void setStopAtMaxRadius(const bool _b) {
      stop_max_radius = _b;
    }

    template <class T>  void direct(const T * _pos0, const int _nneib, std::vector<CDistanceId> * _neib);

private:
    double pos0[3];   // particle to find
    int nneib;        // #neibors to retreive
    std::vector<CDistanceId> * neib; //vector of neibors found

    CTree  * tree;   // ptr to the tree
    double rneib;    // searching radius
    double rnewneib; // new searching radius
    int total;       // total neibors found for a particle
    bool stop_max_radius; // stop if max radius reached
    double max_radius;

    // method
    void countPartInRadius();    // count particles in a given radius

    void searchTree(nodeptr p, double * cpos, double d);           // recursive tree walk descent

    bool openTreeNode(double * cpos, double d);         // open a tree node if necessary
};

}

#endif // CNEIBORS_H
