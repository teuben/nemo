#ifndef CTREE_H
#define CTREE_H

#include <vector>

// ----------------------------------------------------------------------------
// NOTES
// ----------------------------------------------------------------------------
// this work is an adaptation of the tree used by Jun's Makino hackdens method
// from NEMO project. This version is C++ ized and works also on data set with
// particles having same positions. Allocator has been modified to allocate more
// memory space if required.
// ----------------------------------------------------------------------------
namespace jcltree {

/*
 * BODY and CELL data structures are used to represent the tree:
 *
 *         +-----------------------------------------------+
 * root--> | CELL: mass, pos, quad, /, o, /, /, /, /, o, / |
 *         +---------------------------|--------------|----+
 *                                     |              |
 *    +--------------------------------+              |
 *    |                                               |
 *    |    +--------------------------------+         |
 *    +--> | BODY: mass, pos, vel, acc, phi |         |
 *         +--------------------------------+         |
 *                                                    |
 *    +-----------------------------------------------+
 *    |
 *    |    +-----------------------------------------------+
 *    +--> | CELL: mass, pos, quad, o, /, /, o, /, /, o, / |
 *         +------------------------|--------|--------|----+
 *                                 etc      etc      etc
 */


 // NODE: data common to BODY and CELL structures.


typedef struct {
    short type;                   // code for node type
    double mass;                  // total mass of node
    double pos[3];	          // position of node
} node, *nodeptr;

#define Type(x) (((nodeptr) (x))->type)
#define Mass(x) (((nodeptr) (x))->mass)
#define Pos(x)  (((nodeptr) (x))->pos)


 // BODY: data structure used to represent particles.

typedef int int_hack; // int_hack
#define BODY 01                 // type code for bodies

typedef struct {
    short type;
    double mass;                  // mass of body
    double pos[3];                // position of body
    int id;                       // body's id
    int_hack level;               // tree's level
  /* double d2; */                // squared distance from the p0 particle
} body, *bodyptr;

#define Body    body
#define Vel(x)   (((bodyptr) (x))->vel)
#define Dist2(x) (((bodyptr) (x))->d2)
#define Rho(x)   (((bodyptr) (x))->rho)
#define Hsml(x)  (((bodyptr) (x))->hsml)
#define Id(x)    (((bodyptr) (x))->id)
#define Level(x) (((bodyptr) (x))->level)


 // CELL: structure used to represent internal nodes of tree.

#define CELL 02                 // type code for cells

#define NSUB (1 << 3)           // subcells per cell

typedef struct {
    short type;
    double mass;                // total mass of cell
    double pos[3];              // cm. position of cell
    nodeptr subp[NSUB];         // descendents of cell
} cell, *cellptr;

#define Subp(x) (((cellptr) (x))->subp)

// structure to save ids of particles
// having same position
class CSamePos {
public:
  CSamePos(const int _i1, const int _i2) {
    id1 = _i1; id2 = _i2;
  }
  int getId1() { return id1;}
  int getId2() { return id2;}
private:
  int id1,id2;
};


#define IMAX (1 << (8 * sizeof(int) - 2))       /* highest bit */
#define NDIM 3
//------------------------------------------------------------------------------
// class CTree, store a set of 3D particles into a tree structure
//------------------------------------------------------------------------------
class CTree
{
public:
    template <class T> CTree(const int _nbody, const T * pos, const T * mass,
          const double _fcells=0.9, const double _rsize=4.0);
    ~CTree();
  nodeptr    getRoot()       const { return troot;       }  // tree node root
  double *   getRmin()             { return rmin;        }  // tree rmin (lower left coordinates)
  double     getRsize()      const { return rsize;       }  // tree size
  int        getLevelMax()   const { return level_max;   }  // tree depth max
  bodyptr    getBodyData()   const { return btab;     }  // tree stored data
  int        getNbody()      const { return nbody;    }  // #bodies stored in bodytab
  int        getTotalCells() const { return totalcell;}  // total cells allocated
  std::vector<CSamePos>
             * getSamePos()        { return &samepos; }  // return vector of particles with same position
private:
    // data
    int nbody;          // #bodies
    nodeptr troot;      // ROOT: origin of tree; declared as nodeptr for tree with only 1 body
    double  fcells;	// ratio of cells/bodies allocated
    double rmin[3];     // lower-left corner of coord. box
    double rsize;       // side-length of int. coord. box
    bodyptr btab;       // body tree structure
    int level_max;      // tree depth max

    // cells
    int ncell,maxcell;  // count cells in use, max available
    int totalcell;
    std::vector<cellptr> ctab;
    std::vector<CSamePos> samepos; // store particles with same position
    // functions
    void makeTree();

    void expandBox(bodyptr p);
    void loadTree(bodyptr p);
    bool intCoord(int_hack xp[3], double rp[3]);
    int_hack subIndex(int_hack x[3], int_hack l);
    void hackCofm( nodeptr q, int_hack l);
    cellptr makeCell();
};
}

#endif // CTREE_H
