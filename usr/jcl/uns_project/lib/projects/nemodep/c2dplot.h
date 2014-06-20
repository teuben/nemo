// ============================================================================
// Copyright Jean-Charles LAMBERT - 2009-2011                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

/*
  @author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
*/

#ifndef C2DPLOT_H
#define C2DPLOT_H

#include <string>
#include <algorithm>    // std::sort
#include <vector>
#include <iostream>
#ifndef NOBOOST
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#endif
#include "cgaussian.h"

#define NTHREAD_MAX 256

#ifndef NO_CUDA
#define NB_BLOCK 2

#endif
namespace uns_proj {
  typedef struct {
     int devid;
     int nblock, offset;
   } t_cuda_grid;
  
//------------------------------------------------------------------------------
// CPartProp, class to store image properties foreach particles displayed
  class CPartProp {
  public:
    CPartProp(int _x, int _y, float _prop, float _hsml) {
      x    = _x;
      y    = _y;
      prop = _prop;
      hsml = _hsml;
    }
    static bool mySort(const CPartProp&a, const CPartProp &b) {
      return a.prop < b.prop; // ascending tree
    }
    static bool mySortR2(const CPartProp&a, const CPartProp &b) {
      return ((a.x*a.x + a.y*a.y) < (b.x*b.x+b.y*b.y)); // ascending tree
    }
    int x,y; // x,y image coordinates
    float prop, hsml; // properties, hsml
  };

  template <class T> class C2dplot {
  public:
    C2dplot(const int,const int,const int, const int, const T);
    
    void compute(std::string pic, const int _no_frame,const int _nbody, T * _pos , 
		 float _range[3][2], std::string _title, 
         std::string _sel_comp, std::string _filename, const float _time,
                 bool _xy, bool _xz, bool _zy, bool _sview, T * _weight, const int psort, T * hsml,
    const int _itf, const bool wedge, std::string legend,const int _cmap);


  private:
    int psort; // control properties sorting
    int nthreads, dimx, dimy, pixel;
    int itf; // image transfer function
    bool wedge; // toogle on/off color bar display
    std::string legend;
    int cmap; // color map
    jclut::CGaussian<float> * gaussian;
    std::vector <CPartProp> pvec; // vector of particles properties
    T g;
    float * tab[NTHREAD_MAX];
    std::vector <int> indexes;
    // drawImage
    void drawImage(const bool disp,const int xaxis, const int yaxis, const int nbview, int &showtext);
    // computeImage
    void computeImage(const int xaxis, const int yaxis);
    // buildFrameName
    std::string buildFrameName(std::string label,const int idx);
    // display text
    void displayText(bool sview);
    void findIndexes(const int xaxis, const int yaxis);
    void displayIndexes() {
      for (unsigned int i=0; i<indexes.size(); i++) {
	std::cerr << indexes[i] << "\n";
      }
    }
    // Parallel tasks
    void startWorkers(const int nbody, T * data, const int xaxis, const int yaxis,float& zmin,float& zmax);
    void worker(const int ithread, const int offset, const int nbody, T * data, const int xaxis, const int yaxis);
#ifndef NOBOOST
    std::vector<boost::shared_ptr<boost::thread> > v_threads; 
    boost::mutex io_mutex;
#endif
    float xmin,xmax,ymin,ymax;
    
    std::string dev, title, sel_comp, filename;
    int no_frame, nbody;
    float time,range[3][2];
    T * pos;
    bool xy, xz, zy, sview;
    T * weight, * hsml;
#ifndef NO_CUDA
    // CUDA stuffs
    int GPU_N;
    t_cuda_grid * device_data;
    float 
        * cuda_pos;      // particles's positions
    int 
        c_nbody;         // #cuda bodies particles
    void initCuda();
    void startCuda(const int nbody, T * data, const int xaxis, const int yaxis,float& zmin,float& zmax);
    void workerCudaThread(const int id, T * data, const int xaxis, const int yaxis);
#endif
  };
}
#endif
//
