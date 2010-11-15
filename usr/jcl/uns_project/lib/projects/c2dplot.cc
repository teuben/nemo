// ============================================================================
// Copyright Jean-Charles LAMBERT - 2009-2010                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
#include <cpgplot.h>
#include <iostream>
#include <boost/timer.hpp>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <cstring>
#include "c2dplot.h"
#include "csnaptools.h"
#include "ctimer.h"

using namespace uns_proj;
using namespace jclut;
using namespace std;
// ----------------------------------------------------------------------------
// contructor
template <class T> C2dplot<T>::C2dplot(const int _nthreads, const int _pixel, const int _dimx, const int _dimy, const T _g)
{
  nthreads = _nthreads;
  dimx     = _dimx;
  dimy     = _dimy;
  pixel    = _pixel;
  g        = _g;

  // create gaussian
  gaussian = new CGaussian<float>(pixel,g);

  // allocate memory for threads
  for (int i=0; i<nthreads; i++) {
    tab[i] = new float[dimx*dimy];
  }
}
template C2dplot<float >::C2dplot(const int _nthreads, const int pixel, const int _dimx, const int _dimy, const float _g);
template C2dplot<double>::C2dplot(const int _nthreads, const int pixel, const int _dimx, const int _dimy, const double _g);
// ----------------------------------------------------------------------------
// compute
template <class T> void C2dplot<T>::compute(std::string _dev, const int _no_frame,const int _nbody, T * _pos , float _range[3][2], 
                                            std::string _title,std::string _sel_comp,std::string _filename, float _timu,
                                            bool _xy, bool _xz, bool _zy, bool _sview)
{
  // get paramaters
  dev       = _dev;
  no_frame  = _no_frame;
  nbody     = _nbody;
  pos       = _pos;
  title     = _title;
  sel_comp  = _sel_comp;
  filename  = _filename;
  time      = _timu;
  xy        = _xy;
  xz        = _xz;
  zy        = _zy;
  sview     = _sview;
  
  memcpy(range,_range,sizeof(float)*6);
  
  //
  std::string outdev=dev;
  int showtext=0;
  int nbview=1;
  // count how many views
  if (xy) nbview++;
  if (xz) nbview++;
  if (zy) nbview++;
  nbview--;
  if (nbview==1) sview=true;
  int nbxview=1;
  if (sview) nbxview=nbview;
  
  if (sview) { // one single view
    outdev = buildFrameName("",no_frame);
    cpgopen(outdev.c_str());
    cpgsubp(nbxview,1);
  }  
  
  // Process XY view
  drawImage(xy,0,1,nbview,showtext);
  // Process XZ view
  drawImage(xz,0,2,nbview,showtext);
  // Process ZY view
  drawImage(zy,2,1,nbview,showtext);
  
  if (sview) {
    cpgask(1);
    cpgend(); 
  }
}
template void C2dplot<float >::compute(std::string pic,const int no_frame, const int nbody, float  * pos , float range[3][2], 
                                       std::string title,std::string sel_comp, std::string filename, 
                                       float timu,bool xy, bool xz, bool yz, bool sview);
template void C2dplot<double>::compute(std::string pic,const int no_frame, const int nbody, double * pos , float range[3][2], 
                                       std::string title,std::string sel_comp, std::string filename, 
                                       float timu,bool xy, bool xz, bool yz, bool sview);
// ----------------------------------------------------------------------------
// drawImage
template <class T> void C2dplot<T>::drawImage(const bool disp,const int xaxis, const int yaxis, const int nbview, int &showtext)
{
  std::string outdev = dev;
  std::string label[3] = { "X","Y","Z"};
  if (disp) {
    if (!sview) {
      // build extension name
      std::string ext="_"+label[xaxis]+label[yaxis]; // example "_XY"
      outdev = buildFrameName(ext,no_frame);
      cpgopen(outdev.c_str());
      cpgsubp(1,1);
    }
    // compute image in parallel
    computeImage(xaxis,yaxis);
    cpgsci(1);
    cpglab(label[xaxis].c_str(),label[yaxis].c_str(),""); // example "X" "Y"
    if (!sview||nbview==1||(sview&&!showtext)) {
      displayText((nbview==1?false:sview));
      showtext++;
    }
    if (!sview) {
      cpgclos();
    }
  }
}
// ----------------------------------------------------------------------------
// computeImage
template <class T> void C2dplot<T>::computeImage(const int xaxis, const int yaxis)
{
  float zmin,zmax;  
  
  float bright=0.5;
  float contrast=1.0;

  float RL[9] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
  float RR[9] = { 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
  float RG[9] = { 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
  float RB[9] = { 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};

  // look for indexes in the range axis
  findIndexes(xaxis,yaxis);
  std::cerr << "XYZ particles => Nb indexes found = " << indexes.size() << "\n";

  xmin=ymin=min(range[xaxis][0],range[yaxis][0]);
  xmax=ymax=max(range[xaxis][1],range[yaxis][1]);

  CTimer timing;
  // launch computation on threads
  startWorkers(nbody,pos,xaxis,yaxis,zmin,zmax);
  std::cerr << "Work done cpu time : "<< timing.cpu() << "\n";
  std::cerr << "Work done elapsed  : "<< timing.elapsed() << "\n";
  
  // transfert function
  float tr[6] = { xmin, (xmax-xmin)/(float)dimx, 
		  0, ymin, 
		  0, (ymax-ymin)/(float)dimy};
  
  // compute the boundary min/max
  xmin=range[xaxis][0];
  xmax=range[xaxis][1];
  ymin=range[yaxis][0];
  ymax=range[yaxis][1];
  // viewpart according to the boundary
  cpgenv(xmin,xmax,ymin,ymax,1,0);
  // create color image gradient
  cpgctab(RL, RR, RG, RB, 9, contrast, bright);
  cpgimag(tab[0],dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);  
}
// ----------------------------------------------------------------------------
//  buildFrameName
template <class T> std::string C2dplot<T>::buildFrameName(std::string label,const int idx)
{
  std::string outdev=dev;
  if (outdev != "?") {
    outdev="";
    stringstream ss;
    ss << dev << label<< "." << setw(5) << setfill('0') << idx << ".gif/gif";
    outdev = ss.str();
  }
  return outdev;
}
// ----------------------------------------------------------------------------
// displayText()
template <class T> void C2dplot<T>::displayText(bool sview)
{
  float factor=1.0;
  if (!sview) factor=2.0;
  // title
  cpgstbg(15);
  cpgsch(2.5/factor);
  cpgmtxt("t",2.,0.,0., title.c_str());
  cpgstbg(0);
  // filename (on the left)
  cpgsci(1);
  cpgsch(1.5/factor);
  std::string basename=CSnaptools::basename(filename);
  cpgmtxt("t",1.8,0.,0.,basename.c_str() );
  
  // component
  cpgsci(1);
  cpgsch(2.0/factor);
  cpgmtxt("t",2.5,1.,1., sel_comp.c_str());
  // Time (on the right)
  cpgsci(1);
  cpgsch(1.5/factor);
  std::stringstream ss;
  ss << "time: " << std::setw(7) << std::fixed << std::setprecision(3) << time;
  cpgmtxt("t",0.5,1.,1., (ss.str()).c_str());
  // nbody (on the left)
  cpgsci(1);
  cpgsch(1.5/factor);
  ss.str(""); // clear sstream;
  ss  <<"nbody: " << std::setw(9) << nbody;
  cpgmtxt("t",0.5,0.,0., (ss.str()).c_str());
  
  cpgsch(1.0); // set back the font size
  cpgsci(1);   // set back to default color
}
// ----------------------------------------------------------------------------
// findIndexes()
template <class T> void C2dplot<T>::findIndexes(const int xaxis, const int yaxis)
{
  indexes.clear();
  T * pp=pos;
  int r=xaxis;
  int f=yaxis;
  /* loop on all particles */
  for (int i=0 ; i<nbody; i++, pp+=3) {
    if (
	pp[r] >= range[r][0] &&
	pp[r] <= range[r][1] &&
	pp[f] >= range[f][0] &&
	pp[f] <= range[f][1] ) {
      indexes.push_back(i);
    }
  }
}
// ----------------------------------------------------------------------------
// startWorkers
// start all the thread in parallel
//
template <class T> void C2dplot<T>::startWorkers(const int nbody, T * data, const int xaxis, const int yaxis, float& zmin, float& zmax)
{
  int npart=indexes.size()/nthreads;
  int offset=0;
  // build vector of threads and start them
  for (int x=0; x<nthreads; x++) {
    if (x==nthreads-1) { // last threads
      npart=indexes.size()-offset;
    }
    //std::cerr << "startWorkers: "<<x<<" npart="<<npart<<"\n";
    boost::shared_ptr<boost::thread>  vv(new boost::thread(boost::bind(&C2dplot::worker, this, x,offset,npart,data,xaxis,yaxis)));
    v_threads.push_back(vv) ;

    offset+=npart;
  }
  // wait the end of each threads
  for (std::vector<boost::shared_ptr<boost::thread> >::iterator 
	 it=v_threads.begin();
       it != v_threads.end(); it++) {
    (*it)->join();
  }
  // accumulate all frames computed in parallel on the first frame
  for (int x=1; x<nthreads; x++) {
    for (int i=0;i<dimy; i++)
      for (int j=0;j<dimx; j++) {
	tab[0][i*dimx+j]+=tab[x][i*dimx+j] ;
      }
  }
  // find zmin zmax
  zmax=-1000000.0;
  zmin=10000000.0;
  for (int i=0;i<dimy; i++)
    for (int j=0;j<dimx; j++) {
      zmax=std::max(zmax, tab[0][i*dimx+j]);
      float zzmin=std::min(zmin, tab[0][i*dimx+j]);
      if (zzmin != 0.0)  zmin=zzmin;

    }

  for (int i=0;i<dimy; i++)
    for (int j=0;j<dimx; j++) {
      if (tab[0][i*dimx+j]==0.0) tab[0][i*dimx+j]=zmin ;

    }
  
  std::cerr << "FIRST zmax="<<zmax<<" zmin="<<zmin<<"\n";

  // compute LOG
  zmax=log(zmax);
  zmin=log(zmin);
  for (int i=0;i<dimy; i++)
    for (int j=0;j<dimx; j++) {
      if (tab[0][i*dimx+j] != 0.0) {
	tab[0][i*dimx+j] = log(tab[0][i*dimx+j]);///zmax; // normalize
      }
    }
    
}
// ----------------------------------------------------------------------------
// worker
// function executed in parallel by each threads
template <class T> void C2dplot<T>::worker(const int ithread, const int offset, const int npart, T * data,const int xaxis, const int yaxis)
{
  // Reset each local array
  for (int i=0;i<dimy; i++)
    for (int j=0;j<dimx; j++)
      tab[ithread][i*dimx+j] = 0.0;

  //displayIndexes();

  // loop on all sub particles
  for (int i=0; i<npart; i++) {
    int    ii = indexes[(offset+i)];
    float  xx = data[ii*3+xaxis];
    int     x = ((xx-xmin)/(xmax-xmin))*(dimx-1);
    float  yy = data[ii*3+yaxis];
    int     y = ((yy-ymin)/(ymax-ymin))*(dimy-1);

    assert(x<dimx);
    assert(y<dimy);
    // apply gaussian on XY coordinates
    gaussian->applyOnArrayXY(tab[ithread],dimx,dimy,x,y);
  }
}
//
// ----------------------------------------------------------------------------
