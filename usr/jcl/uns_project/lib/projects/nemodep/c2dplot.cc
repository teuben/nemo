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
#include <cpgplot.h>
#include <iostream>
#include <ctime>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <assert.h>
#include <cmath>
//#include <algorithm>    // std::sort
#include "c2dplot.h"
#include "csnaptools.h"
#include "ctimer.h"
#include <limits>
#include <cutilpgplot.h>

using namespace uns_proj;
using namespace jclut;
using namespace std;

#define WEIGHT(X) ((weight==NULL)?1.0:weight[(X)])
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

#ifndef NO_CUDA
  initCuda();
  nthreads = GPU_N;
#endif
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
                                            bool _xy, bool _xz, bool _zy, bool _sview, T * _weight,
                                            const int _psort, T * _hsml, const int _itf, const bool _wedge, std::string _legend, const int _cmap)
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
  weight    = _weight;
  hsml      = _hsml;
  psort     = _psort;
  itf       = _itf;
  wedge     = _wedge;
  legend    = _legend;
  cmap      = _cmap;

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
    if (outdev=="?" && filename=="-") {
      outdev="/xs";
    }
    cpgopen(outdev.c_str());
    cpgsubp(nbxview,1);
  }  
  
  // Process XY view
  if (xy) drawImage(xy,0,1,nbview,showtext);
  // Process XZ view
  if (xz) drawImage(xz,0,2,nbview,showtext);
  // Process ZY view
  if (zy) drawImage(zy,2,1,nbview,showtext);
  
  if (sview) {
    cpgask(1);
    cpgend(); 
  }
}
template void C2dplot<float >::compute(std::string pic,const int no_frame, const int nbody, float  * pos , float range[3][2], 
                                       std::string title,std::string sel_comp, std::string filename, 
                                       float timu,bool xy, bool xz, bool yz, bool sview, float * _weight, const int _psort, float * hsml, const int _itf, const bool _wedge, std::string _legend, const int _cmap);
template void C2dplot<double>::compute(std::string pic,const int no_frame, const int nbody, double * pos , float range[3][2], 
                                       std::string title,std::string sel_comp, std::string filename, 
                                       float timu,bool xy, bool xz, bool yz, bool sview, double * _weight, const int _psort, double * hsml, const int _itf, const bool _wedge, std::string _legend, const int _cmap);
// ----------------------------------------------------------------------------
// drawImage
template <class T> void C2dplot<T>::drawImage(const bool disp,const int xaxis, const int yaxis, const int nbview, int &showtext)
{
  std::string outdev = dev;
  std::string label[3] = { "X","Y","Z"};
  static int nplot=0;
  if (disp) {
    if (!sview) {
      // build extension name
      std::string ext="_"+label[xaxis]+label[yaxis]; // example "_XY"
      outdev = buildFrameName(ext,no_frame);
      if (outdev=="?" && filename=="-") {
        nplot=(nplot+1)%nbview;
        std::stringstream ss;
        ss << nplot+1 << "/xs";
        ss >> outdev;
      }
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
  
  // look for indexes in the range axis
  findIndexes(xaxis,yaxis);
  std::cerr << "XYZ particles => Nb indexes found = " << indexes.size() << "\n";

  xmin=ymin=min(range[xaxis][0],range[yaxis][0]);
  xmax=ymax=max(range[xaxis][1],range[yaxis][1]);

  CTimer timing;
#ifndef NO_CUDA // use CUDA computation
  startCuda(nbody,pos,xaxis,yaxis,zmin,zmax);
#else
  // launch computation on threads
  startWorkers(nbody,pos,xaxis,yaxis,zmin,zmax);
#endif
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
  cpgsitf(itf);

  CUtilPgplot cutpg;
  cutpg.selectCMap(cmap);

  cpgimag(tab[0],dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);  
  if (wedge) {
    cpgwedg("BI",4.,5.,zmin,zmax,legend.c_str());
  }
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
  std::string txt=basename + " " + legend;
  cpgmtxt("t",1.8,0.,0.,txt.c_str() );
  
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
template <class T> void C2dplot<T>::startWorkers(const int nbody, T * data, const int xaxis, const int yaxis, float & zmin, float& zmax)
{
  int npart=indexes.size()/nthreads;
  int offset=0;
  // build vector of threads and start them
  for (int x=0; x<nthreads; x++) {
    if (x==nthreads-1) { // last threads
      npart=indexes.size()-offset;
    }
#ifndef NOBOOST
    //std::cerr << "startWorkers: "<<x<<" npart="<<npart<<"\n";
    boost::shared_ptr<boost::thread>  vv(new boost::thread(boost::bind(&C2dplot::worker, this, x,offset,npart,data,xaxis,yaxis)));
    v_threads.push_back(vv) ;
#else
    // this part is not parallel
    worker(x,offset,npart,data,xaxis,yaxis);
#endif
    offset+=npart;
  }
#ifndef NOBOOST
  // wait the end of each threads
  for (std::vector<boost::shared_ptr<boost::thread> >::iterator 
	 it=v_threads.begin();
       it != v_threads.end(); it++) {
    (*it)->join();
  }
#endif
  // accumulate all frames computed in parallel on the first frame
  for (int x=1; x<nthreads; x++) {
    for (int i=0;i<dimy; i++)
      for (int j=0;j<dimx; j++) {
        tab[0][i*dimx+j]+=tab[x][i*dimx+j] ;
      }
  }
  // find zmin zmax
  zmin=numeric_limits<float>::max();
  zmax=-zmin;
  if (0 &&  weight) {
    for (int i=0;i<nbody;i++) {
      zmin = std::min(zmin,(float) weight[i]);
      zmax = std::max(zmax,(float) weight[i]);
    }
  }
  else {
    for (int i=0;i<dimy; i++)
      for (int j=0;j<dimx; j++) {
        zmax=std::max(zmax, tab[0][i*dimx+j]);
        float zzmin=std::min(zmin, tab[0][i*dimx+j]);
        //if (zzmin != 0.0)  zmin=zzmin;
        zmin=zzmin;
      }
    float offset2=0.0;
#if 0
    if (zmin <= 0.0) offset2 = -zmin+.1+numeric_limits<float>::min();
    for (int i=0;i<dimy; i++)
      for (int j=0;j<dimx; j++) {
        //if (tab[0][i*dimx+j]==0.0) tab[0][i*dimx+j]=zmin ;
        tab[0][i*dimx+j] += (zmin+offset2);
      }
#endif
    //std::cerr << "FIRST zmax="<<zmax<<" zmin="<<zmin<<" correct=" << zmin+offset2 << "\n" ;

    // compute LOG
    //zmax=log(zmax+offset2);
    //zmin=log(zmin+offset2);
    zmax=zmax+offset2;
    zmin=zmin+offset2;
    //zmin += 0.003*(zmax-zmin);
    //zmax -= 0.9999999*(zmax-zmin);
    for (int i=0;i<dimy; i++)
      for (int j=0;j<dimx; j++) {
        if (tab[0][i*dimx+j] != 0.0) {
          //tab[0][i*dimx+j] = log(tab[0][i*dimx+j]);///zmax; // normalize
        } else {
           //tab[0][i*dimx+j]=numeric_limits<float>::min();
        }
      }
  }

  std::cerr << "[-->] zmin="<<zmin<<" zmax="<<zmax<<" nbody="<<nbody<< "\n";
}
template void C2dplot<float >::startWorkers(const int nbody, float  * data, const int xaxis, const int yaxis, float& zmin, float& zmax);
template void C2dplot<double>::startWorkers(const int nbody, double * data, const int xaxis, const int yaxis, float& zmin, float& zmax);
// ----------------------------------------------------------------------------
// worker
// function executed in parallel by each threads
template <class T> void C2dplot<T>::worker(const int ithread, const int offset, const int npart, T * data,const int xaxis, const int yaxis)
{
  float * tab_hsml = new float [dimx*dimy];
  // Reset each local array
  for (int i=0;i<dimy; i++)
    for (int j=0;j<dimx; j++) {
      tab[ithread][i*dimx+j] = 0.0;
      tab_hsml[i*dimx+j] = numeric_limits<float>::max();
    }

  //displayIndexes();

  // loop on all sub particles
  // to fill array 
  float zmin=numeric_limits<float>::max();
  float zmax=-zmin;
  //std::cerr << ">>>> zmin" << zmin << " zmax" << zmax << "\n";
  for (int i=0; i<npart; i++) {
    int    ii = indexes[(offset+i)];
    float  xx = data[ii*3+xaxis];
    int     x = ((xx-xmin)/(xmax-xmin))*(dimx-1);
    float  yy = data[ii*3+yaxis];
    int     y = ((yy-ymin)/(ymax-ymin))*(dimy-1);
    //std::cerr << xmax << " " << xmin << " "<< dimx << " "<< hsml[ii] << " " << ((hsml[ii])/(xmax-xmin))*(dimx-1) << "\n";
    float hsml_size;
    if (hsml) {
      hsml_size=ceil(((hsml[ii])/(xmax-xmin))*(dimx-1));
      if (hsml_size>700)
        std::cerr << "hsml_size="<<hsml_size << " hsml="<< hsml[ii] << "\n";
    }

    assert(x<dimx);
    assert(y<dimy);
    if (i==-1) { // first time
        tab[ithread][x*dimx+y] = (1.0*WEIGHT(ii));
    } else {
        switch (psort) {
        case 0: // accumulated
            tab[ithread][x*dimx+y] += (1.0*WEIGHT(ii)); // one more particles into the cell
            if (hsml) tab_hsml[x*dimx+y] = std::min(tab_hsml[x*dimx+y],hsml_size);
            else      tab_hsml[x*dimx+y] = pixel;
            break;
        case 1: // sort max properties
            //std::cerr << "max = " << tab[ithread][x*dimx+y] << "  weight="<<WEIGHT(ii)<<" ii="<<ii<<"\n";
            tab[ithread][x*dimx+y] = std::max((double) tab[ithread][x*dimx+y],(double) WEIGHT(ii));
            if (hsml) tab_hsml[x*dimx+y] = std::min(tab_hsml[x*dimx+y],hsml_size);
            else     tab_hsml[x*dimx+y] = pixel;

            break;
        case 2: // sort min properties
            tab[ithread][x*dimx+y] = std::min((double) tab[ithread][x*dimx+y],(double) WEIGHT(ii));
            break;
        default: // algo error
            assert(0);
            std::cerr << "Should not be here...\n";
            std::exit(0);
        }
    }
    zmin = std::min(zmin,tab[ithread][x*dimx+y]);
    zmax = std::max(zmax,tab[ithread][x*dimx+y]);
  }  

  std::cerr << "Before gaussian Zmin = "<<zmin<< "  zmax = "<<zmax << "\n";
  // loop to store cells from the array which have a weight
  std::map<int, int> HSML;
  pvec.clear();
  for (int i=0;i<dimy; i++) {
    for (int j=0;j<dimx; j++) {
      if (tab[ithread][i*dimx+j] != 0.0) {      
        float myhsml;
        if (hsml) {
          myhsml=tab_hsml[i*dimx+j];
        } else {
          myhsml=pixel;
        }
        CPartProp p(i,j,(float) tab[ithread][i*dimx+j],myhsml);
        pvec.push_back(p);
        HSML[(int) tab_hsml[i*dimx+j]]++;
        // reset array
        tab[ithread][i*dimx+j] = 0.0;
      }
      if  (tab[ithread][i*dimx+j] != 0.0) {
        assert(1);
      }
    }

  }
  // sort vector of particles
  //std::sort(pvec.begin(),pvec.end(),CPartProp::mySort);

  std::cerr << "HSML size =" << HSML.size() << "\n";
  int cpt=0;
  for (std::map<int, int>::iterator ii=HSML.begin(); ii !=HSML.end(); ii++) {

    HSML[(int) ((*ii).first)] = cpt++;
    //int nb= (*ii).second ;
    //std::cerr << (*ii).first << ": " << (*ii).second << "/" <<nb << std::endl;
  }

  for (std::vector <CPartProp>::iterator it=pvec.begin(); it != pvec.end(); ++it) {
    //std::cerr << "->"<<(*it).prop << "\n";
    int    x = (*it).x; // x
    int    y = (*it).y; // y
    float  w = (*it).prop; // weight
    float  h = (*it).hsml; // hsml size

    //std::cerr << "x=" << x << " y=" << y << "\n";
    //std::cerr << "->" << h << "\n";

    // add gaussian on XY coordinates
    //gaussian->applyOnArrayXY(tab[ithread],dimx,dimy,x,y,w);
    //h = std::min((double) h, dimx/1.);
    //h=std::min((double)h,0.01*dimx);
    //gaussian->computeOnArrayXY(tab[ithread],dimx,dimy,x,y,1.,20);
    //h = 10.;
    h = std::min(h,(float)150.);
    if (hsml) {
      gaussian->computeOnArrayXY(tab[ithread],dimx,dimy,x,y,w,h*2.);
    }
    else {
      //tab[ithread][y*dimx+x] = 1.;
      gaussian->applyOnArrayXY(tab[ithread],dimx,dimy,x,y,w,psort);
    }
  }
  delete [] tab_hsml;
}
template void C2dplot<float >::worker(const int ithread, const int offset, const int npart, float  * data,const int xaxis, const int yaxis);
template void C2dplot<double>::worker(const int ithread, const int offset, const int npart, double * data,const int xaxis, const int yaxis);

//
// ----------------------------------------------------------------------------
