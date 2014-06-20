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
#include <iostream>                                   // C++ I/O     
#include <fstream>                                    // C++ file I/O
#include <sstream>
#include <cstdio>  
#include <cstdlib> 
#include <assert.h>
#include <cmath>
#include <nemo.h>
#include <iomanip>

// ------------------------------------------------------------
// Include file
#include <cpgplot.h>
#include "uns.h"
#include "csnaptools.h"
#include "c2dplot.h"
using namespace jclut;
using namespace uns_proj;
using namespace std;
// ------------------------------------------------------------
// Nemo variable
const char * defv[] = {
  "in=???\n		                UNS input snapshot",
  "select=???\n               select particles (range, or component name)",
  "prop=\n                    properties (metal,age,temp,rho)",
  "range=\n                   -range < x|y|z > +range",
  "xrange=-35.0:35.0\n        x scaling",
  "yrange=-35.0:35.0\n        y scaling",
  "zrange=-35.0:35.0\n        z scaling",
  "psort=0\n                  sort properties (0:none, 1:max, 2:min)",
  "xy=t\n                     display XY",
  "xz=t\n                     display XZ",
  "zy=f\n                     display ZY",
  "dev=?\n                    output device, ? will ask you, name will save in gif format",
  "com=t\n                    center according to com",
  "hsml=f\n                   use hydro smooth length if it exist ?",
  "cb=f\n                     toggle color bar",
  "cmap=0\n                   colormap, 0:rainbow, 1:heat, 2:gray",
  "sview=t\n                  one single view for all display",
  "title=\n                   simulation title",
  "pfname=t\n                 boolean, print filename ? (default t)",
  "times=all\n		            selected time",
  "no=0\n                     output file's index",
  "pixel=20\n                 size in pixel of the gaussian",
  "dimx=1024\n                internal image size",
  "dimy=1024\n                internal image size",
  "itf=1\n                    image transer function (0:linear, 1:logarithmic, 2:square-root",
  "gp=5.\n                    gaussian parameter",
  "threads=1\n                #threads used",
  "verbose=f\n                verbose on/off",
  "VERSION=2.O\n              compiled on <"__DATE__"> JCL  ",
  NULL,
};
const char * usage="2D Plot of a UNS file, y=f(x) and/or z=f(x) and/or y=f(z)";

int nbody,iter;

// ------------------------------------------------------------
// setrange :
// Conver string like "-2.0:1.5" in 2 reals :
// range[0]=-2.0 et range[1]=1.5
void setrange(float range[],const char * ch_range)
{
  char * p;  
  p = strchr((char *) ch_range,':');
  if (p) {
    range[0] = atof(ch_range);
    range[1] = atof(p+1);
  }
  else {
    range[0] = 0.0;
    range[1] = atof(ch_range);
  }
}
// ------------------------------------------------------------
// main program
int main(int argc, char ** argv )
{
  std::string   outgif;
  std::string   sitf[]={ "linear", "log", "square" };
  bool   first = true;  
  float range[3][2];
  //   start  NEMO
  initparam(const_cast<char**>(argv),const_cast<char**>(defv));
  
  // get input parameters
  char * simname     = getparam ((char *) "in"      );
  char * select_c    = getparam ((char *) "select"  );
  std::string prop   = getparam ((char *) "prop"    );
  char * outfile     = getparam ((char *) "dev"     );
  char * select_t    = getparam ((char *) "times"   );
  char * title       = getparam ((char *) "title"   );
  bool   printfname  = getbparam((char *) "pfname"  );
  int    psort       = getiparam((char *) "psort"   );
  bool   com         = getbparam((char *) "com"     );
  bool   use_hsml    = getbparam((char *) "hsml"    );
  bool   wedge       = getbparam((char *) "cb"      );
  int    cmap        = getiparam((char *) "cmap"    );
  bool   xy          = getbparam((char *) "xy"      );
  bool   xz          = getbparam((char *) "xz"      );
  bool   zy          = getbparam((char *) "zy"      );
  bool   sview       = getbparam((char *) "sview"   );
  bool   verbose     = getbparam((char *) "verbose" );
  int    pixel       = getiparam((char *) "pixel"   );
  int    dimx        = getiparam((char *) "dimx"    );
  int    dimy        = getiparam((char *) "dimy"    );
  int    itf         = getiparam((char *) "itf"     );
  float  gp          = getdparam((char *) "gp"      );
  int    no_frame    = getiparam((char *) "no"      );
  int    threads     = getiparam((char *) "threads" );
  std::string rr     = getparam ((char *) "range"   );
  
  if (rr.length()!=0) {
    std::string srange="-"+rr+":"+rr;
    setrange(range[0],srange.c_str());
    setrange(range[1],srange.c_str());
    setrange(range[2],srange.c_str());
  } 
  else {
    setrange(range[0],getparam((char *)"xrange"));
    setrange(range[1],getparam((char *)"yrange"));
    setrange(range[2],getparam((char *)"zrange"));
  }    
  // instantiate a new uns object
  uns::CunsIn * uns = new uns::CunsIn(simname,select_c,select_t,verbose);  
  // C2dplot object
  C2dplot<float> * c2dplot = new C2dplot<float>(threads,pixel,dimx,dimy,gp);
    
  if (uns->isValid()) {
    while(uns->snapshot->nextFrame()) {
      std::stringstream legend;
      bool ok;
      int cnbody,nbody;      
      float * pos, * mass, time;
      // get the input number of bodies according to the selection
      ok=uns->snapshot->getData("nsel",&nbody);
      // get the simulation time
      ok=uns->snapshot->getData("time",&time);
      // get POS from input snapshot
      ok=uns->snapshot->getData("pos" ,&cnbody,&pos);
      if (!ok) {
        std::cerr << "No positions, aborted !\n";
        std::exit(1);
      }
      // get MASS from input snapshot
      ok=uns->snapshot->getData("mass",&cnbody,&mass);
      if (!ok) {
        std::cerr << "No masses, aborted !\n";
        std::exit(1);
      }
      // get properties
      float * weight=NULL;
      ok = false;
      std::string comp_prop=select_c;
      if (prop != "") {
        ok=uns->snapshot->getData(select_c,prop,&cnbody,&weight);
      }
      if (!ok) {
        std::cerr << "No propertie["<<prop<<"] for the selection\n";
        //std::exit(1);
      } else {
        if (prop=="age") {
          legend << "Time of birth ";
        } else {
          comp_prop += " " + prop;
          legend << prop << " ";
        }
        std::cerr << "Properties = " << prop << "\n";
      }
      // get hsml
      float * hsml=NULL,maxhsml,minhsml;
      ok = false;
      if (use_hsml) {
        ok=uns->snapshot->getData(select_c,"hsml",&cnbody,&hsml);
        if (!ok) {
          std::cerr << "No HSML for the selection\n";
          //std::exit(1);
        } else {
          legend << "H ";
          maxhsml=minhsml=hsml[0];
          for (int i=0; i<nbody; i++) {
            maxhsml = std::max(hsml[i],maxhsml);
            minhsml = std::min(hsml[i],minhsml);

          }
          std::cerr << "max HSML =  " << maxhsml << "\n";
          std::cerr << "min HSML =  " << minhsml << "\n";
        }
      }
      std::string filename=uns->snapshot->getFileName();;
      std::cerr << "nbody=" << nbody << " time="<< time <<"\n";
      std::cerr << "filename = " << filename << "\n";
      if (psort) legend << "sort ";
      else       legend << "add ";
      if (itf>=0 && itf<=3) legend << sitf[itf];
      std::cerr << "LEGEND ="<< legend.str() << "\n";
      if (!printfname) filename = ""; 
      
      if (com) CSnaptools::moveToCom<float>(nbody,pos,mass); // COM centering
      
      /* call engine */
      if (nbody >0) c2dplot->compute(outfile,no_frame,nbody,pos,range,title,
                                     comp_prop,filename,time,xy,xz,zy,sview, weight, psort, hsml, itf, wedge, legend.str(),cmap);
      
      /* iteration */
      no_frame++;      
    }
    if (! first)
      cpgend();    
    //   finish NEMO
    finiparam();
    return 1;
  } 
}	      
//
