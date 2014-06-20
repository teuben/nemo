// ============================================================================
// Copyright Jean-Charles LAMBERT - 2010-2014
//           Centre de donneeS Astrophysiques de Marseille (CeSAM)
// e-mail:   Jean-Charles.Lambert@lam.fr                                      
// address:  Aix Marseille Universite, CNRS, LAM 
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS UMR 7326                                       
// ============================================================================
#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <assert.h>
#include "tipsydefs.std.h"
#include "defs.h"
#include "xdrfuncs.h"
#include <string>
#include <nemo.h>
#include <uns.h>

// ------------------------------------------------------------
// Nemo variable
const char * defv[] = {  // use `::'string because of 'using namespace std'
  "in=???\n           tipsy input file ",
  "out=???\n          gadget2 output file)",
  "scale=1.0\n      scaling factor (pos,hsml)",
  "verbose=f\n        verbose mode                    ",
  "VERSION=1.0\n      compiled on <"__DATE__"> JCL   ",
  NULL
};
const char * usage="Simple converter tipsy to gadget2";

void readtipsystd_(const char * filename,std::string outname,float scale)
{
  XDR xdr;
  FILE *fp;
  int i,j;
  //  char filename[80];
  struct star_particle sp;   
  struct dark_particle dp;
  struct gas_particle gp;
  struct dump head;
  
  uns::CunsOut * unsout = new uns::CunsOut(outname,"gadget2",false);
 
  fp = fopen(filename,"r");
  xdrstdio_create(&xdr,fp,XDR_DECODE);
  
  /* read header */
  xdr_header(&xdr,&head);
  std::cerr << "Scale = "<<scale <<"\n";
  /*printf("  Particles in TIPSY file: %d gas, %d dark, %d star\n",
    head.nsph,head.ndark,head.nstar);*/
  assert(head.ndim==3);
  float time = head.time;
  //int ndim=head.ndim;
  int nbods = head.nbodies;
  int ndark=head.ndark;
  int ngas = head.nsph;
  int nstar=head.nstar;
  /*printf(" Time: %g\n",*time);*/

  fprintf(stderr,"Time=%f nbods=%d nhalo=%d ngas=%d nstars=%d\n",
	  time,nbods,ndark,ngas,nstar);
  /* reads data for gas particles */
  float * mass=NULL, * xx=NULL, * vv=NULL, * rho=NULL, *hsml=NULL; 
  bool ok;

  // set time
  ok=unsout->setData("time",time);
  if (ok) ;  // remove compiler warning

  /* reads data for gas particles */
  i = 0;
  while (i<head.nsph) {
    xdr_gas(&xdr,&gp);

    if (!i) {
      mass = new float[head.nsph];
      xx   = new float[head.nsph*3];
      vv   = new float[head.nsph*3];
      rho  = new float[head.nsph];
      hsml = new float[head.nsph];
    }

    mass[i]  =gp.mass;
    xx[i*3+0]=gp.pos[0]*scale;
    xx[i*3+1]=gp.pos[1]*scale;
    xx[i*3+2]=gp.pos[2]*scale;
    vv[i*3+0]=gp.vel[0];
    vv[i*3+1]=gp.vel[1];
    vv[i*3+2]=gp.vel[2];
    rho[i]   =gp.rho;
    //temp[i]=gp.temp;
    hsml[i]=gp.hsmooth*scale;
    //metals[i]=gp.metals;
    //pot[i]=gp.phi;

    ++i;
  }

  if (head.nsph) {
    ok=unsout->setData("gas","pos" ,head.nsph,xx  ,false);
    ok=unsout->setData("gas","vel" ,head.nsph,vv  ,false);
    ok=unsout->setData("gas","mass",head.nsph,mass,false);
    ok=unsout->setData("gas","rho" ,head.nsph,rho ,false);
    ok=unsout->setData("gas","hsml",head.nsph,hsml,false);
    delete [] mass;
    delete [] xx;
    delete [] vv;
    delete [] rho;
    delete [] hsml;
  }
  /* reads data for dark particles */

  j=0;
  i = 0;
  while (i<head.ndark) {
    xdr_dark(&xdr,&dp);
    
    if (!i) {
      mass = new float[head.ndark];
      xx   = new float[head.ndark*3];
      vv   = new float[head.ndark*3];
    }
    
    mass[i]=dp.mass;
    xx[i*3+0]=dp.pos[0]*scale;
    xx[i*3+1]=dp.pos[1]*scale;
    xx[i*3+2]=dp.pos[2]*scale;
    vv[i*3+0]=dp.vel[0];
    vv[i*3+1]=dp.vel[1];
    vv[i*3+2]=dp.vel[2];
    
    ++i;
    ++j;
  }
  
  if (head.ndark) {
    ok=unsout->setData("halo","pos" ,head.ndark,xx  ,false);
    ok=unsout->setData("halo","vel" ,head.ndark,vv  ,false);
    ok=unsout->setData("halo","mass",head.ndark,mass,false);
    delete [] mass;
    delete [] xx;
    delete [] vv;
  }

  /* reads data for star particles */
  i=0;
  while (i<head.nstar) {
    xdr_star(&xdr,&sp);

    if (!i) {
      mass = new float[head.nstar];
      xx   = new float[head.nstar*3];
      vv   = new float[head.nstar*3];
    }
    
    mass[i]=sp.mass;
    xx[i*3+0]=sp.pos[0]*scale;
    xx[i*3+1]=sp.pos[1]*scale;
    xx[i*3+2]=sp.pos[2]*scale;
    vv[i*3+0]=sp.vel[0];
    vv[i*3+1]=sp.vel[1];
    vv[i*3+2]=sp.vel[2];

    ++i;
 
  }

  if (head.nstar) {
    ok=unsout->setData("stars","pos" ,head.nstar,xx  ,false);
    ok=unsout->setData("stars","vel" ,head.nstar,vv  ,false);
    ok=unsout->setData("stars","mass",head.nstar,mass,false);
    delete [] mass;
    delete [] xx;
    delete [] vv;
  }

  fclose(fp);

  //printf("Data for %d particles read from %s \n",i,filename);

  unsout->save();

}

// ------------------------------------------------------------
// main program
int main(int argc, char ** argv) { 
  //   start  NEMO
  initparam(const_cast<char**>(argv),const_cast<char**>(defv));
  if (argc) {;} // remove compiler warning :)
  
  // Get input parameters
  std::string simname   = (getparam ((char *) "in"      ));
  std::string outname   = (getparam ((char *) "out"     ));
  float       scale     = (getdparam((char *) "scale"   ));
  std::cerr << "Scale = "<<scale <<"\n";
  readtipsystd_(simname.c_str(),outname,scale);

  //   finish NEMO
  finiparam();
}
// ------------------------------------------------------------
		 
		 
		 






