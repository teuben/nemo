// =====================================================================
// Copyright Jean-Charles LAMBERT - 2012-2013
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies       
//           CeSAM (Centre de donneeS Astrophysique de Marseille)                                     
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// =====================================================================
// 
// The following program shows how to use UNSIO library from
// c++ program
//
// This program reads an unsio compatible snapshot from the command line
// and save it in gadget2 format
//
// Syntaxe : unsio_fortran myinput myoutput select_comp select_time
//
// myinput     -> an unsio compatible input snapshot
// myoutput    -> output filename
// select_comp -> the component to be saved
// select_time -> a range of time which select input snapshot
//
// For more information about how to use UNSIO, visit:
// http://projets.oamp.fr/projects/unsio/
//
// =====================================================================
#include <iostream>
#include <cstdlib>
#include <string>
#include <uns.h>

using namespace std;

int main(int argc, char ** argv)
{

  if (argc < 3) {
    cerr << "\nsyntaxe: " << argv[0] << " filename component\n";
    cerr << "example: "   << argv[0] << " list all\n\n";
    exit(1);
  }

  string simname  (argv[1]);
  string component(argv[2]);
  bool verbose=false;

  // -----------------------------------------------
  // instantiate a new UNS input object (for reading)
  uns::CunsIn * uns = new uns::CunsIn(simname,component,"all",verbose);

  if (uns->isValid()) { // input file is known by UNS lib
    while(uns->snapshot->nextFrame()) { // there is a new frame
      string stype = uns->snapshot->getInterfaceType();
      cout << "File name : " << uns->snapshot->getFileName()<<"\n";
      cout << "File type : " << stype                       <<"\n";

      int nbody; float time;
      // get the input number of bodies according to the selection
      uns->snapshot->getData("nsel",&nbody);
      // get the simulation time
      uns->snapshot->getData("time",&time);

      // read position
      float * pos;
      bool ok = uns->snapshot->getData("all","pos",&nbody,&pos); // we get positions here
                                                                 // in one dimensionnal array
      if (ok) {
        for (int i=0; i<nbody; i++) {
          // display positions x,y,z
          cerr << "x="<<pos[i*3+0] << " y="<<pos[i*3+1] << " z="<<pos[i*3+2] <<"\n" ;
        }
      }
    }
  }
  delete uns;
}

//
