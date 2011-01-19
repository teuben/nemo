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

#include <iostream>
#include <cmath>
#include <sstream>
#include <assert.h>
#include "csnaptools.h"

using namespace jclut;

//
// moveToCod
// move particles positions and velocities to center of density
//template <class T> void CSnaptools<T>::moveToCod(const int nbody,T * pos, T * mass)
template <class T> void CSnaptools::moveToCod(const int nbody,T * pos,T * vel, T * mass, T * rho, double cod[6], bool move)
{
  double codp[3] = {0., 0., 0.};
  double codv[3] = {0., 0., 0.};
  double w_i, w_sum=0.0;
  // loop on all the bodies
  for (int i=0; i<nbody;i++) {
    w_i    = rho[i] * mass[i]; // weight = rho * mass
    w_sum += w_i;              // sum
    if (pos) {
      codp[0] +=(pos[i*3  ]*w_i);
      codp[1] +=(pos[i*3+1]*w_i);
      codp[2] +=(pos[i*3+2]*w_i);
    }
    if (vel) {
      codv[0] +=(vel[i*3  ]*w_i);
      codv[1] +=(vel[i*3+1]*w_i);
      codv[2] +=(vel[i*3+2]*w_i);
    }    
  }
  assert(w_sum>0.0); // total weight must be positif
  if (pos) {
    codp[0] /= w_sum;
    codp[1] /= w_sum;
    codp[2] /= w_sum;
  }
  cod[0] = codp[0]; cod[1] = codp[1]; cod[2] = codp[2];
  if (vel) {
    codv[0] /= w_sum;
    codv[1] /= w_sum;
    codv[2] /= w_sum;
    
  }
  cod[3] = codv[0]; cod[4] = codv[1]; cod[5] = codv[2];
  std::cerr << "COD = " << cod[0]<<" "<< cod[1]<<" "<< cod[2]<<" "
      << cod[3]<<" "<< cod[4]<<" "<< cod[5]<<"\n";
  if (move) {
    // move to cod requested
    for (int i=0; i<nbody;i++) {
      for (int j=0;j<3;j++) {
        if (pos)
          pos[i*3+j] -= codp[j]; // pos
        if (vel)
          vel[i*3+j] -= codv[j]; // vel
      }
    }
  }
}
template void CSnaptools::moveToCod(const int nbody,float * pos,float * vel, float * mass, float * rho, double cod[6], bool move);
template void CSnaptools::moveToCod(const int nbody,double * pos,double * vel, double * mass, double * rho, double cod[6], bool move);
//
// moveToCom
// move particles positions to center of mass
//template <class T> void CSnaptools<T>::moveToCom(const int nbody,T * pos, T * mass)
template <class T> void CSnaptools::moveToCom(const int nbody,T * pos, T * mass)
{
  double com[3] = {0., 0., 0.};
  double np=0.,masstot=0;;
  for (int i=0; i<nbody;i++) {
    float massi;
    if (mass) massi = mass[i];
    else      massi = 1.0;
    masstot+=massi;
    np++;
    int jndex= i;
    com[0] +=(pos[jndex*3  ]*massi);
    com[1] +=(pos[jndex*3+1]*massi);
    com[2] +=(pos[jndex*3+2]*massi);
  }
  if (!mass) {
    std::cerr << "No mass in the snapshot, we assum mass=1.0 for each particles...\n";
  }
  std::cerr <<"COM     ="<<com[0]<<" "<<com[1]<<" "<<com[2]<<"\n";
  std::cerr <<"np      ="<<np<<"\n";
  std::cerr <<"mass tot="<<masstot<<"\n";
  // center
  for (int i=0; i<nbody;i++) {
    pos[i*3+0] -= (com[0]/masstot);
    pos[i*3+1] -= (com[1]/masstot);
    pos[i*3+2] -= (com[2]/masstot);
  }  
}
template void CSnaptools::moveToCom<float>(const int nbody,float * pos, float * mass);
template void CSnaptools::moveToCom<double>(const int nbody,double * pos, double * mass);
//
// basename
std::string  CSnaptools::basename(const std::string str)
{
  size_t found=str.find_last_of("/\\");
  return str.substr(found+1);
}

//
// dirname
std::string CSnaptools::dirname(const std::string str)
{
  size_t found=str.find_last_of("/\\");
  return str.substr(0,found);
}
//
// parseString
std::string CSnaptools::parseString(std::string & next_string, const std::string sep)
{
  std::string return_string;
  std::string::size_type coma=next_string.find(sep,0);  // try to find ","
  if (coma != std::string::npos) { // found "separator"
    return_string = next_string.substr(0,coma);
    next_string   = next_string.substr(coma+1,next_string.length());
  } else {                         // not found
    return_string = next_string;
    next_string = "";
  }
  return return_string;
  
}

//
// stringToVector
template <class T>  std::vector<T> CSnaptools::stringToVector(const std::string s, const int min, T val, std::string sep)
{
  std::string current_s,next_s;
  next_s = s;              // string to be parsed
  
  std::vector <T> vec;
  T value;
  // parse 
  while ((current_s=parseString(next_s,sep)) != "") {  
    std::stringstream parse; // string parsed
    parse << current_s;      // read value
    parse >> value;          // convert value
    vec.push_back(value);
  }
  // complete to default value if size < min
  for (int i=vec.size(); i<min; i++) {
    vec.push_back(val); // default value;
  }  
  return vec;
}
template std::vector<float> CSnaptools::stringToVector<float>(const std::string s, const int min, float val, std::string sep);
template std::vector<int  > CSnaptools::stringToVector<int  >(const std::string s, const int min, int   val, std::string sep);
//
