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
#include "ctools.h"
#include <cstdlib>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <algorithm>
namespace tools {
// ----------------------------------------------------------------------------
// isFileExist                                                                 
bool Ctools::isFileExist(std::string _file)
{
  bool status=false;
  std::ifstream ftest;
  ftest.open(_file.c_str(),std::ios::in);
  if ( ftest.is_open()) {
    status=true;
    ftest.close();
  }
  return status;
}
// ----------------------------------------------------------------------------
// fixFortran
std::string Ctools::fixFortran(const char * _ff, const int len, const bool lower)
{
  char * buff = new char[len+1];
  strncpy(buff,_ff,len);
  buff[len]='\0';
  std::string str=buff;
  
  delete [] buff;
  
  size_t found;
  // for backward compatibility
  // correct old program with explicit '\0' character
  found = str.find("\\");
  if (found!=std::string::npos) {
    //std::cerr << "FOUND at "<<found<< "\n";
    str.replace(found,2," ");
  }    
  
  found=str.find_last_not_of(" ");
  if (found!=std::string::npos)
    str.erase(found+1);
  else
    str.clear();            // str is all whitespace
  //std::cerr << "fix_fortran2 =["<<str<<"]\n";
  
  return str;
}
// ----------------------------------------------------------------------------
// fixFortran
std::string Ctools::fixFortran(const char * _ff, const bool lower)
{
  static char buff[200], * p;
  memset( buff, '\0', 200 );

  //std::cerr << "Fortran string ["<<_ff<<"]\n";

  p=(char *) strchr(_ff,'\\');
  if (p) {
    //std::cerr << "Got \\ \n";
    assert (p-_ff<=200);
    strncpy(buff,_ff,p-_ff);
  }
  else {
    p=(char *) strchr(_ff,'#');
    if (p) {
      //std::cerr << "Got #\n";
      assert (p-_ff<=200);
      strncpy(buff,_ff,p-_ff);
    } else {
      //std::cerr << "Got nothing.....\n";
      strcpy(buff,_ff);
    }
  }
  //std::cerr << "Buff ["<<buff<<"]\n";
  if (lower)
    return tolower(std::string(buff));
  else
    return std::string(buff);
}
// ----------------------------------------------------------------------------
// tolower
std::string Ctools::tolower(std::string s)
{
  std::transform(s.begin(),s.end(),s.begin(),(int(*)(int)) std::tolower);
  return s;
}
// ----------------------------------------------------------------------------
// tolupper
std::string Ctools::toupper(std::string s)
{
  std::transform(s.begin(),s.end(),s.begin(),(int(*)(int)) std::toupper);
  return s;
}
}
