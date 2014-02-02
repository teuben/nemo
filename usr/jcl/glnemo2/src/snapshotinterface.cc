// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2014                                  
// e-mail:   Jean-Charles.Lambert@lam.fr                                      
// address:  Centre de donneeS Astrophysique de Marseille (CeSAM)              
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include "snapshotinterface.h"
#include <sstream>
#include <cmath>

namespace glnemo {

  // ============================================================================
  // parseSelectTime                                                             
  void SnapshotInterface::parseSelectTime()
  {
    std::string current_s,next_s;
    next_s = select_time;
    while ((current_s=parseString(next_s)) != "") {  // || 1 force parsing
      getRangeTime(current_s);
    }
  }
  // ============================================================================
  // parseString                                                                 
  // return the string at the position after the next 'coma' otherwise ""        
  std::string SnapshotInterface::parseString(std::string & next_string)
  {
    std::string return_string;
    std::string::size_type coma=next_string.find(",",0);  // try to find "'"
    if (coma != std::string::npos) { // found coma
      return_string = next_string.substr(0,coma);
      next_string   = next_string.substr(coma+1,next_string.length());
    } else {                         // not found
      return_string = next_string;
      next_string = "";
    }
    return return_string;
  }
  // ============================================================================
  // getRangeTime                                                                
  // split a string a:b[:c] into inf, sup,delta                                  
  void SnapshotInterface::getRangeTime(const std::string rtime)
  {
    int status;
    std::vector<float> store;
    int ppos=0;
    bool stop=false;
    int cpt=0;
    while (! stop) {
      size_t found = rtime.find(':',ppos);
      if (found!=std::string::npos) {
        if (found > (size_t) (ppos)) {
          cpt++;
          std::string str=rtime.substr(ppos,found-ppos);
          std::istringstream ss(str);
          float val;
          ss>>val;
          store.push_back(val);
        }
        ppos=found+1; //
      } else { // no more ":"
        std::string str=rtime.substr(ppos);
        if (str == "all") {
          store.push_back(-1.0);
        }
        else {
          std::istringstream ss(str);
          float val;
          ss>>val;
          store.push_back(val);
        }
        stop=true;
      }
    }
    for (std::vector<float>::iterator it=store.begin(); it!=store.end(); it++) {
      //std::cerr << "SnapshotInterface::getRangeTime range i ="<< it-store.begin() << "\n";
      //std::cerr << "SnapshotInterface::getRangeTime value   =" << *it << "\n";
    }
    float inf=store[0];
    float sup=inf;
    float offset=0.0;
    if (store.size()>1) {
      sup = store[1];
    }
    if (store.size()>2) {
      offset = store[2];
    }
    
    assert(sup>=inf);
    status=0;
    // store inf and supp time in vector array
    CSelectTime st(inf,sup,offset,-666.0);
    stv.push_back(st);
  }
  // ============================================================================
  // checkRangeTime                                                              
  // check if current time belong to the interval of the selected time           
  bool SnapshotInterface::checkRangeTime(const float current_time)
  {
    assert(stv.size()>0);
    for (CSelectTimeVector::iterator it=stv.begin(); it!=stv.end(); it++) {
      //std::cerr << "t inf="<<it->inf<<" t sup="<<it->sup
      //      << " c t="<<current_time <<"\n";
      if (it->inf == -1. || it->sup == -1. ||   // "all" time
          ( current_time >= it->inf && current_time <= it->sup)) { // in the range
        bool status=true;
        if (it->offset > 0.) {     // there is an offset time 
          if (it->lastt != 666.) { // it's not the first stime
            //std::cerr << std::setprecision(9) << "ct = " << current_time << " l+o = " <<it->lastt+it->offset<<" fabsdiff = "<< fabs(current_time-it->lastt-it->offset)<<" offset = "<<it->offset<<"\n"; 
            if ( current_time >= (it->lastt+it->offset) || 
                 diffTime(current_time-it->lastt-it->offset)) { // time match!!
              it->lastt = current_time; // save current time
            } else { // time does not match
              status=false;
            }
            
          } else {                 // first time
            it->lastt = current_time; // save current time
          }
        }
        return status;;
      }
    }
    return false;
  }

}
