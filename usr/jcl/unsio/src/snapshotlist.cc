// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008-2013                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

/**
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */
#include "uns.h"
#include "snapshotlist.h"
#include "snapshotgadget.h"
#include "snapshotnemo.h"
#include "ctools.h"
#include <sstream>
#include <iomanip>
#include <iostream>

namespace uns {
  // ----------------------------------------------------------------------------
  // constructor
  CSnapshotList::CSnapshotList(const std::string _name,
                               const std::string _comp, 
                               const std::string _time,
                               const bool        verb)
                                 :CSnapshotInterfaceIn(_name, _comp, _time, verb)
  {
    snapshot = NULL;
    unsin    = NULL;
    nframe   = 0;   // # frames read
    nemosim  = "";
    valid=openFileList();
  }
  // ----------------------------------------------------------------------------
  // constructor
  CSnapshotList::~CSnapshotList()
  {
    if (unsin) delete unsin;
  }
  // ============================================================================
  // getSnapshotRange 
  int CSnapshotList::nextFrame(uns::UserSelection &user_select)
  {
    assert(snapshot != NULL);
    assert(snapshot->isValidData()==true);
    snapshot->setNsel(nsel);   
    return (snapshot->nextFrame(user_select));
  }
  // ============================================================================
  // nextFrame 
  int CSnapshotList::nextFrameSelect(ComponentRangeVector * crvs)
  {
    snapshot->user_select.setSelection(getSelectPart(),crvs);
    setNsel(snapshot->user_select.getNSel());
    snapshot->setReqBits(req_bits);
    snapshot->setNsel(snapshot->user_select.getNSel());
    return(snapshot->nextFrame(snapshot->user_select));
  }
  // ============================================================================
  // isNewFrame() 
  bool CSnapshotList::isNewFrame()
  {
    bool stop=false;
    while (!stop && getLine()) {  
      if (unsin) {
        delete unsin;
      }
      unsin = new uns::CunsIn(snapname.c_str(),select_part.c_str(),select_time.c_str(),verbose);
      float t;
      bool ok=unsin->snapshot->getData("time",&t);
      if (unsin->isValid() && ok && checkRangeTime(t)) {
        snapshot = unsin->snapshot;
        stop=true;        
        interface_type = snapshot->getInterfaceType();
      }
    }
    if (!stop) end_of_data = true; // EOF reached
    return stop;
    
  }
  // ============================================================================
  ComponentRangeVector * CSnapshotList::getSnapshotRange()
  {
    assert(snapshot != NULL);
    assert(snapshot->isValidData());
      if ((tools::Ctools::tolower(simtype) == "nemo") && nemosim != "" && crv.size()>0) {
      return &crv;
    } 
    else {
      return snapshot->getSnapshotRange();
    }
  }
  // ----------------------------------------------------------------------------
  // openFileList
  bool CSnapshotList::openFileList()
  {
    bool status=false;
    if (filename == "-" ) 
      ;
    else 
      fi.open(filename.c_str(),std::ios::in);
    if (! fi.is_open()) {
      std::cerr << "Unable to open file ["<<filename<<"] for reading, aborting...\n";
      status = false;
    }
    else {
      // read magic numberé
      std::string line;
      if (getLine(true)) { // read the first file
        // -----------------------------------------------
        // instantiate a new UNS input object (for reading)
        uns::CunsIn * test_data = new uns::CunsIn(snapname.c_str(),select_part.c_str(),select_time.c_str(),verbose);
        
        if (test_data->isValid()) { // it's a valid snaphot
          delete test_data;
          status = true;
          fi.seekg(0, std::ios::beg); // go back to the beginning
        }
      } 
      else {        
        status=false;
        fi.close();
      }
    }
    return status;
  }
  // ----------------------------------------------------------------------------
  // 
  bool CSnapshotList::getLine(const bool force)
  {
    bool status=false,stop=false;
    if (valid || force) {
      while (!stop && ! fi.eof()) {           // while ! eof
        std::string line;
        getline(fi,line);
        if ( ! fi.eof()) {
          std::istringstream str(line);  // stream line
          std::string parse;
          // following loop parse each lines previously read   
          //
          int cpt=0;
          while (  str >> parse   &&              // something to read 
                   parse[0] != '#' &&             // not commented out 
                   parse[0] != '!' &&             // not commented out 
                   parse[0] != '\n'               // not a blank line
                   ) {
            cpt++;
            if (cpt==1) snapname=parse;
          }
          if (cpt > 0 ) {
            unsigned int i=0;
            while(i<snapname.length() && snapname[i]==' ') i++; // search first non blank
            if (i<snapname.length() && snapname[i]!='/')        // first char not a '/'  
            {;}//snapname = dirpath.toStdString() + snapname;      // append to dirpath     
            stop   = true; // we have a snapname
            status = true; // so we can stop reading            
          }
        }
        else { // end of file
          stop   = true;
          status = false;
        }
      }
    }
    else status=false;
    return status;
  }
}


//
