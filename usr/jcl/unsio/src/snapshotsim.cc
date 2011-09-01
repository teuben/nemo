// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008-2010                                   
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
#ifndef NOSQLITE3  // do not compite if no sqlite3 lib
#include "snapshotsim.h"
#include "snapshotgadget.h"
#include "snapshotnemo.h"
#include "ctools.h"
#include <sstream>
#include <iomanip>
#include <iostream>

#define DEBUG 0
#include "unsdebug.h"
#include "uns.h"

namespace uns {
  // ASCII database
  std::string uns::CSnapshotInterfaceIn::sim_db_file="/pil/programs/DB/sim_info.txt";
  std::string uns::CSnapshotInterfaceIn::nemo_range_file="/pil/programs/DB/nemo_range.txt";
  std::string uns::CSnapshotInterfaceIn::eps_db_file="/pil/programs/DB/sim_eps.txt";
  // SQLITE database
  
// ----------------------------------------------------------------------------
// constructor
CSnapshotSimIn::CSnapshotSimIn(const std::string _name,
			   const std::string _comp, 
			   const std::string _time,
			   const bool        verb)
  :CSnapshotInterfaceIn(_name, _comp, _time, verb)
{
  snapshot = NULL;
  sql      = NULL;
  nframe   = 0;   // # frames read
  nemosim  = "";
  verbose  = verb;
  if (0) {
    valid=openDbFile();
  } else { // SQLite3
    valid=openSqlDb();
  }
}
// ----------------------------------------------------------------------------
// constructor
CSnapshotSimIn::~CSnapshotSimIn()
{
  if (snapshot) delete snapshot;
  if (sql) delete sql;
}
// ============================================================================
// getSnapshotRange                                                            
ComponentRangeVector * CSnapshotSimIn::getSnapshotRange()
{
  assert(snapshot != NULL);
  assert(snapshot->isValidData());
    if ((simtype == "Nemo") && nemosim != "" && crv.size()>0) {
    return &crv;
  } 
  else {
    return snapshot->getSnapshotRange();
  }
}
// ============================================================================
// getSnapshotRange 
int CSnapshotSimIn::nextFrame(uns::UserSelection &user_select)
{
  assert(snapshot != NULL);
  assert(snapshot->isValidData()==true);
  snapshot->setNsel(nsel);
  return (snapshot->nextFrame(user_select));
}
// ============================================================================
// nextFrame 
int CSnapshotSimIn::nextFrameSelect(ComponentRangeVector * crvs)
{
  snapshot->user_select.setSelection(getSelectPart(),crvs);
  setNsel(snapshot->user_select.getNSel());
  snapshot->setReqBits(req_bits);
  snapshot->setNsel(snapshot->user_select.getNSel());
  return(snapshot->nextFrame(snapshot->user_select));
}
//                     - - - - - - - - - - - - - - 
//                           SQlite database       
//                     - - - - - - - - - - - - - - 

// ============================================================================
// openSqlDb                                                                   
bool CSnapshotSimIn::openSqlDb(std::string db)
{
  sqlite_db = db;
  sql = new jclt::CSQLite3(sqlite_db);
  bool status=sql->isOpen();
  if (! status) {
    //std::cerr << __FILE__<< " " << __LINE__ << "Aborting ....\n";
    //std::exit(1);
  } else {
    status = findSqlSim();
    if (status) {
      eps_exist = readSqlEps();
    } else {
      eps_exist = false;
    }
  }
  return status;
}
// ============================================================================
// findSqlSim                                                                  
bool CSnapshotSimIn::findSqlSim()
{
  std::string select="select * from info where name='"+filename+"'";
  if (verbose) std::cerr << "select = "<<select <<"\n";
  int status = sql->exe(select);
  if (status) {
    if (verbose) sql->display();
    assert(sql->vdata[0]==filename);
    simname        = sql->vdata[0];
    simtype        = sql->vdata[1];
    dirname        = sql->vdata[2];
    basename       = sql->vdata[3];
    
    interface_type = simtype;
    if (interface_type == "Gadget") interface_index=1;
    else
      if (interface_type == "Nemo") interface_index=0;
      else {
	std::cerr <<"CSnapshotSimIn::findSqlSim => Unknown interface type....\n";
      }
  }
  return status;
}
// ============================================================================
// readSqlEps                                                                  
bool CSnapshotSimIn::readSqlEps()
{
  std::string select="select * from eps where name='"+filename+"'";
  if (verbose) std::cerr << "select = "<<select <<"\n";
  int status = sql->exe(select);
  if (status) {
    if (verbose) sql->display();
    assert(sql->vdata[0]==filename);
    std::stringstream str;
    for (unsigned int i=1; i< sql->vdata.size(); i++) {
      str << sql->vdata[i];  // convert string to stream string
      str >> eps[i-1];       // convert to float
    }
  }
  return status;
}
// ============================================================================
// fillSqlNemoRange                                                            
// look for simulation in Sql Database file                                    
bool CSnapshotSimIn::fillSqlNemoRange()
{
  std::string select="select * from nemorange where name='"+filename+"'";
  if (verbose) std::cerr << "select = "<<select <<"\n";
  int status = sql->exe(select);
  if (status) {
    if (verbose) sql->display();
    int offset=0;
    assert(sql->vdata[0]==filename);
    addNemoComponent(offset,sql->vdata[1],"all");
    addNemoComponent(offset,sql->vdata[2],"disk");
    addNemoComponent(offset,sql->vdata[3],"bulge");
    addNemoComponent(offset,sql->vdata[4],"halo");
    addNemoComponent(offset,sql->vdata[5],"halo2");
  }
  return status;
}
//                       - - - - - - - - - - - - - - 
//                             ASCII database        
//                       - - - - - - - - - - - - - - 

// ============================================================================
// opendbFile                                                                  
bool CSnapshotSimIn::openDbFile()
{
  bool status=true;
  fi.open(sim_db_file.c_str(),std::ios::in);
  if (! fi.is_open()) {
    std::cerr << "Unable to open file ["<<filename<<"] for reading, aborting...\n";
    status = false;
  }
  if (status) {
    status = findSim();
    if (status) {
      eps_exist = readEpsFile();
    } else {
      eps_exist = false;
    }
  }
  return status;
}
// ============================================================================
// findSim                                                                     
// look for simulation in Database file                                        
bool CSnapshotSimIn::findSim()
{
  bool status=false;
  bool stop = false;
  while (!stop && ! fi.eof()) {           // while ! eof
    std::string line;
    getline(fi,line); // read on eline
    if ( ! fi.eof()) {
      std::istringstream str(line);  // stream line
      std::string parse;
      // following loop parse each lines previously read   
      //
      int cpt=0;
      while (  str >> parse   &&              // something to read 
	       parse[0] != '#' &&             // not commented out 
	       parse[0] != '!'                // not commented out 
	       ) {
	cpt++;
	if (cpt==1) { // simname
	  simname=parse;
	}
	if (cpt==2) { // sim type
	  //std::istringstream ss(parse);
	  //ss >> simtype;
	  simtype=parse;
	  interface_type = simtype;
	  if (interface_type == "Gadget") interface_index=1;
	  else
	    if (interface_type == "Nemo") interface_index=0;
	    else {
	      std::cerr <<"CSnapshotSimIn::findSim => Unknown interface type....\n";
	    }
	}
	if (cpt==3) { // sim's dirname
	  dirname=parse;
	}
	if (cpt==4) { // sim's basename
	  basename=parse;
	}
      }
      if (simname == filename) { // we found simulation
	stop   = true; // we have a snapshot
	status = true; // so we can stop reading
	std::cerr << "SIM DB:Found simulation ["<<simname<<"] in database !\n";
      }
      if (cpt != 4) {
	std::cerr << "\n\nWarning, bad #strings ["<<cpt<<"] parsed\n"
		  << "during CSnapshotSimIn::findSim()....\n";
      }
    }
    else { // end of file
      stop   = true;
      status = false;
    }
  }
  return status;
}
// ============================================================================
// readEpsFile                                                                 
// Read Eps database file                                                      
bool CSnapshotSimIn::readEpsFile()
{
  bool stop = false;
  bool status=true;

  std::ifstream fi;
  std::string simname;
  fi.open(eps_db_file.c_str(),std::ios::in);
  if (! fi.is_open()) {
    std::cerr << "Warning !!! Unable to open file ["<<filename<<"] for reading...\n";
    status = false;
  }
  if (status) {
    while (!stop && ! fi.eof()) {           // while ! eof
      std::string line;
      getline(fi,line); // read on eline
      if ( ! fi.eof()) {
	std::istringstream str(line);  // stream line
	std::string parse;
	// following loop parse each lines previously read   
	//
	int cpt=0;
	while (  str >> parse    &&             // something to read 
		 parse[0] != '#' &&             // not commented out 
		 parse[0] != '!'                // not commented out 
		 ) {
	  cpt++;
	  if (cpt==1) { // simname
	    simname=parse;
	    if (simname == filename) { // we found simulation
	      //stop   = true; // we have a snapshot
	      status = true; // we have all components	      
	      std::cerr << "EPS:Found simulation ["<<simname<<"] in database !\n";
	    }
	  }
	  if (simname == filename) { // we found simulation
	     std::istringstream ss(parse);
	    if (cpt < MAX_EPS+2) { // 
	      ss >> eps[cpt-2];    // file EPS array
	    }
	  }
	} // while ( str ....
	// 
	if (simname == filename) {
	  stop=true;     // simulation has been found
	  assert(cpt>1); // we must have read at least one eps

	  // copy last eps read to next eps   
	  // it's a trick for NEMO simulations
	  // which have only one eps          
	  for (int i=cpt-1; i<MAX_EPS; i++) {
	    std::cerr << "eps shift i="<<i<<" cpt="<<cpt<<" eps="<<eps[cpt-2]<<"\n";
	    eps[i] = eps[cpt-2]; 
	  }
	}
      } // if !eof ...
      else { // end of file
	stop   = true;
	status = false;
      }
    } // while (!stop ...
  } // if (status....
  if (! status) {
    std::cerr<<"\n\nWARNING, simulation ["<<filename<<"] has no entry in the"
	     <<"EPS datafile ["<<uns::CSnapshotInterfaceIn::eps_db_file<<"]\n\n";
  }
  return status;
}

// ============================================================================
// getEps                                                                      
// return the component according to the component requested                   
// if eps does not exist return -1                                             
float CSnapshotSimIn::getEps(const std::string comp)
{
  float status=-1;
  if (eps_exist) {
    if (comp == "gas"  ) status=eps[0];
    if (comp == "halo" ) status=eps[1];
    if (comp == "disk" ) status=eps[2];
    if (comp == "bulge") status=eps[3];
    if (comp == "stars") status=eps[4];
  }
  //std::cerr << "comp ="<<comp<<" status="<<status<<" eps_exist="<<eps_exist<<"\n";
  return status;
}
// ============================================================================
// fillNemoRange                                                               
// look for simulation in Database file                                        
bool CSnapshotSimIn::fillNemoRange()
{
  bool stop = false;
  bool status=true;

  std::ifstream fi;
  int offset;
  fi.open(nemo_range_file.c_str(),std::ios::in);
  if (! fi.is_open()) {
    std::cerr << "Unable to open file ["<<filename<<"] for reading, aborting...\n";
    status = false;
  }
  if (status) {
    while (!stop && ! fi.eof()) {           // while ! eof
      std::string line;
      getline(fi,line); // read on eline
      if ( ! fi.eof()) {
	std::istringstream str(line);  // stream line
	std::string parse;
	// following loop parse each lines previously read   
	//
	int cpt=0;
	while (  str >> parse    &&             // something to read 
		 parse[0] != '#' &&             // not commented out 
		 parse[0] != '!'                // not commented out 
		 ) {
	  cpt++;
	  if (cpt==1) { // simname
	    simname=parse;
	    if (simname == filename) { // we found simulation
	      stop   = true; // we have a snapshot
	      status = true; // we have all components	      
	      std::cerr << "Found simulation ["<<simname<<"] in database !\n";
	      crv.clear();
	      offset=0;
	    }
	  }
	  if (simname == filename) { // we found simulation
	    if (cpt==2) { // #total
	      addNemoComponent(offset,parse,"all");
	    }
	    if (cpt==3) { // #disk
	      addNemoComponent(offset,parse,"disk");
	    }
	    if (cpt==4) { // #bulge
	      addNemoComponent(offset,parse,"bulge");
	    }
	    if (cpt==5) { // #halo
	      addNemoComponent(offset,parse,"halo");
	    }
	    if (cpt==6) { // #halo2
	      addNemoComponent(offset,parse,"halo2");
	    }
	  }
	} // while ( str ....
      } // if !eof ...
      else { // end of file
	stop   = true;
	status = false;
      }
    } // while (!stop ...
  } // if (status....

  return status;
}
// ============================================================================
// addNemoComponent                                                            
int CSnapshotSimIn::addNemoComponent(int& offset, std::string parse,
				    std::string comp )
{
  int nbody;
  std::istringstream ss(parse);
  ss >> nbody;
  if (nbody) {
    uns::ComponentRange cr;
    cr.setData(offset,nbody-1,comp);
    crv.push_back(cr);
    if (comp!="all") {
      offset+=nbody;
    }
  }
  return offset;
}
// ============================================================================
// isNewFrame                                                                  
bool CSnapshotSimIn::isNewFrame()
{
  bool status=false;
  if (valid) {
    if (simtype=="Gadget") {
      status=buildGadgetFile();
    }
    else {
      if (simtype=="Nemo") {
	status=buildNemoFile();
      }
      else {
	std::cerr <<"\nUnknown simulation type ["<<simtype<<"]\n";
      }
    }
    if (status) {
      interface_type  = snapshot->getInterfaceType();
      interface_index = snapshot->getInterfaceIndex();
    }
  }
  return status;
}
// ============================================================================
// buildNemoFile                                                             
bool CSnapshotSimIn::buildNemoFile()
{
  bool status=false;
  if (nemosim != "") {
    status = true;
  }
  else {
    std::string myfile=dirname+'/'+basename;
    if (snapshot) delete snapshot;
    if (0) { // ASCII database
      if (fillNemoRange()) {
	if (verbose) uns::ComponentRange::list(&crv);
      }
    } else {
      if (fillSqlNemoRange()) {
	if (verbose) uns::ComponentRange::list(&crv);
      }
    }
    // try to open NEMO sim
    PRINT("trying top instantiate  CSnapshotNemo("<<myfile<<") verbose="<<verbose<<"\n";)
      snapshot = new CSnapshotNemoIn(myfile, select_part, select_time,verbose);
    if (snapshot->isValidData()) {
      status=true;
      nemosim=myfile;
    } else {
      status=false;
    }
  }
  return status;
}
// ============================================================================
// buildGadgetFile                                                             
bool CSnapshotSimIn::buildGadgetFile()
{
  bool stop=false,status=false;
  int cpt=1;
  // loop on all the possibility of file
  // dirname+basename+nframe
  // ex : gas001_0 gas001_00 gas001_000
  while (!stop && cpt<5) {
    std::ostringstream ss;
    ss << std::setw(cpt) << std::setfill('0') << nframe;
    std::string myfile = dirname+'/'+basename+'_'+ss.str();
    PRINT("CSnapshotSimIn::buildGadgetFile()  myfile=["<<myfile<<"]\n";)

    if (snapshot) delete snapshot;
    // try to open file
    snapshot = new CSnapshotGadgetIn(myfile, select_part, select_time, verbose);
    if (snapshot->isValidData()) {                // file exist         
      float t;
      bool ok=snapshot->getData("time",&t);
      if (ok && checkRangeTime(t)) {              //  time in range     
	status=true;                              //   valid snap       
	stop=true;                                //   get out loop     
      } else {                                    //  time out of range 
	delete snapshot;                          //   del object       
	snapshot = NULL;                          //   NULL for the next
	nframe++;                                 //   try next frame   
      }
    } 
    else {                                        // file does not exist
      delete snapshot;
      snapshot = NULL;
      cpt++;
    }
  }

  if (status) {
    nframe++;  // next frame index
  }
  return status;
}
// ============================================================================
// getCod
// returns: 
// -2 file exist but can't open
// -1 file does not exist
// 0  time does not exist
// 1  time found 
                                                                   
int CSnapshotSimIn::getCod(const std::string select,
			 const float time, float * tcod,
			 const std::string base, const std::string ext)
{
  int status=-3;  // sim not valid
  if (valid) {
    std::string codfile=dirname+'/'+base+'/'+simname+'.'+select+'.'+ext;
    if (tools::Ctools::isFileExist(codfile)) { // cod file exist
      std::cerr << "cod file = "<< codfile << "\n";
      std::ifstream fi;
      status = 0;
      fi.open(codfile.c_str(),std::ios::in);
      if (! fi.is_open()) {
	std::cerr << "Unable to open file ["<<codfile<<"] for reading...\n";
	status = -2;
      }
      else {
	bool stop=false;
	while (!stop && ! fi.eof()) {           // while ! eof
	  status = 0; // time not match
	  std::string line;
	  getline(fi,line); // read line by line
	  if ( ! fi.eof()) {
	    std::istringstream str(line);  // stream line
	    std::string parse;
	    int cpt=0;
	    // get time
	    str >> parse;                 // read time        
	    std::stringstream str2; // convert to stream
	    str2 << parse;
	    str2 >> tcod[cpt++];              // convert to float 
	    if (tcod[0]-0.00001 < time && tcod[0]+0.00001 > time) {
	      while (  str >> parse    &&             // something to read 
		       parse[0] != '#' &&             // not commented out 
		       parse[0] != '!'                // not commented out 
		       ) {
		assert(cpt < 7);
		std::stringstream str2(parse); // read cod data       
		str2 >> tcod[cpt++];          // store in float array
	      } // while str >> ...
	      assert(cpt==7); // bc cpt+1
	      status=1;  // match cod time      
	      stop=true; //  we can stop so read
	    } // if (tcod[0]-0.00
	  } // !fi.eof
	} // while !stop....
      } // else 
      fi.close(); // close cod file
    } else { // cod file does not exist
      status=-1;
      std::cerr << "cod file = "<< codfile << " does not exist\n";
    }
  }
  return status;
}
}
#endif // NOSQLITE3
