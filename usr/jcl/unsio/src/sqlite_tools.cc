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

#ifndef NOSQLITE3
#include <sqlite_tools.h>
namespace jclt {

// ----------------------------------------------------------------------------
// constructor
CSQLite3::CSQLite3 (std::string tablename): zErrMsg(0), rc(0),db_open(false) 
{
  rc = sqlite3_open(tablename.c_str(), &db);
  if( rc ){
    std::cerr << "Can't open database: " << sqlite3_errmsg(db) << "\n";
    sqlite3_close(db);
  } else {
    db_open=true;
  }
}
// ----------------------------------------------------------------------------
//
CSQLite3::~CSQLite3()
{
  sqlite3_close(db);
  vcol_head.clear();
  vdata.clear();
} 

// ----------------------------------------------------------------------------
//
int CSQLite3::exe(std::string s_exe) 
{
  rc = sqlite3_get_table(
			 db,            // An open database 
			 s_exe.c_str(), // SQL to be executed
			 &result,       // Result written to a char *[] that this points to 
			 &nrow,         // Number of result rows written here 
			 &ncol,         // Number of result columns written here
			 &zErrMsg       // Error msg written here 
			 );

  if(vcol_head.size() > 0) {vcol_head.clear();}
  if(vdata.size()>0) {vdata.clear();}

  if( rc == SQLITE_OK ){
    for(int i=0; i < ncol; ++i)
      vcol_head.push_back(result[i]); // First row heading
    for(int i=0; i < ncol*nrow; ++i)
      vdata.push_back(result[ncol+i]);
  }
  sqlite3_free_table(result);
  return ((rc==SQLITE_OK&&ncol>1)?1:0);
}
// ----------------------------------------------------------------------------
//
void CSQLite3::display()
{
  if(vcol_head.size() > 0 ) {
    copy(vcol_head.begin(),vcol_head.end(),
	 std::ostream_iterator<std::string>(std::cerr,"\t"));
    std::cerr << "\n";
    for (unsigned int j=0; j<vdata.size();) {
      for (unsigned int i=0; i<vcol_head.size() ; i++) {
	std::cerr << vdata[j] << "\t";
	j++;
      }
      std::cerr << "\n";
    }
    //copy(vdata.begin(),vdata.end(),
    //     std::ostream_iterator<std::string>(std::cout,"\n")); 
    //std::cout << std::endl;
  } 
}
};
#endif
