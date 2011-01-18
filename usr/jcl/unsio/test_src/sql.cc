
/*
  Very simple C program.

   Compile:
     gcc -o simplesqlite3cpp2 simplesqlite3cpp2.cc  -Wall -W -O2 -Wl,-R/usr/local/lib -lsqlite3

   Usage:

       ./simplesqlite3cpp2 test.db "create table emp (name varchar(15),age int,weight double)"
       ./simplesqlite3cpp2 test.db "insert into emp (name,age,weight) values ('Bob',47,172)"
       ./simplesqlite3cpp2 test.db "insert into emp (name,age,weight) values ('Sue',38,134)"

       ./simplesqlite3cpp2 test.db "select * from emp"

       Headings                                        
       name    age     weight                          
                                                       
       Data                                            
       Bob     47      172     Sue     38      134     





  
   Note sqlite3 shared library, by default, installs in /usr/local/lib. 
   The compile command above will directly link the full path of 
   this library into this program.

*/

#include "sqlite_tools.h"
#include <iostream>
#include <cstdlib>
int main(int argc, char **argv){
  if( argc!=3 ){
   std::cerr << "Usage: " << argv[0] << " DATABASE SQL-STATEMENT" << std::endl;
   exit(1);
  }
#ifndef NOSQLITE3
  jclt::CSQLite3 sql(argv[1]);
  sql.exe(argv[2]);

  if( sql.vcol_head.size() > 0 )
    {
      std::cout << "Headings " <<  sql.vcol_head.size()<< " " << std::endl;
      copy(sql.vcol_head.begin(),sql.vcol_head.end(),std::ostream_iterator<std::string>(std::cout,"\t")); 
      std::cout << std::endl << std::endl;
      std::cout << "Data " << sql.vdata.size() << " " << std::endl;
      for (unsigned int j=0; j<sql.vdata.size();) {
	for (unsigned int i=0; i<sql.vcol_head.size() ; i++) {
	  std::cerr << sql.vdata[j] << "\t";
	  j++;
	}
	std::cerr << "\n";
      }
      //copy(sql.vdata.begin(),sql.vdata.end(),std::ostream_iterator<std::string>(std::cout,"\n")); 
      //std::cout << std::endl;
    }

  return 0;
#endif
}


