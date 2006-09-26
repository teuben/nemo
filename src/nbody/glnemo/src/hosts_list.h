// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2008                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#ifndef HOSTS_LIST_H
#define HOSTS_LIST_H
static const char* hosts_list[]={ 
          "node1"  ,"node2"  ,"node3"   , "node4"   ,
          "node5"  ,"node6"  ,"node7"   , "node8"   , 
          "paxi"   ,"ouzo"   ,"beaver"  , "batis"   ,
          "teddy"  ,"koala"  ,"amos"    , "anixi"   ,
          "loutre" ,"slinky" ,"pentium" , "dunk"    ,
          "oniro"  ,"panda"  ,"ouki"    , "heron"   ,
          "meli"   ,"ionio"  ,"kerkira" , "pilos"   ,
          "pirgos" ,"raccoon","tripa"   , "zeus"    ,
          "apolon" ,"athena" ,"hra"     , "poseidon",
          "elia"   ,0         };
        
	  class ServerList {
	    public:
	      ServerList():h(hosts_list) {};
	      const char ** h;
	  };
#endif // HOSTS_LIST_H
