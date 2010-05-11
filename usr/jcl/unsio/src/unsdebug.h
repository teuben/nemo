// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008                                       
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

#ifndef UNSDEBUG_H
#define UNSDEBUG_H

#if DEBUG
#define PRINT(A) std::cerr << A
#else
#define PRINT(A) ;
#endif

#endif
