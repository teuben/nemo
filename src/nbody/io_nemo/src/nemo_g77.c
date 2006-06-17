/* =================================================================
|  Copyright Jean-Charles LAMBERT - 2005                            
|  e-mail:   Jean-Charles.Lambert@oamp.fr                           
|  address:  Dynamique des galaxies                                 
|            Laboratoire d'Astrophysique de Marseille               
|            2, place Le Verrier                                    
|            13248 Marseille Cedex 4, France                        
|            CNRS U.M.R 6110                                        
| ==================================================================
|* Main program necessary with g77 compiler                         
+----------------------------------------------------------------  */

void f_setarg(int argc, char ** argv);
void f_setsig();
void f_init();

extern int MAIN__();
int main(int argc, char ** argv)
{
#if G77
  f_setarg(argc, argv);
  f_setsig();
  f_init();
#endif
  MAIN__();
  return 0; /* For compilers that complain of missing return values; */
}
/* ----------------------------------------------------------------
|  End of nemo_g77.c                                               
+---------------------------------------------------------------- */
