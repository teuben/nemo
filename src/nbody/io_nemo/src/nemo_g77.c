/* -------------------------------------------------------------- *\
|* $Id$
|* 
|* Main program for g77 program
|*
\* -------------------------------------------------------------- */

extern int MAIN__();
int main(int argc, char ** argv)
{
  f_setarg(argc, argv);
  f_setsig();
  f_init();
  MAIN__();
  return 0;    /* For compilers that complain of missing return values; */
}
