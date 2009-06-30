#include <iostream>


typedef union intbite
{
  char           ch4[4] ;
  int            n ;
}t_intbite;

main()
{
char cmd_string[] = "test_maschine_big_...";	

   t_intbite      intbe;

  intbe.n  = 1;

  if ( intbe.ch4[0]  == 1 ) {
      printf("  %-22s : this maschine is little endian \n", cmd_string );
     
    }
  else {
      printf("  %-22s : this machine is big endian \n", cmd_string );
     
    }
  int qq=1;
  char * p;
   p = (char *)&qq;
  for (int i=0; i<4;i++) {
  if ( p[0]  == 1 ) {
      printf("  %-22s : this maschine is little endian \n", cmd_string );
     
    }
  else {
      printf("  %-22s : this machine is big endian \n", cmd_string );
     
    }
    std::cerr << "p["<<i<<"]="<<p[i]<<"\n";
    p++;

  }
}
//
