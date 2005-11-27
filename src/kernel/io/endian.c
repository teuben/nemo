/*
 *  endian:      report endianism of this machine
 */

#include <nemo.h>

static short one = 1;


bool BigEndian(void)
{
  char *c_one = (char *) &one;

  if (c_one[0] == 0 && c_one[1] == 1)
    return TRUE;
  else if (c_one[0] == 1 && c_one[1] == 0)
    return FALSE;
  else 
    error("BigEndian: saw ONE = {%d,%d} : should be (0,1) or (1,0); check 2 ?= sizeof(short)= %d",
	  c_one[0],c_one[1],sizeof(short));
  return TRUE;
}

bool LittleEndian(void)
{
  return !BigEndian();
}


#if defined(TOOLBOX)

string defv[] = {
    "verbose=f\n    verbose display?",
    "VERSION=0.2\n  26-nov-2005 PJT",
    NULL,
};

string usage="display endianism (1=big, like G5,Sparc    0=little, like i386)";

string cvsid="$Id$";

void nemo_main()
{
    bool q = BigEndian();
    bool v = getbparam("verbose");
    if (v)
       printf("%s\n", q ? "BigEndian" : "LittleEndian");
    else
       printf("%d\n", q ? 1 : 0);

    if (v) {
#ifdef WORDS_BIGENDIAN
      printf("autoconf ->  big endian\n");
#else
      printf("autoconf ->  little endian\n");
#endif
    }

}

#endif
