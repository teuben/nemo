/*
 *	NEMOPASS:  a  curious silly program that simply
 *      passes input to output without command line processing
 *      or other NEMO error checking.
 *
 *
 *      30-oct-03   an experiment, could also solve it with --
 */

#include <stdio.h>

#define MAXBUF   8192

int main(int argc,char *argv[])
{
  char buf[MAXBUF];
  size_t n;

  for(;;) {
    n = fread(buf,1,MAXBUF,stdin);
    if (n<=0) break;
    n = fwrite(buf,1,n,stdout);
    if (n<0) {
      fprintf(stderr,"### Fatal error: NEMOPASS: Could not write\n");
      exit(1);
    }
  } 
  return 0;
}
