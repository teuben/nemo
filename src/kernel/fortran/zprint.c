

/*
 * BSD Unix fortran callable print-but-do-not-go-to-next-line 
 *      Routine also strips off trailing blanks
 *
 *      example fortran call:
 *
 *		CALL ZPRINT('Hello World')
 *
 *	6-jan-00	converted to F77_FUNC macros
 */

#include <nemo.h>


#define MAXLEN  254

#define zprint F77_FUNC(zprint,ZPRINT)

zprint(char *mesg, int mesg_len)
{
  char line[MAXLEN+2];      /* need extra 2 padding for \r and \0 */
  int n;

  for (n=mesg_len-1; n>=0; n--)     /* work backwards to find real strlen */
    if (line[n] != ' ') break;      /* a non-blank is end of line */

  if (n>MAXLEN) n=MAXLEN;           /* if too long, chop rest of line ... */
  if (n<0) return;                  /* empty line: do nothing, return now */

  strncpy(line,mesg,n+1);
  line[++n] = '\r';
  line[++n] = '\0';
  fputs(line,stdout);
  fflush(stdout);
}

#if defined(TESTBED)
main()
{
  int i;
  char line[100];

  printf("Going to work\n");
  for (i=1; i<100; i++) {
    usleep(100);
    sprintf(line,"Working on plane %d.",i);
    zprint(line,strlen(line));
  }
  printf("\nAll done\n");
}
#endif
