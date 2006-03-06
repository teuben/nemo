/* get fortran I/O header size (4 or 8) */
/* this code is used by configure during NEMO's install */

#include <stdio.h>

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

int file_size(char *name)
{
  struct stat buf;
  
  if (stat(name,&buf) == 0)
    return buf.st_size;
  else {
    fprintf(stderr,"file_size: stat returned errno=%d for %s\n",errno,name);
    return -1;
  }
}

int main(int argc, char *argv[])
{
    int s = file_size(argv[1]);
    
    if (s == 12)
        printf("4\n");
    else if (s == 20)
        printf("8\n");
    else
        printf("%d\n",s);
}

