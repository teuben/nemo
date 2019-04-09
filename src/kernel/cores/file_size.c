/* file_size, file_lines:
 *
 *      FILE_SIZE:   get size (in bytes) of a file,
 *                   return -1 if file not found or some other error
 *                   find out by looking at 'errno'
 *
 *      xx-xxx-xx   Original version        PJT
 *      25-nov-90   doc written - nemo_main PJT
 *	25-feb-92   happy gcc 2.0	    PJT
 *	 5-nov-93   added second arg to file_lines!!!	PJT
 *	26-feb-94   ansi and better TESTBED defaults	pjt
 *      20-jun-01   prototypes 
 *	12-sep-01   nemo_
 *      21-sep-03   changed meaning of 'deflen' in nemo_file_lines    PJT
 *       4-apr-19   nemo_file_size::stat() now returns 0 for pipes?   PJT
 */

#include <stdinc.h>

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

int nemo_file_size(char *name)
{
  struct stat buf;
  
  if (stat(name,&buf) == 0) {
    dprintf(9,"nemo_file_size: %d\n",buf.st_size);
    return buf.st_size;
  } else {
    dprintf(1,"nemo_file_size: stat returned errno=%d\n",errno);
    return -1;
  }
}

int nemo_file_time(char *name)
{
  struct stat buf;
  
  if (stat(name,&buf) == 0)
    return (int)buf.st_mtime;
  else {
    dprintf(1,"nemo_file_time: stat returned errno=%d\n",errno);
    return -1;
  }
}

/*
 *      FILE_LINES:   get size (in lines) of a file,
 *                    return -1 if file not found or some other error
 *                    find out by looking at 'errno'
 */

#define BUFSIZE  8192

int nemo_file_lines(char *name, int deflen)
{
    int len, n, cnt=0;
    char *cp, *buf;
    stream str;

    len = nemo_file_size(name);
    if (len<=0) {       /* probably a pipe or special NEMO thing, so we need hints */
      if (deflen == 0) return MAXLINES;
      return ABS(deflen);
    }
    if (len==0) return 0;    /* file has 0 length, so no lines */
    if (deflen > 0) return deflen;
    
#if 0
    if (isatty(fileno(str))) {
        warning("nemo_file_lines: %s is not a true file",name);
        return -1;
    }
#endif

    buf = (char *) allocate(BUFSIZE);
    str = stropen(name,"r");

    /* this is a very expensive operation, preferred it to pass
     * an estimated length to the program 
     */

    while ( (n=fread(buf,1,BUFSIZE,str)) > 0) {
        cp = buf;
        while (n--)
            if (*cp++ == '\n')
                cnt++;
    }
    dprintf(1,"Expensive count=%d on %s\n",cnt,name);
    free(buf);
    strclose(str);
    return cnt;
}


#if defined(TESTBED)

#include <getparam.h>

string defv[] = {
    "in=file_size.c\n       File to test file_size/lines on",
    "nmax=-10000\n          Default size for files on pipes",
    "VERSION=1.5a\n         4-apr-2019 PJT",
    NULL,
};
string usage="file_size.c testbed";

void nemo_main(void)
{
    string fname;
    int deflen;
    char cmd[128];
    extern int errno;

    fname = getparam("in");
    deflen = getiparam("nmax");
    printf("%s file_size() = %d %d\n",
        fname, nemo_file_size(fname), nemo_file_time(fname));
    if (errno != 0)
        perror("perror::Error in file_size");

    printf("%s file_lines() = %d\nwc: ",
        fname, nemo_file_lines(fname,deflen));
    sprintf(cmd,"wc %s",fname);
    system(cmd);
        
}
#endif
