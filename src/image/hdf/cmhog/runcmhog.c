/*
 *  runcmhog: preprocessor to run (a series of) CMHOG:
 *  (the original version is in NEMO/src/image/hdf/cmhog)
 *
 *  it expects the 'cmhogin' input file to be present in the local
 *  directory, uses that, but any variables can be overridden via
 *  the commandline.
 *
 *  The first argument must be a (non-existent) directory name, in
 *  which a new 'cmhogin' file will be written, and cmhog will be 
 *  run. The 'cmhog' program is allowed to be in the unix PATH
 *
 *  Note: no provision for same key in multiple namelists yet
 *        solution:  precede keyword with "name/"
 *        no new keywords to be added, must be in replace mode
 *
 *  14-may-2003   also be able to create a directory hierarchy in one call
 *  21-jul-2003   fix namelist/key=val documentation
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>


#define MAXCARDS 256

static char *database[MAXCARDS];    /* full lines */
static char *names[MAXCARDS];       /* individual namelist */
static int ncards = 0;

static char *haskey(char *, char *);

static char *program_version = "runcmhog:  Version 1.3a 19-nov-2019";

static int Qdebug = 0;                  /* -d (debug toggle) */

// forward references

int usage(char *name);
void read_namelist(char *filename);
void patch_namelist(char *keyval);
char *haskey(char *line, char *key);
void strinsert(char *a, char *b, int n);
void goto_rundir(char *name);
void write_namelist(char *name);
void run_program(char *exe);


int main(int argc, char *argv[])
{
    int i;
    char *cp;
    char *namelist = "cmhogin";		/* -n name */
    char *exefile = "cmhog";		/* -e name */
    char *rundir;			/* no flag */

    if (argc < 2) usage(argv[0]);
    i = 1;
    while (*argv[i] == '-') {
      cp = argv[i];
      cp++;
      switch (*cp) {
      case 'h':
	usage(argv[0]); break;
      case 'd':
	i++;
	Qdebug = 1;
	break;
      case 'n':
	i++;  if (i==argc) usage(argv[0]);
	namelist = argv[i++];
	break;
      case 'e' :
	i++;  if (i==argc) usage(argv[0]);
	exefile = argv[i++];
	break;
      default:
	break;
      }
    }
    rundir = argv[i++];

    if (Qdebug)
      fprintf(stderr,"%s -n %s %s\n",exefile,namelist,rundir);

    read_namelist(namelist);
    for (; i<argc; i++)
        patch_namelist(argv[i]);
    goto_rundir(rundir);
    write_namelist("cmhogin");
    run_program(exefile);
    exit(0);        
}

int usage(char *name)
{
  fprintf(stderr,"%s\n", program_version);
  fprintf(stderr,"Usage: %s [-d] [-e exe] [-n namelist] run_directory [[namelist/][key=val]]\n",name);
  exit(0);
}

void read_namelist(char *filename)
{
    FILE *fp;
    char line[256], *cp1, *cp2;

    fp = fopen(filename,"r");
    if (fp==NULL) {
        fprintf(stderr,"File %s could not be opened\n",filename);
        exit(0);
    }
    if (Qdebug) fprintf(stderr,"Opening namelist %s\n",filename);
    while (fgets(line,256,fp)) {
        /* allocate a new line, with extra space for later insertions */
        database[ncards] = (char *) malloc(sizeof(line)+1+256);
        strcpy(database[ncards],line);

        cp1 = &line[2];  /* start of namelist */
        cp2 = strchr(cp1,' ');
        if (cp2==0) {
            fprintf(stderr,"%s: Syntax error on line: %s\n",filename,line);
            exit(0);
        }
        *cp2 = 0;
        names[ncards] = (char *) malloc(strlen(cp1)+1);
        strcpy(names[ncards],cp1);
        ncards++;
    }
}

void patch_namelist(char *keyval)
{
    char key[128], val[128], *cp, *cp1;
    int i, j, k = -1;

    /*  bit unsafe, no error checking */

    strcpy(key,keyval);
    cp = strchr(key,'/');
    if (cp) {                   /* ok, need to search for specific namelist */
        *cp++ = 0;
        for (i=0; i<ncards; i++)
            if (strcmp(key,names[i])==0) {
                k = i;
                break;
            }
        if (k<0) {
            fprintf(stderr,"%s: not a valid namelist\n",key);
            exit(1);
        }
        strcpy(val,cp);     /* silly, but need to copy key=val */
        strcpy(key,val);    /* back into key for the next step */
    } 

    cp = strchr(key,'=');
    if (cp==NULL) {
      fprintf(stderr,"argument \"%s\" : did not find '='\n",key);
      exit(1);
    }
    *cp++ = 0;
    strcpy(val,cp);

    for (i=0; i<ncards; i++) {
        if (k>=0 && k!=i) continue;
        cp = haskey(database[i],key);
        if (cp) {                           /* patch a key */
            if (k<0) {     /* check for duplicated if no namelist suppplied */
                for (j=i+1; j<ncards; j++) {
                    if (haskey(database[j],key)) {
                        fprintf(stderr,"%s=: namelist %s and %s have the same key\n",
                        key, names[i], names[j]);
                        exit(0);
                    }
                }
            }
            /* now insert */
            cp1 = strchr(cp,',');
            if (cp1==0) cp1 = strchr(cp,' ');
            if (cp1==0) {
                fprintf(stderr,"Namelist syntax error on: %s",database[i]);
                exit(0);
            }
            strinsert(cp,val,cp1-cp);
        } 
    }
}

char *haskey(char *line,char *key)
{
    char *cp, *w = line, *s = key;

    for (;;) {
        w = strchr(w,key[0]);
        if (w) {
            cp = w;
            s = key;
            while (cp) {
                cp++;
                s++;
                if (*cp==0) return NULL;
                if (*s==0)
                    return (*cp=='=') ? ++cp : NULL;
                if (*s != *cp) break;
            }
            w++;
        } else
            return NULL;
    }
    return NULL;
}


/*
 * insert a string 'b' into 'a' replacing the first 'n' positions into 'a'
 *      (see als hsh.h in ..../hermes/lib  --  PJT)
 */

void strinsert(char *a, char *b, int n)
{
    int    i, idiff, alen, blen;

    alen = strlen(a);
    blen = strlen(b);
    idiff = blen - n;

    if (idiff < 0)          /* shift a to the left */
        for (i = n; i <= alen; i++)
            a[i + idiff] = a[i];
    else if (idiff > 0)     /* shift a to the right */
        for (i = alen; i >= n; i--)
            a[i + idiff] = a[i];

    for (i = 0; i < blen; i++)
        a[i] = b[i];
}





void goto_rundir(char *name)
{
#if 1
  char *cmd = (char *)malloc(strlen(name) + 30);
  struct stat fbuf;
  int retval;
  extern int errno;

  if (name == 0 || *name == 0) {
    exit(1);
  }
  sprintf(cmd,"mkdir -p %s",name);
  if (system(cmd)) {
    fprintf(stderr,"Problem executing \"%s\"\n",cmd);
    exit(1);    
  }
  sprintf(cmd,"%s/cmhogin",name);
  retval = stat(cmd,&fbuf);
  if (retval == 0) {
    fprintf(stderr,"Directory %s already used\n",cmd);
    exit(1);    
  } else if (errno != ENOENT) {
    fprintf(stderr,"Some problem with/in directory %s [%d]\n",cmd,errno);
    exit(1);    
  }
#else
  if (mkdir(name, 0755)) {
    fprintf(stderr,"Run directory %s already exists\n",name);
    exit(1);
  }
#endif
  if (chdir(name)) exit(2);
}

void write_namelist(char *name)
{
    int i;
    FILE *fp;

    fp = fopen(name,"w");
    if (fp==NULL) {
        fprintf(stderr,"File %s could not be written\n",name);
        exit(1);
    }
    for (i=0; i<ncards; i++)
        fprintf(fp,"%s",database[i]);
    fclose(fp);
}



void run_program(char *exe)
{
#if 1
    system(exe);
#else
    if (execlp(exe,NULL)) {
        fprintf(stderr,"Problem executing %s\n",exe);
        exit(1);
    }
#endif
}
