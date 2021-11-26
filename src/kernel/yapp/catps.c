/*
 * CATPS.C: concat a series of PostScript files produced by yapp_ps.
 *
 *	25-jun-92 fixed error() declaration
 *	22-jun-97 error -> myerror
 */

#include <stdinc.h>

#define STARTPAGE  "%%Page: 1 1\n"
#define SHOWPAGE   "showpage\n"

#define LBUF 256

void myerror(string msg,string p1);
bool readline(string buf, stream str);

int main(int argc, string argv[])
{
    int i;
    stream s;
    bool inprolog;
    char buf[LBUF];

    if (argc < 2)
	myerror("Usage: %s file1 file2 ...\n", argv[0]);
    for (i = 1; i < argc; i++) {
	s = fopen(argv[i], "r");
	if (s == NULL)
	    myerror("### Fatal error: File %s not found\n", argv[i]);
	inprolog = TRUE;
	while (readline(buf, s))		/* loop reading lines */
	    if (inprolog) {			/*   reading file header? */
		if (i == 1)			/*     first input file? */
		    printf("%s", buf);		/*       copy to output */
		if (streq(buf, STARTPAGE))	/*     end of header? */
		    inprolog = FALSE;		/*       new input state */
	    } else {
		if (streq(buf, SHOWPAGE) && i < argc-1) {
		    				/*     end of frame? */
		    (void) readline(buf, s);	/*       get restore cmmd */
		    printf("%s", buf);		/*       and copy it out */
		    break;			/*       exit while loop */
		}
		printf("%s", buf);		/*     copy to output */
	    }
	fclose(s);
    }
	return 0;
}

bool readline(string buf, stream str)
{
    return (fgets(buf, LBUF, str) != NULL);
}

/*
 * Cannot use the NEMO error() call, since it assumes
 *	the user interface initparam()
 */
void myerror(string msg,string p1)
{
    fprintf(stderr,msg,p1);
    exit(-1);
}


