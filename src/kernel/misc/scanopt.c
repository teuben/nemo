/*
 * SCANOPT: scan string of the form "word1,word2,..." for match. Warning:
 * words must be separated by exactly one comma -- no spaces allowed!
 */

#include <stdinc.h>

bool scanopt(string opt, string key)
{
    char *op, *kp;

    op = (char *) opt;				/* start scan of options */
    while (*op != 0) 	{			/* loop over words */
	kp = key;
	while ((*op != ',' ? *op : 0) == *kp) {
						/*   compare with key */
	    if (*kp++ == 0)
		return TRUE;			/*     found keyword */
	    op++;
	}
	while (*op != 0 && *op++ != ',') ;	/*   scan for next word */
    }
    return FALSE;				/* keyword not found */
}

#ifdef TOOLBOX

#include <getparam.h>

string defv[] = {
    "opt=foo,bar,fum\n	String to an option from",
    "key=foo\n		Option, check if it's in the string",
    "verbose=true\n     Verbose output (1/0) or use shell status",
    "VERSION=2.0        30-dec-97 PJT",
    NULL,
};

string usage="TESTBED scanopt";

nemo_main()
{
    string opt, key;
    bool verbose = getbparam("verbose");
    int retval;

    opt = getparam("opt");
    key = getparam("key");
    dprintf(1,"scanopt(\"%s\", \"%s\") returns %s\n",
	   opt, key, scanopt(opt, key) ? "true" : "false");
    retval = scanopt(opt,key);
    if (verbose)
        printf("%d\n",retval);
    else
        stop(retval);
}

#endif
