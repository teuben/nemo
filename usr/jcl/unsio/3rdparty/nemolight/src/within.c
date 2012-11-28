/*
 * WITHIN.C: determine if a floating point number is within specified range,
 * represented as a string of the form "<subrange1>,<subrange2>,..." where
 * each <range> is either a single floating point number, or a pair of
 * numbers seperated by a ":".  To allow for small uncertainties in the
 * values tested, floating-point comparison is done with a specified
 * fuzzyness parameter.
 *	dark-ages	created					JEB
 *	19-oct-90	index() -> strchr()			PJT
 *      20-feb-94       ansi
 *	22-jan-95	ansi prototypes, also all in real	 pjt
 *      10-aug-95       allow the #N syntax for snapshots        pjt  
 *      29-nov-06       real->double for more general usage      pjt
 */

#include <stdinc.h>

local int within_count = 0;

bool within(double val, string range, double fuzz)
{
    char *endptr, *subptr, *sepptr, *colptr;
    double sublow, subhi;
    int count;

    subptr = range;
    if (*subptr++ == '#') {   /* special select Nth (first=1) occurence */
        within_count++;
        count = atoi(subptr);
        return count == within_count;
    }

    endptr = range + strlen(range);		/* point to term. NULL */
    for (subptr = range; subptr != endptr; ) {	/* for each subrange */
        sepptr = strchr(subptr, ',');		/*   pnt to subrange end */
	if (sepptr == NULL)			/*   last subrange listed? */
	    sepptr = endptr;			/*     fix up subend ptr */
	colptr = strchr(subptr, ':');		/*   scan subrange for : */
	if (colptr > sepptr)			/*   in another subrange? */
	    colptr = NULL;			/*     then dont use it */
	sublow = atof(subptr) - fuzz/2.0;	/*   set low end of range */
	if (colptr != NULL)			/*   high end specified? */
	    subhi = atof(colptr+1) + fuzz/2.0;	/*     set high end */
	else
	    subhi = sublow + fuzz;		/*     just use low end */
	if (sublow <= val && val <= subhi)	/*   within subrange? */
	    return (TRUE);
	subptr = sepptr;			/*   advance subrange ptr */
	if (*subptr == ',')			/*   more ranges to do? */
	    subptr++;				/*     move on to next */
    }
    return FALSE;
}

#ifdef TESTBED

#include <getparam.h>

string defv[] = {
    "val=1.0\n			 Value to be tested",
    "range=0.5:0.7,1.0,1.2:1.3\n Ranges to test value in",
    "fuzz=0.001\n		 Fuzziness allowed",
    "VERSION=1.2\n               10-aug-95 PJT",
    NULL,
};

string usage="test within()";

nemo_main()
{
    if (within(getdparam("val"), getparam("range"), getdparam("fuzz")))
        printf("within returns TRUE\n");
    else
        printf("within returns FALSE\n");
}

#endif
