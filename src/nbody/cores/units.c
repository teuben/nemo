/*
 * UNITS.C: code associated with symbolic names for units.
 */

#include <stdinc.h>
#include <ctype.h>
#include "units.h"

local string unittab[] = _UNIT_TABLE;

/*
 * A2UNIT: function to convert symbolic name for a system of units
 * to its corresponding numerical value, by position in _UNIT_TABLE.
 */

int a2unit(s)
string s;
{
    int i;
    char sup[128];

    strncpy(sup, s, 127);
    for (i = 0; sup[i] != 0; i++)
	if (islower(sup[i]))
	    sup[i] = toupper(sup[i]);
    for (i = 0; unittab[i] != NULL; i++)
	if (streq(sup, unittab[i]))
	    return (i + 1);

    /* Error section: */

    printf("Units are: ");
    for (i = 0; unittab[i] != NULL; i++)
        printf(" %s",unittab[i]);
    printf("\n");
    warning("a2unit: bad unit \"%s\"", s);
    return(-1); /* fool lint and friends */
}

#ifdef TESTBED

#include <getparam.h>

string defv[] = {
    "units=VIRIAL\n	Test",
    "VERSION=1.1\n      13-oct-90 PJT",
    NULL,
};

main(argc, argv)
int argc;
string argv[];
{
    string units;

    initparam(argv, defv);
    units = getparam("units");
    printf("a2unit(\"%s\") = %d\n", units, a2unit(units));
}

#endif
