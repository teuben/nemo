#include <stdio.h>
#include <string.h>

#if defined(INT)
typedef int bool;
#else
typedef short bool;
#endif
typedef char *string;

/* local functions */

static bool matchname(string bind, string name, bool exact);
static string bindpar(string name, string value);

#ifdef PROTO
static bool matchname(string bind, string name, bool exact)
#else
static bool matchname(bind, name, exact)
string bind, name;
bool exact;
#endif
{
    char *bp, *np;
    
    if (name==NULL) return 1;

    bp = bind;
    np = name;
    while (*bp == ' ')          /* skip initial blanks */
        bp++;
    while (*np == ' ')          /* skip initial blanks */
        np++;
    while (*bp == *np) {
        bp++;
        np++;
    }
    if (exact)
    	return (*bp == '=' || *bp == ' ') && *np == NULL;
    else
        return *np==NULL || *np==' ';
}

static string bindpar(name, value)
string name;
string value;
{
    return "yes";
}
