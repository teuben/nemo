/*
 * CSF: copy structured file.
 *
 *      xx-xxx-87   V1.0 Original version - Josh Barnes
 *       9-dec-89   V1.1 added select= keyword	PJT
 *	 4-apr-90   V1.2 fixed 'select=all' bug	PJT
 *	 9-oct-90   V1.3 select= now item= ; must now be matched by case too
 *			 select= now counts which ordinal ones
 *	 5-mar-91   V1.4 allow type switching
 *	18-jul-91       a    last printf is now dprintf(0,
 *	19-feb-92       b    deleted two unused variables
 *	19-aug-92       c    removed lint complaint about free()
 *      20-feb-94       d    ANSI
 *	26-mar-95   V1.5  fixed counting bug in select= and item=, 1 is now
 *			  the first item in all item= selected (as meant)
 *			a fixed NULL vs. 0 warning
 *      11-dec-09   V1.6  experimenting with half precision
 */

#include <stdinc.h>
#include <getparam.h>
#include <extstring.h>
#include <ctype.h>
#include <filestruct.h>

string defv[] = {
    "in=???\n		Input file",
    "out=???\n		Output file",
    "item=\n		Top level selection items [default: all]",
    "select=\n          Selection numbers (1...) [default: all]",
    "convert=\n		Conversion options {d2f,f2d,i2f,f2i,d2i,i2d,h2d,d2h}",
    "VERSION=1.6\n	11-dec-09 PJT",
    NULL,
};

string usage = "Copy a binary structured file";

#define  MAXSEL 512

extern string *burststring(string, string);

local bool check(string *, string);

void nemo_main(void)
{
    stream instr, outstr, nullstr;
    string item, *items;
    string *tags, *sels, *cvt;
    int    i, nsel, select[MAXSEL], cntrd, cntwr;
    bool all;

    instr = stropen(getparam("in"), "r");
    nullstr = stropen("/dev/null","w+");        /* kludge: force write */
    item = getparam("item");
    sels = burststring(getparam("select"),", ");
    nsel = xstrlen(sels, sizeof(string)) - 1;
    for (i=0; i<nsel; i++) {		/* note: array should be sorted */
        select[i] = atoi(sels[i]);
	dprintf(1,"#selected #%d\n",select[i]);
	if (select[i]<1)
            error("%d: bad number for select=, must be > 0",select[i]);
        if (i>0 && (select[i]<select[i-1]))
            error("%d: bad number for select=, array must be sorted",select[i]);
    }
    cvt = burststring(getparam("convert"),", ");

    if (*item == 0) {
        all = TRUE;
	items = (string *) allocate(sizeof(string));    /* PPAP */
	items[0] = NULL;
    } else {
        all = FALSE;
        items = burststring(item,", ");
        if (items==NULL)
            error("error parsing item=%s\n",item);
        dprintf(1,"selected: ");
        for (i=0; items[i]!=NULL; i++)
            dprintf(1," %s",items[i]);
        dprintf(1,"\n\n");
    }
    outstr = stropen(getparam("out"), "w");
    cntrd = 0;      /* keep track of items read */
    cntwr = 0;      /* and written */
    while ((tags = list_tags(instr)) != NULL) {
        if (check(items, *tags)) {
            cntrd++;
            dprintf(1," %d: %s",cntrd,*tags);
	    if (nsel==0 || (nsel>0 && select[cntwr] == cntrd)) {
               copy_item_cvt(outstr, instr, *tags, cvt);
               cntwr++;
               dprintf(1," copied\n");
	    } else {
	       copy_item(nullstr, instr, *tags);
               dprintf(1," skipped\n");
            }
        } else {
            dprintf(1," *: %s",*tags);        	
	    copy_item(nullstr, instr, *tags);
            dprintf(1," skipped\n");
        }
	free((char *)*tags);
	free((char *)tags);
        if (nsel>0 && cntwr >= nsel)     /* early bailout ? */
            break;
    }
    dprintf(0,"On top-level: %d items  read, %d items written\n",cntrd,cntwr);
    strclose(instr);
    strclose(outstr);
}

/* check if name is in list of items - with minimum match in the items list */

local bool check(
string *items,      /* list if items, user wants to have selected */
string name         /* the set name, as found in file */
) {
    int i;
    char lname[32];

    if (items[0] == NULL) return TRUE;		/* designate all */
    strcpy(lname,name);
#if 0
    for (i=0; lname[i]!=NULL; i++)		/* convert to lower case */
        if (isupper(lname[i]))
            lname[i] = tolower(lname[i]);
#endif

    for (i=0; items[i]!=NULL; i++)
        if (streq(items[i],lname))
            return TRUE;
    return FALSE;
}
