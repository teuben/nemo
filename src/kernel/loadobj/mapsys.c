/*
 * MAPSYS:  this routine translates a generic symbol name into the
 *           one which can be found in the symbol table of the
 *           host loader. See TESTBED at end
 *
 *  Note: the unicos version has not been tested
 * 
 *	xx-xxx-89	written for bodytrans() 	            PJT
 *	17-sep-90	isolated file and moved in with loadobj/    PJT
 *	25-oct-90	cosmetic				    PJT
 *	25-feb-92	happy gcc2.0			  	    PJT
 *      20-jan-94       solaris 2.x				    pjt
 *	10-jan-95	linux					    pjt
 *      15-dec-95       some cleanup and TESTBED additions          pjt
 * TODO:
 	this routine should use -Dmapsys_bsd etc. that are to be
 	set in e.g.. $NEMOBIN/cc; too many systems do not behave
 	from what you think
 	example: gcc on sunos 4.1.4 really needs -Dsparc -Dbsd
	since the gcc compilers doesn't supply enough to figure
	out the proper action below


	system			name	maps to:
	------			----	-------
	sun3 (sunos 3.x,4.x)	f	_f
	sun4 (sunos 3.x,4.x)	f	_f
	linux (a.out!)		f	_f
	
	sun5 (sunos 5.x)	f	f
	sgi (irix 5.x, 6.x)	f	f
	linux (elf!)		f	f
 	
 */

#include <stdinc.h>

/*
 *  Next define may have to be extended to find more bsd machines
 *  bsd machines prepend an underscore to the symbol
 */

#if !defined(bsd) && !defined(SYSV)
#if defined(sun) || defined(ultrix) || defined(_trace_) || defined(NeXT) || defined(linux)
#define bsd 1
#endif
#endif


void mapsys(char *func)
{
    int i, n;

    n = strlen(func);

#if defined(bsd)
    for (i=n+1; i>0; i--)
        func[i] = func[i-1];        /* shift right, also the '\0' */
    func[0] = '_';                  /* and prepend the underscore */
    n++;
#endif

#if defined(unicos)
    for (i=0; i<n; i++)		    /* all caps, no underscores */
        if (islower(func[i]))
            func[i] = toupper(func[i]);
#endif
}

#ifdef TESTBED

#include <getparam.h>

/*
 * Test function
 */

string defv[] = {
    "name=name\n	Name to text",
    "VERSION=1.0a\n	25-oct-90 PJT",
    NULL,
};

nemo_main()
{
    char *name, newname[80];
    stream fstr;

    name = getparam("name");
    strcpy(newname,name);
    mapsys(newname);
    printf("Function name: \"%s\" On this system it is mapped to: \"%s\"\n",
		name,newname);
    printf("Compiling \"int i_%s\" and \"double d_%s()\" to check with nm:\n",
                name,name);

    fstr = stropen("ms-tmp.c","w!");
    fprintf(fstr,"int i_%s;\n",name);
    fprintf(fstr,"double d_%s() { return(1.2345); }\n",name);
    fprintf(fstr,"static int static_%s;\n",name);
    strclose(fstr);
    system("cc -c ms-tmp.c ; nm ms-tmp.o");
}

#endif


