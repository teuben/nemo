/***************************************************************/
/* File: loadobj.c                                             */
/* Last modified on Wed Aug 20 09:22:36 1986 by roberts        */
/* Last modified on Sat Dec 06          1986 by josh           */
/*   This is an ATTEMPT to get it running on Alliant ***********/
/*   has #ifdef'd code for 'sun' and 'alliant'       **********/
/* ----------------------------------------------------------- */
/*     This package is used to implement dynamic loading of    */
/* functions from object files.  Any references in the         */
/* object file to previously defined objects are allowed,      */
/* but no additional searching or unresolved references        */
/* are allowed.  The contents of the package are:              */
/*                                                             */
/*     loadobj(pathname)   -- loads object file into memory    */
/*     findfn(fnname)      -- return function pointer          */
/*     mysymbols(progname) -- declare current symbols          */
/***************************************************************/

#include <stdinc.h>
#include <a.out.h>
#include <strlib.h>
#include <filefn.h>
#include <loadobj.h>

/***************************************************************/
/* Types for binary tree symbol table                          */
/***************************************************************/

typedef struct _ste {
   struct nlist stdata;
   struct _ste *before,*after;
} *symtree;

typedef struct nlist *symbol;

/***************************************************************/
/* Package variables                                           */
/*                                                             */
/*     The use of these variables describes their use to a     */
/* reasonably sufficient level of detail.  The only subtlety   */
/* is that localtable is used as a dynamically-allocated       */
/* array of symbol blocks.                                     */
/***************************************************************/

static symtree symbase;

static stream infile;
static struct exec header;
static char *stringtable;
static struct nlist *localtable;
static long nsymbols;

static char *textsegment;
static char *datasegment;
static char *bsssegment;

/***************************************************************/
/* Local function declarations                                 */
/***************************************************************/

static symbol lookup(/* name */);
static symbol enter(/* sym */);
static symbol findentry(/* name, eptr, value */);
static void readstrings();
static void readsymbols(/* reloc */);
static void processrelocation(/* size, segment */);



/***************************************************************/
/* loadobj(pathname);                                          */
/*                                                             */
/*     Loads the object file indicated by pathname.  If        */
/* successful, nothing is returned.  Errors are handled        */
/* by internal calls on the error routine.  The global         */
/* symbols from this file are added to the runtime database.   */
/***************************************************************/

void loadobj(pathname)
string pathname;
{
    dprintf(1,"loadobj: %s\n",pathname);	
    infile = pathopen("", pathname, "r");
    if (infile == NULL)
	error("loadobj: can't find input file\n");
    fread((char *) &header, sizeof header, 1, infile);
#if defined(sun)    
    if (header.a_magic != OMAGIC)
	error("loadobj: sun file must be .o format\n");
#endif
#if defined(alliant)	
    if (header.a_magic != ACMAGIC)
	error("loadobj: alliant file must be .o format\n");
#endif
    fseek(infile, (long) N_TXTOFF(header), 0);
    textsegment =
      getmem((int) (header.a_text + header.a_data + header.a_bss));
    datasegment = textsegment + header.a_text;
    bsssegment = datasegment + header.a_data;
    fread(textsegment, 1, (int) header.a_text, infile);
    fread(datasegment, 1, (int) header.a_data, infile);
    dprintf(5,"   textsegment \@ 0x%x   size=0x%x\n",
					textsegment, header.a_text);
    dprintf(5,"   datasegment \@ 0x%x   size=0x%x\n",
					datasegment, header.a_data);
    dprintf(5,"    bsssegment \@ 0x%x   size=0x%x\n", 
					 bsssegment, header.a_bss);
    readstrings();
    readsymbols(TRUE);
    fseek(infile, (long) (N_TXTOFF(header)+header.a_text+header.a_data), 0);
    processrelocation(header.a_trsize, textsegment);
    processrelocation(header.a_drsize, datasegment);
    free((char *) stringtable);
    free((char *) localtable);
    fclose(infile);
}



/***************************************************************/
/* fn = findfn(fnname);                                        */
/* (*fn)();                                                    */
/*                                                             */
/*      The findfn routine looks up fnname in the symbol table */
/* and returns a function pointer by which that function may   */
/* be called.  If the name is not defined, NULL is returned.   */
/***************************************************************/

proc findfn(fnname)
string fnname;
{
    symbol sym;

    sym = lookup(fnname);
    return ((sym == NULL) ? NULL : (proc) sym->n_value);
}



/***************************************************************/
/* mysymbols(progname);                                        */
/*                                                             */
/*     Loads the symbols defined in progname, which is usually */
/* argv[0] from the main program.  In case argv[0] is not the  */
/* complete name, the path environment variable is used.       */
/***************************************************************/

void mysymbols(progname)
string progname;
{
    dprintf (1,"mysymbols: %s\n",progname);
    infile = pathopen(getenv("PATH"), progname, "r");
    if (infile == NULL)
	error("mysymbols: can't open %s\n", progname);
    fread((char *) &header, sizeof header, 1, infile);
#if defined(sun)    
    if (header.a_magic != ZMAGIC)
	error("mysymbols: sun file must be executable\n");
#endif
#if defined(alliant)
    if (header.a_magic != AVMAGIC)
	error("mysymbols: alliant file must be executable\n");
#endif	
    readstrings();
    readsymbols(FALSE);
    free(stringtable);
    fclose(infile);
}



/***************************************************************/
/* readstrings();                                              */
/*                                                             */
/*     Reads in the complete string table from the a.out file. */
/* This storage is freed at the end of the loadfn call.        */
/*	potential bug: strings may be stripped off	       */
/***************************************************************/

static void readstrings()
{
    int size;
    extern int errno;
    if (fseek(infile, (long) N_STROFF(header), 0) < 0)
    	dprintf(0,"   readstrings: fseek returned errno=%d\n",errno);
    fread((char *) &size, sizeof size, 1, infile);
    dprintf(1,"   readstrings: stringtable @ 0x%x size=%d\n",
    			 N_STROFF(header), size);
    fseek(infile, (long) N_STROFF(header), 0);
    stringtable = (string) getmem(size);
    fread(stringtable, 1, size, infile);
    dprintf(1,"             copied -> @ 0x%x\n",stringtable);
}



/***************************************************************/
/* readsymbols(reloc);                                         */
/*                                                             */
/*     Reads in all of the symbols and defines the external    */
/* symbols in the symtab tree.  If reloc is TRUE, the value    */
/* of each symbol is relocated relative to the start of the    */
/* text segment, and all the symbols are stored in localtable  */
/* for relocation.  The localtable storage should is released  */
/* at the end of each loadfn.                                  */
/***************************************************************/

static void readsymbols(reloc)
bool reloc;
{
    struct nlist entryblk;
    symbol entry;
    int i;

    nsymbols = header.a_syms / sizeof *entry;
    if (nsymbols==0)
        dprintf(0,"readsymbols: warning: no symbols in object file\n");

    if (reloc) localtable = (symbol) getmem((int) header.a_syms);
    fseek(infile, (long) N_SYMOFF(header), 0);
    for (i = 0; i < nsymbols; i++) {
        entry = (reloc) ? &localtable[i] : &entryblk;
	fread((char *) entry, sizeof *entry, 1, infile);
	if (entry->n_un.n_strx) {
	    entry->n_un.n_name = stringtable + entry->n_un.n_strx;
	    switch (entry->n_type) {
	      case N_ABS|N_EXT:
		enter(entry);
		break;
	      case N_TEXT|N_EXT:
	      case N_BSS|N_EXT:
	      case N_DATA|N_EXT:
		if (reloc) entry->n_value += (long) textsegment;
		enter(entry);
		break;
	    }
        }
    }
}



/***************************************************************/
/* processrelocation(size, segment);                           */
/*                                                             */
/*     Processes the relocation information contained in the   */
/* next size bytes of the input file and makes the necessary   */
/* adjustments to the memory beginning at the segment pointer. */
/***************************************************************/

static void processrelocation(size, segment)
unsigned long size;
char *segment;
{
    struct relocation_info item;
    string name;
    symbol extsym;
    long offset;
    int i;

    for (i = 0; i < size; i += sizeof item) {
	fread((char *) &item, sizeof item, 1, infile);
	if (item.r_extern) {
            name = localtable[item.r_symbolnum].n_un.n_name;
            extsym = lookup(name);
            if (extsym == NULL)
		error("loadobj: undefined symbol %s\n", name);
	    offset = extsym->n_value;
	} else {
            offset = (long) textsegment;
	}
        if (item.r_pcrel) offset -= (long) textsegment;
        switch (item.r_length) {
	    case 0 : segment[item.r_address] += offset; break;
            case 1 : error("loadobj: word relocation not supported\n");
            case 2 : *((long *) &segment[item.r_address]) += offset; break;
	}
    }
}



/***************************************************************/
/* Symbol table functions:                                     */
/*      sym = lookup("str")  -- looks up string value in tree  */
/*      sym = enter(sym)     -- enters copy of sym in tree     */
/***************************************************************/

static symbol lookup(name)
string name;
{
    register symbol tmp;
    

    tmp=findentry(name, &symbase, (symbol) NULL);
    dprintf (9,"    lookup symbol %s @ 0x%x\n",name,tmp->n_value);
    return(tmp);
}

static symbol enter(sym)
symbol sym;
{
    dprintf (9,"    enter symbol %s      @ 0x%x\n",
                    sym->n_un.n_name, sym->n_value);
    return (findentry(sym->n_un.n_name, &symbase, sym));
}



/***************************************************************/
/* Local routine that does the work for lookup and enter.      */
/* This is separated out because this level represents the     */
/* appropriate recursive formulation.                          */
/***************************************************************/

static symbol findentry(name, eptr, value)
string name;
symtree *eptr;
symbol value;
{
    int cmp;
    symtree entry;
    string savename;

    if ((entry = *eptr) == NULL) {
	if (value == NULL) return NULL;
	*eptr = entry = (symtree) getmem(sizeof *entry);
	entry->stdata.n_un.n_name = scopy(value->n_un.n_name);
	entry->before = entry->after = NULL;
	cmp = 0;
    } else {
	cmp = strcmp(name, entry->stdata.n_un.n_name);
    }
    if (cmp == 0) {
	if (value != NULL) {
	    savename = entry->stdata.n_un.n_name;
	    entry->stdata = *value;
	    entry->stdata.n_un.n_name = savename;
	}
	return (&entry->stdata);
    }
    return (findentry(name, (cmp<0) ? &entry->before : &entry->after, value));
}

#ifdef TESTBED

string defv[] = {
    "func=(n==0?1:n*f(n-1))",
    "low=1",
    "high=10",
    NULL,
};

double sin(), cos();

main(argc, argv)
int argc;
string argv[];
{
    string getparam(), fdef;
    int getiparam(), low, high, i;
    stream tmpfile;
    iproc fn;

    (void) (sin(1.0) + cos(1.0));	/* insure sin(), cos() are loaded */
    initparam(argv, defv);
    fdef = getparam("func");
    low = getiparam("low");
    high = getiparam("high");
    mysymbols(getparam("argv0"));
    tmpfile = fopen("ld-tmp.c", "w");
    fprintf(tmpfile, "double sin(), cos();\n");
    fprintf(tmpfile, "int f(n)\n");
    fprintf(tmpfile, "int n;\n");
    fprintf(tmpfile, "{\n");
    fprintf(tmpfile, "    return (%s);\n", fdef);
    fprintf(tmpfile, "}\n");
    fclose(tmpfile);
    if (system("cc -f68881 -c ld-tmp.c") != 0)
	error("function %s does not parse\n", fdef);
    loadobj("ld-tmp.o");
    if (system("rm -f ld-tmp.c ld-tmp.o") != 0)
	error("cannot rm ld-tmp.*\n");
    fn = (iproc) findfn("_f");
    if (fn == NULL)
	error("function not correctly defined");
    for (i = low; i <= high; i++)
	printf("f(%d) = %d\n", i, (*fn)(i));
    exit(0);
}

#endif
