/***************************************************************/
/* File: loadobj.c                                             */
/* Last modified on Wed Aug 20 09:22:36 1986 by roberts        */
/* Last modified on Sat Dec 06          1986 by josh           */
/* Last modified on Wed May 17 03:17:00 1989 by Peter (COFF)   */
/*	See PJT's comments in code for some portability	issues */
/*	although this code is purely System 5 (ATT-3b1)	       */
/*	It needs to be linked with standard System 5 library   */
/*	/usr/lib/libld.a				       */
/* ----------------------------------------------------------- */
/*     This package is used to implement dynamic loading of    */
/* functions from object files.  Any references in the         */
/* object file to previously defined objects are allowed,      */
/* but no additional searching or unresolved references        */
/* are allowed.  The public contents of the package are:       */
/*                                                             */
/*     loadobj(pathname)   -- loads object file into memory    */
/*     findfn(fnname)      -- return function pointer          */
/*     mysymbols(progname) -- declare current symbols          */
/***************************************************************/

#include <stdinc.h>
#include <a.out.h>	/* also includes things as filehdr.h etc... in SYS5 */
#include <ldfcn.h>      /* SYS5 */
#include <strlib.h>
#include <filefn.h>
#include <loadobj.h>

/* next one is for silly 3b1: see a.out.h -> syms.h and nlist.h */
#undef n_name

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

#if 0
static stream infile;				/* not used in SYS5 */
static struct exec header;			/* not used in SYS5 */
#endif
static LDFILE *ldptr;             /* replaces infile/header in SYS5 */
static FILHDR   header;                          /* this is in SYS5 */
static AOUTHDR  oheader;                         /* this is in SYS5 */

static char *stringtable;    /* first 4 bytes = length ; then strings */
static struct nlist *localtable;
static long nsymbols;

static char *textsegment;	/* points to text of object file */
static char *datasegment;	/*	     data                */
static char *bsssegment;	/*	     bss		 */

/***************************************************************/
/* Local function declarations                                 */
/***************************************************************/

static symbol lookup(/* name */);
static symbol enter(/* sym */);
static symbol findentry(/* name, eptr, value */);
static void readstrings();
static void emptystrings();
static void readsymbols(/* reloc */);
static void processrelocation(/* size, segment */);
static string savename( /* char[8] */);



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
    char *newname;
    SCNHDR texthead, datahead, bsshead;

    newname = pathfind("", pathname);
    if (newname == NULL)
        error("mysymbols: can't find %s\n",pathname);
    ldptr = ldopen(newname,NULL);
    if (ldptr == NULL)
	error("mysymbols: cant ldopen object file %s\n", newname);
    dprintf(4,"\nloadobj : object file = %s\n",newname);
    dprintf(5,"              type = 0x%x\n",TYPE(ldptr));
    dprintf(5,"            offset = 0x%x\n",OFFSET(ldptr));
    dprintf(5,"         #sections = %d\n",HEADER(ldptr).f_nscns);
    dprintf(5,"         #symbols  = %d\n",HEADER(ldptr).f_nsyms);
    dprintf(5,"  aout opthdr size = %d\n",HEADER(ldptr).f_opthdr);

    if (OFFSET(ldptr) != 0)     /* don't allow archives files for now */
        error("mysymbols: COFF nonzero offset for filehdr??? \n");
    if (HEADER(ldptr).f_opthdr != 0)
        error("mysymbols: COFF has a Unix System header??\n");

    ldnshread (ldptr, ".text", &texthead);
    ldnshread (ldptr, ".data", &datahead);
    ldnshread (ldptr, ".bss",   &bsshead);
    dprintf(5,"   textsegment \@ 0x%x   size=0x%x\n",
					texthead.s_paddr,texthead.s_size);
    dprintf(5,"   datasegment \@ 0x%x   size=0x%x\n",
					datahead.s_paddr,datahead.s_size);
    dprintf(5,"    bsssegment \@ 0x%x   size=0x%x\n", 
					bsshead.s_paddr, bsshead.s_size);
    textsegment = 
       getmem( (int) texthead.s_size + datahead.s_size + bsshead.s_size);
    datasegment = textsegment + texthead.s_size;
    bsssegment = datasegment + datahead.s_size;
    if (ldnsseek(ldptr, ".text") == FAILURE) /* there's always text, right? */
        error("loadobj: cannot ldnsseek to .text\n");
    FREAD(textsegment, 1, texthead.s_size,ldptr);
    if (datahead.s_size > 0) {           /* there may be no data though ... */
        if (ldnsseek(ldptr, ".data") == FAILURE)
            error("loadobj: cannot ldnsseek to .data\n");
        FREAD(datasegment, 1, datahead.s_size,ldptr);
    }
    readstrings();
    readsymbols(TRUE);

    ldnrseek(ldptr,".text");                        /* seek to text part */
    processrelocation(texthead.s_nreloc, textsegment);/* and process it */
    ldnrseek(ldptr,".data");                        /* seek to data part */
    processrelocation(datahead.s_nreloc, datasegment);/* and process it */

    free((char *) localtable);                      /* free locally used -- */
    free((char *) stringtable);                     /* -- tables */
}



/***************************************************************/
/* fn = findfn(fnname);                                        */
/* (*fn)();                                                    */
/*                                                             */
/*      The findfn routine looks up fnname in the symbol table */
/* and returns a function pointer by which that function may   */
/* be called.  If the name is not defined, NULL is returned.   */
/*      PJT comments:                                          */ 
/*      Shouldn't findfn() know about the naminbg convention   */
/*      of symbols of your host system (underscores, case)     */
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
    char *newname;

    newname = pathfind(getenv("PATH"), progname);
    if (newname == NULL)
        error("mysymbols: can't find %s along PATH\n",progname);
    ldptr = ldopen(newname,NULL);
    if (ldptr == NULL)
	error("mysymbols: cant ldopen executable %s\n", newname);
    dprintf(4,"\nmysymbols : exec file = %s\n",newname);
    dprintf(5,"              type = 0x%x\n",TYPE(ldptr));
    dprintf(5,"            offset = 0x%x\n",OFFSET(ldptr));
    dprintf(5,"         #sections = %d\n",HEADER(ldptr).f_nscns);
    dprintf(5,"         #symbols  = %d\n",HEADER(ldptr).f_nsyms);
    dprintf(5,"  aout opthdr size = %d\n",HEADER(ldptr).f_opthdr);

    if (OFFSET(ldptr) != 0)
        error("mysymbols: COFF nonzero offset for filehdr??? \n");
    if (HEADER(ldptr).f_opthdr == 0)
        error("mysymbols: COFF has no Unix System header??\n");
    if (ldohseek(ldptr) == FAILURE)
        error("mysymbols: cannot ldohseek to Unix System header\n");
    FREAD((char *) &oheader, sizeof(AOUTHDR), 1, ldptr);
    dprintf(5,"     AOUTHDR :  magic = 0x%x\n",oheader.magic);
    dprintf(5,"             textsize = %d\n",oheader.tsize);
    dprintf(5,"             datasize = %d\n",oheader.dsize);
    dprintf(5,"              bsssize = %d\n",oheader.bsize);

    readstrings();
    readsymbols(FALSE);
    free(stringtable);
}



/***************************************************************/
/* readstrings();                                              */
/*                                                             */
/*     Reads in the complete string table from the a.out file. */
/* This storage is freed at the end of the loadfn call.        */
/*                          (mysymbols or loadobj)             */
/***************************************************************/

static void readstrings()
{
    int size, err;

    dprintf(5,"readstrings:\n");
    if (SYMESZ != AUXESZ)
        error("readstrings: how can I find strings????? check COFF code\n");
    if (ldtbseek(ldptr) == FAILURE) { /* I'm not sure if we should allow this */
					/* maybe just error() out here */
        dprintf(5,"   warning: (COFF) cannot even seek to symbol section\n");
        emptystrings();
        return;
    }
    dprintf(5,"   seeking strings: symbols      \@ 0x%x\n",FTELL(ldptr));
    err = FSEEK(ldptr, SYMESZ * HEADER(ldptr).f_nsyms, CURRENT);
    dprintf(5,"   skipping symbols, get size \@ 0x%x fseek-> %d\n",
				FTELL(ldptr),err);
    if (err != 0) {			/* expected fseek to work */
        emptystrings();
        return;
    }
    err=FREAD((char *) &size, sizeof size, 1, ldptr);
    dprintf(5,"   size of stringtable = %d  fread-> %d\n",size,err);
    if (err != 1) {			/* expected only 1 item to read */
        emptystrings();         /* assume empty */
        return;
    }
    stringtable = (string) getmem(size);
    dprintf(5,"   stringtable declared \@ 0x%x\n",stringtable);
    err=FSEEK(ldptr, -(sizeof size), CURRENT);      /* go back */
    dprintf(5,"   reading strings \@ 0x%x fseek-> %d\n",FTELL(ldptr),err);
    err=FREAD(stringtable, 1, size, ldptr);
    dprintf(5,"   read stringtable, fread-> %d\n",err);
    if (err != size)
        error("readstrings: read incomplete stringtable, size,err=%d %d\n",
                size,err);
}

static void emptystrings()	/* declare empty stringtable */
{
    stringtable = (string) getmem(5);
    dprintf(5,"   emptystrings : assume dummy empty stringable anyhow\n");
    dprintf(5,"   stringtable declared \@ 0x%x\n",stringtable);
    *(stringtable+0) = 0;     /* PJT: is this always true??? */
    *(stringtable+1) = 0;     /*      I mean 4 bytes ... */
    *(stringtable+2) = 0;
    *(stringtable+3) = 4;     /* hence size = 4  */
    *(stringtable+4) = 0;     /* but really be sure have a null string in case */
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
    struct syment symentry;         /* SYS5 temporary symbol entry */
    struct nlist entryblk;
    symbol entry;
    int i;
    static int numaux;

    nsymbols = HEADER(ldptr).f_nsyms;
    dprintf(5,"\nreadsymbols: nsymbols=%d\n",nsymbols);
    if (reloc) localtable = (symbol) getmem((int) nsymbols*SYMESZ);
    ldtbseek(ldptr);                        /* get to symbol table */
    numaux = 0;                             /* set to no auxent */
    for (i = 0; i < nsymbols; i++) {
        entry = (reloc) ? &localtable[i] : &entryblk;
	FREAD((char *) &symentry, sizeof symentry, 1, ldptr);/* syment/auxent */
        if (numaux>0) {
            numaux = 0;    /* at most one numaux ; so reset instead decrement */
            continue;                         /* no further processing needed */
        }
        if (symentry.n_zeroes == 0) 
            entry->n_name = stringtable + symentry.n_offset;
        else
	    entry->n_name = savename(symentry._n._n_name);
        entry->n_value = symentry.n_value;
        entry->n_scnum = symentry.n_scnum;
        entry->n_type = symentry.n_type;
        entry->n_sclass = symentry.n_sclass;
        entry->n_numaux = symentry.n_numaux;
        dprintf( (reloc) ? 5 : 6 ,
		"   entering symbol %-20s v=0x%x sn=%d ty=%d sc=%d\n",
	            entry->n_name, entry->n_value, entry->n_scnum, 
                                   entry->n_type,  entry->n_sclass);
	numaux = entry->n_numaux;		/* if there symaux ? */
#if 1
	if (reloc && entry->n_scnum==0)		/* terrible hack */
		continue;		/* got to find a way to skip these */
#endif
                                /* THIS SWITCH STATEMENT DOES TOO MUCH */
	switch (entry->n_type) {
	case N_ABS:
	    enter(entry);
	    break;
	default:
	    if (reloc) entry->n_value += (long) textsegment;
	    enter(entry);
	    break;
	}
    }
}



/***************************************************************/
/* processrelocation(nreloc, segment);                         */
/*                                                             */
/*     Processes the relocation information contained in the   */   
/* next nreloc bytes of the input file and makes the necessary */
/* adjustments to the memory beginning at the segment pointer. */
/***************************************************************/

static void processrelocation(nreloc, segment)
unsigned long nreloc;
char *segment;
{
#if 0
    struct relocation_info item;		NOT IN SYS5 */
#endif
    struct reloc item;                      /* SYS5 part */
    string name;
    symbol extsym;
    long offset;
    int i;


    dprintf(5,"processrelocation: nreloc = %d start \@ 0x%x\n",
               nreloc,FTELL(ldptr));
    for (i = 0; i < nreloc; i++) {
	FREAD((char *) &item, sizeof item, 1, ldptr);
        name = localtable[item.r_symndx].n_name;
	dprintf(5,"   looking for %-20s type=%d ",name,item.r_type);

            /* NEXT SECTION ONLY WORKS FOR EXTERN SYMBOLS - CAREFUL RECODE */
        extsym = lookup(name);
        if (extsym == NULL)
	    error("loadobj: undefined symbol %s\n", name);
	offset = extsym->n_value;
	dprintf(5," offset=0x%x\n",offset);
            /* ELSE (i.e. if local = static variable ) */
            /* offset = textsegment; */

	switch (item.r_type) {
            case R_ABS:
                break;          /* no relocation needed */
            case R_RELLONG:
                *((long *) &segment[item.r_vaddr]) += offset; break;
            case R_RELBYTE:
                segment[item.r_vaddr] += offset; break;
            default:
                error("loadobj: this type of reloc (%d) not yet supported.\n",
                        item.r_type);
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
    return (findentry(name, &symbase, (symbol) NULL));
}

static symbol enter(sym)
symbol sym;
{
    return (findentry(sym->n_name, &symbase, sym));
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

    if ((entry = *eptr) == NULL) {      /* if at end a of branch */
	if (value == NULL) return NULL;     /* if looking -- nothing found */
	*eptr = entry = (symtree) getmem(sizeof *entry);  /* if entering */
	entry->stdata.n_name = scopy(value->n_name);    /* allocate and copy */
	entry->before = entry->after = NULL;
	cmp = 0;
    } else {
	cmp = strcmp(name, entry->stdata.n_name);       /* still checking */
    }
    if (cmp == 0) {                 /* found one or initializing */
	if (value != NULL) {   /* entering new symbol or overwriting old one */
	    savename = entry->stdata.n_name;
	    entry->stdata = *value;
	    entry->stdata.n_name = savename;
	}
	return (&entry->stdata);
    }
    return (findentry(name, (cmp<0) ? &entry->before : &entry->after, value));

}

/***************************************************************/
/* Local routine that saves a short (<= 8 char) symbol name in */
/* safe memory, as our information in 'nlist' items can only   */
/* have pointers to the name. Long names live in the string-   */
/* table which was already declared in safe memory in          */
/* readstrings()                                               */
/***************************************************************/

static string savename(shortname)
char shortname[8];
{
    int len;
    string cp;
    
    len = strlen(shortname);
    if (len>8)
        len=8;
    cp = getmem(len+1);
    strncpy(cp,shortname,len);
    *(cp+len) = '\0';
    return(cp);
}


#ifdef TESTBED

string defv[] = {
    "func=(n==0?1:n*f(n-1))",
    "low=1",
    "high=10",
     NULL,
};

double sin(), cos();
static int nmain = 0;
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
    tmpfile = fopen("ld-tmp.c", "w");               /****** TMPFILE ******/
    fprintf(tmpfile, "double sin(), cos();\n");
    fprintf(tmpfile, "extern int nmain;\n");
    fprintf(tmpfile, "static int nlocal=0;\n");
    fprintf(tmpfile, "int f(n)\n");
    fprintf(tmpfile, "int n;\n");
    fprintf(tmpfile, "{\n");
    fprintf(tmpfile, "    printf(\" calling f(n=%%d)\",n);\n");
    fprintf(tmpfile, "    nmain += n;\n");
    fprintf(tmpfile, "    return (%s);\n", fdef);
    fprintf(tmpfile, "}\n");
    fclose(tmpfile);                                /**********************/
    if (system("cc -c ld-tmp.c") != 0)             
	error("function %s does not parse\n", fdef);
    loadobj("ld-tmp.o");
#if 0
    if (system("rm -f ld-tmp.c ld-tmp.o") != 0)
	error("cannot rm ld-tmp.*\n");
#endif
    fn = (iproc) findfn("f");	        /* SYS5 has no leading underscore */
    if (fn == NULL)
	error("function not correctly defined\n");
    for (i = low; i <= high; i++)
	printf("f(%d) = %d,  nmain=%d\n", i, (*fn)(i),nmain);
    exit(0);
}

#endif
