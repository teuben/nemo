/***************************************************************/
/* File: loadobj.c                                             */
/* Last modified on Wed Aug 20 09:22:36 1986 by roberts        */
/* Last modified on Sat Dec 06          1986 by josh           */
/* Last modified on Wed May 17 03:17:00 1989 by Peter (COFF)   */
/* Last modified on Tue Jul 24 16:17:00 1990 by Peter (MIPS)   */
/*      Nov 28 @ Hat Creek: continued debugging. Current OS    */
/*      Ultrix 4.1 is different; e.g. macro FREAD -> FREADM    */
/*	See PJT's comments in code for some portability	issues */
/*	This code is cloned off the original COFF version: MIPS*/
/*	is not quite COFF, symbols and relocation is somewhat  */
/*	different, by the filehdr and sections remain          */
/*	It needs to be linked with standard System 5 library   */
/*	/usr/lib/libmld.a (-lmld)               	       */
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

#include <stdinc.h>	/* NEMO's stdio-stuff */

#include <a.out.h>	/* also includes things as filehdr.h etc... in SYS5 */
			/* but no in MIPS */
#include <ldfcn.h>      /* general SYS5 */

#include <strlib.h>
#include <filefn.h>
#include <loadobj.h>

/* some versions need this, others not ... */
#define FREADM FREAD

/* next one is for silly 3b1: see a.out.h -> syms.h and nlist.h */
#if u3b 
#undef n_name
#endif

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

static LDFILE *ldptr;             /* replaces infile/header in SYS5 */
static FILHDR   header;
static AOUTHDR  oheader;

static int   stringsize;	/* size of stringtable */
static char *stringtable;       /*  strings */
static struct nlist *localtable;
static long nsymbols;	        /* total number of symbols */


#define MAXSECT 10
static char *ssegment[MAXSECT]; /* location of segments in our alloc memory */
static int isegment[MAXSECT];   /* pointer list for relocation */
static int nsect;               /* actual # segments currently in loadobj */

static char *textsegment;	/* points to text of object file */
#if 0
static char *datasegment;	/*	     data                */
static char *bsssegment;	/*	     bss		 */
#endif

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
static string savename( /* *char */);



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
    SCNHDR shead[MAXSECT];  /* texthead, datahead, bsshead -> array */
    int i, totsiz;

    newname = pathfind("", pathname);
    if (newname == NULL)
        error("mysymbols: can't find %s\n",pathname);

    ldptr = ldopen(newname,NULL);   /* Open the object file for COFF access */
    if (ldptr == NULL)
	error("mysymbols: cant ldopen object file %s\n", newname);

    dprintf(4,"\nloadobj : object file = %s\n",newname);
    dprintf(5,"              type = 0x%x\n",TYPE(ldptr));
    dprintf(5,"            offset = 0x%x\n",OFFSET(ldptr));
    dprintf(5,"         #sections = %d\n",HEADER(ldptr).f_nscns);
    dprintf(5,"   symbolic header = @0x%x (size 0x%x\n",
		HEADER(ldptr).f_symptr, HEADER(ldptr).f_nsyms);
    dprintf(5,"  aout opthdr size = %d\n",HEADER(ldptr).f_opthdr);

    if (OFFSET(ldptr))     /* don't allow archives files for now */
        error("mysymbols: COFF nonzero offset for filehdr: no archives yet.\n");
    if (HEADER(ldptr).f_opthdr) {    /* only check if optional header present */
        if (ldohseek(ldptr) == FAILURE)
            error("mysymbols: cannot ldohseek to Unix System header");
        FREADM((char *) &oheader, sizeof(AOUTHDR), 1, ldptr);
        dprintf(5,"     AOUTHDR :  magic = 0x%x\n",oheader.magic);
        if (oheader.magic != OMAGIC)
            error("Object file %s must be OMAGIC.",pathname);
    } 

    nsect = HEADER(ldptr).f_nscns;       /* #sections in COFF file */
    if (nsect > MAXSECT)
        error("Too many sections (%d), only room for %d",nsect,MAXSECT);
    totsiz = 0;

    for (i=0; i<MAXSECT; i++) isegment[i] = -1;
    
    for (i=0; i<nsect; i++) {       /* read all section headers */
        if (ldshread(ldptr,i+1,&shead[i]) == FAILURE)
            dprintf(0,"   failure reading section header %d\n",i+1);
        dprintf(5,"  section %d = %8s @ 0x%x dat@ 0x%x (0x%x) rel@ 0x%x (%d)\n",
		    i+1,shead[i].s_name,shead[i].s_paddr,
		    shead[i].s_scnptr, shead[i].s_size,
		    shead[i].s_relptr, shead[i].s_nreloc);
        totsiz += shead[i].s_size;
        if (streq(shead[i].s_name,_TEXT))       isegment[R_SN_TEXT] = i;
        else if (streq(shead[i].s_name,_DATA))  isegment[R_SN_DATA] = i;
        else if (streq(shead[i].s_name,_BSS))   isegment[R_SN_BSS] = i;
        else if (streq(shead[i].s_name,_RDATA)) isegment[R_SN_RDATA] = i;
/*        else if (streq(shead[i].s_name,_SDATA)) isegment[R_SN_SDATA] = i;*/
        else warning("Cannot index segment %s for relocation",shead[i].s_name);
    }
    dprintf(5,"     Need to allocate %d of memory for new obj.file\n",totsiz);
    textsegment = getmem(totsiz);   /* allocate space for the object file */

    for (i=0; i<nsect; i++) {
        ssegment[i] = textsegment + shead[i].s_paddr;
        if (shead[i].s_paddr + shead[i].s_size > totsiz)
            error("Assumed object file was contigues in memory, fix it\n");
        if (ldsseek(ldptr,i+1) == FAILURE)
            error("Couldn't seek to section %d\n",i+1);
        FREADM(ssegment[i],1,shead[i].s_size,ldptr);
    }

    readstrings();
    readsymbols(TRUE);

    for (i=0; i<nsect; i++) {		/* for each section */
        if (shead[i].s_nreloc == 0) 
            continue;
        if (ldrseek(ldptr,i+1) == FAILURE)	/* seek to reloction */
            error("Failure to seek to relocation of section %d\n",i+1);
        else
        processrelocation(shead[i].s_nreloc, ssegment[i]);  /* process it */
    }
#if 0
    free((char *) localtable);                      /* free locally used -- */
#endif
    free((char *) stringtable);                     /* -- tables */
    nsect=0;
    ldclose(ldptr);
} /* loadobj */

/***************************************************************/
/* fn = findfn(fnname);                                        */
/* (*fn)();                                                    */
/*                                                             */
/*      The findfn routine looks up fnname in the symbol table */
/* and returns a function pointer by which that function may   */
/* be called.  If the name is not defined, NULL is returned.   */
/*      PJT comments:                                          */ 
/*      Shouldn't findfn() know about the naming convention    */
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

    ldptr = ldopen(newname,NULL);       /* open file for COFF access */
    if (ldptr == NULL)
	error("mysymbols: cant ldopen executable %s\n", newname);

    dprintf(4,"\nmysymbols : exec file = %s\n",newname);
    dprintf(5,"              type = 0x%x\n",TYPE(ldptr));
    dprintf(5,"            offset = 0x%x\n",OFFSET(ldptr));
    dprintf(5,"         #sections = %d\n",HEADER(ldptr).f_nscns);
    dprintf(5,"   symbolic header = @%d size %d\n",
		HEADER(ldptr).f_symptr, HEADER(ldptr).f_nsyms);
    dprintf(5,"  aout opthdr size = %d\n",HEADER(ldptr).f_opthdr);

    if (OFFSET(ldptr) != 0)
        error("mysymbols: COFF archive files not supported\n");
    if (HEADER(ldptr).f_opthdr == 0)
        error("mysymbols: COFF has no Unix System header??\n");
    if (ldohseek(ldptr) == FAILURE)
        error("mysymbols: cannot ldohseek to Unix System header\n");

    FREADM((char *) &oheader, sizeof(AOUTHDR), 1, ldptr);
    dprintf(5,"     AOUTHDR :  magic = 0x%x\n",oheader.magic);
    dprintf(5,"             textsize = %d\n",oheader.tsize);
    dprintf(5,"             datasize = %d\n",oheader.dsize);
    dprintf(5,"              bsssize = %d\n",oheader.bsize);

    if (oheader.magic != OMAGIC)
        dprintf(0,"Warning: executable file is not OMAGIC!\n");

    readstrings();
    readsymbols(FALSE);
    free(stringtable);
} /* mysymbols */



/***************************************************************/
/* readstrings();                                              */
/*                                                             */
/*     Reads in the complete string table from the a.out file. */
/* This storage is freed at the end of any of the loadfn calls */
/*                          (mysymbols or loadobj)             */
/***************************************************************/

static void readstrings()
{
    int size,err;
    HDRR symhdr;    /* the symbolic header; this is main diff with COFF */
    char *cp;

    dprintf(5,"readstrings:  experimental\n");

    err = FSEEK(ldptr, HEADER(ldptr).f_symptr, BEGINNING);
    if (err) error("Cannot seek to symbolic header\n");
    FREADM(&symhdr, 1, sizeof(symhdr), ldptr);
    dprintf(5," Symbolic header: \n");
    dprintf(5,"     local strings @ 0x%x size 0x%x\n",
            symhdr.cbSsOffset, symhdr.issMax);
    dprintf(5,"  external strings @ 0x%x size 0x%x\n",
            symhdr.cbSsExtOffset, symhdr.issExtMax);

    nsymbols = symhdr.isymMax + symhdr.iextMax;
    dprintf(5,"   symbol table claims %d (aux: %d) + %d = %d\n",
            symhdr.isymMax, symhdr.iauxMax, symhdr.iextMax, nsymbols);

    stringsize = symhdr.issMax + symhdr.issExtMax; /* this s.table does NOT */
    stringtable = (string) getmem(stringsize);	/* contain the size !!!      */

    err = FSEEK(ldptr, symhdr.cbSsOffset, BEGINNING); /* seek local strings */
    cp = stringtable;                            /* point to destination */
    FREADM(cp, 1, symhdr.issMax, ldptr);             /* read it */

    err = FSEEK(ldptr, symhdr.cbSsExtOffset, BEGINNING); /* seek ext. strings */
    cp += symhdr.issMax;                         /* increase destination */
    FREADM(cp, 1, symhdr.issExtMax, ldptr);          /* read it */
} /* readstrings */

static void emptystrings()	/* declare an empty stringtable */
{
    stringsize = 0;
    stringtable = NULL;
    dprintf(5,"   emptystrings : assume dummy empty stringable anyhow\n");
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
    SYMR psymbol;			/* MIPS def. of a symbol */
    struct nlist entryblk;
    symbol entry;
    int i;
    static int numaux;
    char *cp, *ldgetname(), *newp, *malloc();

		/* nsymbols was computed in: readstrings() */
    dprintf(5,"\nreadsymbols: %d\n",nsymbols);

    if (reloc) localtable = (symbol) getmem(nsymbols*sizeof(symbol));

    if (ldtbseek(ldptr) == FAILURE)         /* get to symbol table */
	error("Could not seek to symbol table\n");
    numaux = 0;                             /* set to no auxent */

    for (i=0; i<nsymbols ;i++) {
        if (ldtbread(ldptr,i,&psymbol) == FAILURE)
            break;
        cp = ldgetname(ldptr,&psymbol);     /* get name of symbol */
	/* cp = &stringtable[psymbol.iss];  /* other way tgo get symbol name */
        dprintf( (reloc ? 5 : 9) ,  
	    "Reading symbol, iss=%d val=%d st=%d sc=%d idx=%d, name=%s (%s)\n",
            psymbol.iss,psymbol.value,psymbol.st, psymbol.sc,
	    psymbol.index,cp,&stringtable[psymbol.iss]);

	entry = (reloc) ? &localtable[i] : &entryblk ;
        entry->n_name = savename(cp);       /* save the name from ldgetname */
	entry->n_value = psymbol.value; 
	entry->n_type = (psymbol.st << 16) + psymbol.sc;  /* just save it */
                /* st = (n_type | 0xFF00) >> 16  */
                /* sc = (n_type | 0x00FF)        */

        switch (psymbol.sc) {
          case scAbs:
                enter(entry); break;
          case scText:
          case scData:
          case scBss:
                if (reloc) entry->n_value += (long) textsegment;
                enter(entry);
                break;
          default:
                break;
        }
    }  /* for */
} /* readsymbols */



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
    struct reloc item;                      /* SYS5 part */
    string name;
    symbol extsym;
    long offset;
    int i, idx;


    dprintf(5,"processrelocation: nreloc = %d start \@ 0x%x\n",
               nreloc,FTELL(ldptr));

    for (i = 0; i < nreloc; i++) {
	FREADM((char *) &item, sizeof item, 1, ldptr);
	dprintf(5,"   item %d:  @ 0x%x  symidx=%d type %d ext %d\n",
                i,item.r_vaddr, item.r_symndx, item.r_type, item.r_extern);
        if (item.r_extern) {	/* if index in external symbol table */
           name = localtable[item.r_symndx].n_name;
           dprintf(5,"   looking for %s type=0x%x ",name,item.r_type);
           extsym = lookup(name);
           if (extsym == NULL) {
              dprintf(5,"loadobj: undefined symbol %s\n", name);
              offset = 0;
           } else
              offset = extsym->n_value;
        } else {		/* else point to segment */
            idx = isegment[item.r_symndx];
            if(idx<0 || idx>=nsect) {
                warning("invalid segment # %d ?",idx);
		offset = 0;
            } else
                offset = (long) ssegment[idx];
        }

        dprintf(5," offset=0x%x\n",offset);
        switch (item.r_type) {
            case R_ABS:
                break;          /* no relocation needed */
	    case R_REFWORD:
            case R_RELLONG:
                *((long *) &segment[item.r_vaddr]) += offset; break;
            case R_RELBYTE:
                segment[item.r_vaddr] += offset; break;
            case R_REFHI:
            case R_REFLO:
            default:
                warning("loadobj: this type of reloc (%d) not yet supported.\n",
                        item.r_type);
	}  /* switch */
    } /* for */
} /* processrelocation */



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
    if (nsymbols < 100)
        dprintf(7,"Entering symbol : %s\n",sym->n_name);

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
    string sname;

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
	    sname = entry->stdata.n_name;
	    entry->stdata = *value;
	    entry->stdata.n_name = sname;
	}
	return (&entry->stdata);
    }
    return (findentry(name, (cmp<0) ? &entry->before : &entry->after, value));

} /* findentry */

/***************************************************************/
/* Local routine that saves a symbol name in                   */
/* safe memory, as our information in 'nlist' items can only   */
/* have pointers to the name. The names live in the string-    */
/* table which was temporary declared in safe memory in        */
/* readstrings()                                               */
/***************************************************************/

static string savename(name)
char *name;
{
    int len;
    char *cp, *getmem();

    len = strlen(name);
    cp = getmem(len+1);
    strncpy(cp,name,len);
    *(cp+len) = '\0';
    return(cp);
}
