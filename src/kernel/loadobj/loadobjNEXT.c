/***************************************************************/
/* File: loadobj.c                                             */
/* NeXT hacking: Peter Teuben  - 6 nov 91                      */
/*       make loadobjtest MACH=NEXT			       */
/* ----------------------------------------------------------- */
/*     loadobj(pathname)   -- loads object file into memory    */
/*     findfn(fnname)      -- return function pointer          */
/*     mysymbols(progname) -- declare current symbols          */
/***************************************************************/

/* Nemo headers */
#include <stdinc.h>
#include <strlib.h>
#include <filefn.h>
#include <loadobj.h>
/* NeXT headers */
#include <sys/loader.h>
#include <nlist.h>
#include <stab.h>
#include <reloc.h>

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
static struct mach_header header;
static char *stringtable;
static struct nlist *localtable;
static long nsymbols;

#define MAXSECTIONS 10

static char *code[MAXSECTIONS];
static char *textsegment;           /* these are aliases that point */
static char *datasegment;           /* into the *code[] array */
static char *bsssegment;            /* to aid older code and readability ?*/
static int textid, dataid, bssid;   /* and these are the *code id's ?*/
static struct section codesect[MAXSECTIONS];


/***************************************************************/
/* Local function declarations                                 */
/***************************************************************/

static symbol lookup(/* name */);
static symbol enter(/* sym */);
static symbol findentry(/* name, eptr, value */);
static void readstrings(/* symtab */);
static void readsymbols(/* symtab, reloc */);
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
    int cmds, sects;
    long fileptr;
    struct load_command    command;
    struct segment_command segment;
    struct section         section;
    struct symtab_command  symtab;

    infile = pathopen("", pathname, "r");
    if (infile == NULL)
	error("loadobj: can't find input file %s",pathname);
    fread((char *) &header, sizeof header, 1, infile);
    if (header.magic != MH_MAGIC)
	error("loadobj: file %s must be MACH .o format\n",pathname);
    dprintf(1,"loadobj: %s; MACH cpu type=%d subtype=%d\n",
                pathname, header.cputype, header.cpusubtype);
    if (header.filetype != MH_OBJECT)
	error("loadobj: file %s is not a simple .o file",pathname);

    for (cmds=0; cmds<header.ncmds; cmds++) {      /* read all load cmds */
        fread((char *)&command, sizeof command, 1, infile);  /* read */
        fseek(infile, (long) (-sizeof(command)), 1);    /* and seed back */
        if(command.cmd == LC_SEGMENT) {
            dprintf(1,"Reading %d Load Command LC_SEGMENT\n",cmds+1);
            fread((char *)&segment, sizeof segment, 1, infile);
            if (segment.nsects > MAXSECTIONS) 
                error("loadobj: Too many sections (%d/%d)",
                            segment.nsects,MAXSECTIONS);
            textid = dataid = bssid = 0;
            for(sects=0; sects<segment.nsects; sects++) { /*get all sections */
                fread((char *)&section, sizeof section, 1, infile); /* read */
                code[sects] = getmem(section.size);         /* allocate */
                bcopy(&section,&codesect[sects],sizeof section); /* backup */
                if(streq(section.sectname,SECT_TEXT)) {    /* make aliases */
                    textsegment = code[sects];
                    textid = sects+1;
                } else if (streq(section.sectname,SECT_DATA)) {
                    datasegment = code[sects];
                    dataid = sects+1;
                } else if (streq(section.sectname,SECT_BSS)) {
                    bsssegment = code[sects];
                    bssid = sects+1;
                } 
                dprintf(1,"  Read Section %s @ 0x%x (0x%x) 0x%x %d Allocated @0x%x\n",
                     section.sectname,
                     section.addr,section.size,section.offset,section.align,
                     code[sects]);
                if (section.size) {    /* BUG: alignment still wrong */
                  fileptr = ftell(infile);
                  fseek(infile,section.offset, 0);
                  fread(code[sects], 1, section.size, infile);
                  fseek(infile,fileptr,0);
                }
            }
        } else if (command.cmd == LC_SYMTAB) {
            dprintf(1,"Reading %d Load Command LC_SYMTAB\n",cmds+1);
            fread((char *)&symtab, sizeof symtab, 1, infile);
            readstrings(&symtab);
            readsymbols(&symtab,TRUE);
            for (sects=0; sects<segment.nsects; sects++) { /* do prev again */
                if(codesect[sects].nreloc) {
                    fileptr = ftell(infile);
                    fseek(infile,codesect[sects].reloff,0);
                    dprintf(1,"processrelo on section %s\n",
                        codesect[sects].sectname);
                    processrelocation(codesect[sects].nreloc,
                            code[sects]);
                    fseek(infile,fileptr,0);
                }
            }
            textid = dataid = bssid = 0;        /* reset to prevent silly */
        } else {
            dprintf(1,"Skipping unknown %d Load Command %d @ 0x%x\n",
                      cmds+1,command.cmd,ftell(infile));
            fseek(infile, (long) command.cmdsize, 1);
        }
    }
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
    int cmds, sects;
    long cmdptr;
    struct load_command    command;
    struct segment_command segment;
    struct section         section;
    struct symtab_command  symtab;

    infile = pathopen(getenv("PATH"), progname, "r");
    if (infile == NULL)
	error("mysymbols: can't open %s\n", progname);
    fread((char *) &header, sizeof header, 1, infile);
    if (header.magic != MH_MAGIC)
	error("mysymbols: executable file %s not MACH",progname);
    if (header.filetype != MH_EXECUTE)
	error("mysymbols: file %s not an executable",progname);
    dprintf(1,"mysymbols: Opened %s\n",progname);

    cmdptr = ftell(infile);
    for (cmds=0; cmds<header.ncmds; cmds++) {      /* read all load cmds */
        fseek(infile, cmdptr, 0);
        fread((char *)&command, sizeof command, 1, infile);  /* read */
        fseek(infile, (long) (-sizeof(command)), 1);    /* and seed back */
        if(command.cmd == LC_SEGMENT) {
            dprintf(1,"Reading %d Load Command LC_SEGMENT @ 0x%x\n",
                                cmds+1, ftell(infile));
            fread((char *)&segment, sizeof segment, 1, infile);
            for(sects=0; sects<segment.nsects; sects++) { /* get all sections */
                fread((char *)&section, sizeof section, 1, infile);
                dprintf(1,"  Reading Section %s @ 0x%x (0x%x) 0x%x %d\n",
                     section.sectname,
                     section.addr,section.size,section.offset,section.align);
            }
            cmdptr = ftell(infile);
        } else if (command.cmd == LC_SYMTAB) {
            dprintf(1,"Reading %d Load Command LC_SYMTAB @ 0x%x\n",
                            cmds+1,ftell(infile));
            fread((char *)&symtab, sizeof symtab, 1, infile);
            cmdptr = ftell(infile);
            readstrings(&symtab);
            readsymbols(&symtab,FALSE);
        } else {
            dprintf(1,"Skipping unknown load command %d @ 0x%x\n",
                      command.cmd,ftell(infile));
            fseek(infile, (long) command.cmdsize, 1);
            cmdptr = ftell(infile);
        }
    }
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

static void readstrings(symtab)
struct symtab_command *symtab;
{
    extern int errno;

    if (fseek(infile, (long) symtab->stroff, 0) < 0)
    	dprintf(1,"   readstrings: fseek returned errno=%d\n",errno);
    dprintf(1,"   readstrings: stringtable @ 0x%x size=%d\n",
    			 symtab->stroff, symtab->strsize);
    stringtable = (string) getmem(symtab->strsize);
    fread(stringtable, 1, symtab->strsize, infile);
}


/***************************************************************/
/* readsymbols(symtab,reloc);                                  */
/*                                                             */
/*     Reads in all of the symbols and defines the external    */
/* symbols in the symtab tree.  If reloc is TRUE, the value    */
/* of each symbol is relocated relative to the start of the    */
/* text segment, and all the symbols are stored in localtable  */
/* for relocation.  The localtable storage should is released  */
/* at the end of each loadfn.                                  */
/***************************************************************/

static void readsymbols(symtab,reloc)
struct symtab_command *symtab;
bool reloc;
{
    struct nlist entryblk;
    symbol entry;           /* struct nlist *entry */
    int i;

    nsymbols = symtab->nsyms;
    if (nsymbols==0) {
        dprintf(1,"readsymbols: warning: no symbols in object file\n");
        return;  /* ? */
    }

    if (reloc) 
        localtable = (symbol) getmem((int) sizeof(entryblk) * nsymbols);

    fseek(infile, (long) symtab->symoff, 0);    /* go to symbols */
    for (i = 0; i < nsymbols; i++) {        /* and read them all */
        entry = (reloc) ? &localtable[i] : &entryblk;
	fread((char *) entry, sizeof *entry, 1, infile);
	if (entry->n_un.n_strx) {
	    entry->n_un.n_name = stringtable + entry->n_un.n_strx;

            if (entry->n_type&N_STAB) continue;

            dprintf(2,"  SYMBOL %s type %d section %d value=0x%x\n",
                    entry->n_un.n_name, entry->n_type & N_TYPE,
                    entry->n_sect, entry->n_value);
            switch(entry->n_type & (N_EXT|N_TYPE)) {
                case N_UNDF|N_EXT:
                    /* in !reloc these can be skipped 
                       in reloc they need to be looked up
                       from the symbol table, so entering them
                       is not needed
                     */
                    break;
                case N_ABS|N_EXT:
                    /* BAD NEWS: cannot do the !reloc case, 
                       these are in a VM shared library
                       Though perhaps their VM value is OK ??? 
                     */
                    if(reloc)warning("N_ABS and reloc");
#if 0
		    if(!reloc && entry->n_value!=0) 
                        enter(entry);           /* try it anyhow ?? */
                    else
                        warning("skipping N_ABS==0");
#endif
                    break;
                case N_SECT|N_EXT:
                    if (reloc) entry->n_value += (long) textsegment;
                    enter(entry);
                    break;
                case N_INDR|N_EXT:
                    warning("N_INDR");
                    break;
            } /* switch */
        } /* if (entry) */
    } /* i */
}


/***************************************************************/
/* processrelocation(size, segment);                           */
/*                                                             */
/*     Processes the relocation information contained in the   */
/* next size items of the input file and makes the necessary   */
/* adjustments to the memory beginning at the segment pointer. */
/***************************************************************/
/*     THIS FUNCTION IS NOT FREE OF BSD-isms yet               */
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

    size *= sizeof item;        /* note difference between BSD convention */
                                /* now size is in the BSD convention */
    for (i = 0; i < size; i += sizeof item) {
	fread((char *) &item, sizeof item, 1, infile);
		dprintf(2,"r_len= %d\n",item.r_length);
		dprintf(2,"r_ext= %d\n",item.r_extern);
		dprintf(2,"r_sym= %d\n",item.r_symbolnum);
	if (item.r_extern) {
	    dprintf(2,"r_symbolnum = %d\n",item.r_symbolnum);
            name = localtable[item.r_symbolnum].n_un.n_name;
            dprintf(2," Looking up %s\n",name);
            extsym = lookup(name);
            if (extsym == NULL)
		error("loadobj: undefined symbol %s", name);
	    offset = extsym->n_value;
	} else {
            offset = (long) textsegment;
	}
        if (item.r_pcrel) offset -= (long) textsegment;
        switch (item.r_length) {
	    case 0 : segment[item.r_address] += offset; break;
            case 1 : error("loadobj: word relocation not supported");
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
    return (findentry(name, &symbase, (symbol) NULL));
}

static symbol enter(sym)
symbol sym;
{
    dprintf(3,"ENTER: %s\n",sym->n_un.n_name);
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
