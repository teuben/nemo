/* 
 *
 * layout: ascii interpreter of basic YAPP commands
 *      20-oct-92       original version                pjt
 *	22-jan-95	ansi prototypes			pjt
 *	12-apr-95	fixed potential stack problem, no more ARGS
 *      16-sep-96  0.1  added init= keyword             pjt 
 *       7-apr-01       gcc warnings                    pjt
 *      14-apr-01       added color
 *      14-mar-04       skip blank line
 *      17-sep-05       use MAXL_INELEN, added pl_readlines    PJT
 *
 */

#include <stdinc.h>
#include <yapp.h>
#include <strlib.h>
#include <ctype.h>
#include <layout.h>
#include <extstring.h>
#if HAVE_LIBREADLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif

local struct www {                  /*  List of all basic YAPP commands */
    pl_id id;                       /* YAPP function id; see layout.h */
    char *command;                  /* name of YAPP function command */
    char *args;                     /* type encoding, one char for each par */
}  www[] = {
    { Swap,   "swap",     ""     }, 
    { Xscale, "xscale",   "rr"   },
    { Yscale, "yscale",   "rr"   },
    { Ltype,  "ltype",    "ii"   },
    { Line,   "line",     "rr"   },
    { Move,   "move",     "rr"   },
    { Color,  "color",    "i"    },
    { Point,  "point",    "rr"   },
    { Circle, "circle",   "rrr"  },
    { Cross,  "cross",    "rrr"  },
    { Box,    "box",      "rrr"  },
    { Just,   "just",     "i"    },
    { Text,   "text",     "srrrr"},
    { Flush,  "flush",    ""     },
    { Frame,  "frame",    ""     },
    { Init,   "init",     "srrrr"},
    { Stop,   "stop",     ""     },
    { Help,   "?",        ""     },
    { NOP,    "#",        ""     },
    { NOP,    NULL,       NULL   },
};

local string *split(string);
local void layout_help(void);


/*
 * PL_FREAD: read a text file with layout commands, and return
 *           a pointer to a linked list of plcommand's that can
 *           be executed at later times.
 */

plcommand *pl_fread(string file)
{
    stream fp=stropen(file,"r");
    char line[MAX_LINELEN];
    plcommand *p, *start;

    start = (plcommand *) allocate(sizeof(plcommand));  /* alloc first */
    start->id = End;
    p = start;

    while (fgets(line,MAX_LINELEN,fp)) {        /* read all lines in file */
	dprintf(1,"%s",line);
        if (line[0] == '#') continue;
        if (line[0] == '\n') continue;
        if (pl_lread(line,p)>0) {
            p->next = (plcommand *) allocate(sizeof(plcommand));
            p = p->next;
            p->id = End;
        }
    }
    strclose(fp);
    return start;
}

/* 
 * PL_READLINES: read lines interactively, and execute them
 */

void pl_readlines(void)
{
  plcommand p;
  char line[MAX_LINELEN], *cmd;

  p.next = NULL;
#if HAVE_LIBREADLINE 
  warning("** Experimental readline based READLINE layout **");
  while(1) {
    cmd = readline("LAYOUT>");
    if (cmd == 0) break;
    if (strlen(cmd) == 0) continue;
    add_history(cmd);
    if (pl_lread(cmd,&p)>0) {
      pl_exec(&p);
    }
    free(cmd);
  }
#else
  printf("LAYOUT> "); fflush(stdout);
  warning("** Experimental simple stdin READLINE layout **");
  while (fgets(line,MAX_LINELEN,stdin)) {        /* read all lines from stdin */
    if (strlen(line) == 0) continue;
    if (pl_lread(line,&p)>0) {
      pl_exec(&p);
    }
    printf("LAYOUT> "); fflush(stdout);
  }
#endif
  warning(" *** All done with readlines input ***\n");
}



/*
 * PL_LREAD:  read one line, parse it into a plcommand
 *            returns the number of arguments ('argc'), including
 *            the command itself.
 *          
 */

int pl_lread(string line, plcommand *p)
{
    int i, j, n;
    char *cp;
    string *sp;

    sp = split(line);
    n = xstrlen(sp,sizeof(string))-2;   /* number of cmd arguments */

    p->id = NOP;               /* Default is unknown command (NOP) */
    p->next = NULL;                 /* and make sure no more known */

    for (i=0; www[i].command; i++)      /* find command from table */
        if (streq(www[i].command,sp[0])) {            /* if found: */
            dprintf(1,"pl_parse: %s\n",sp[0]);
            cp = www[i].args;
            if (strlen(cp) != n) {
	        warning("pl_parse: Need %d arguments for %s, got %d",strlen(cp),sp[0],n);
                break;
            }
            p->id = www[i].id;
            for (j=0; *cp; j++, cp++)
                switch(*cp) {
                  case 'r': p->pars[j].r = atof(sp[j+1]);   break;
                  case 'i': p->pars[j].i = atoi(sp[j+1]);   break;
                  case 's': p->pars[j].s = scopy(sp[j+1]);  break;
                  default:  warning("type %c not implemented",*cp);
                }
            break;
        }
    if (!www[i].command)
        warning("Command %s not understood ?\n",sp[0]);
    return n+1;
}

/* 
 * PL_EXEC:    execute a linked list of plcommand's through your
 *              current yapp implementation
 */

void pl_exec(plcommand *p)
{	
    for(;;) {                       /* loop over all commands */
        if (p==NULL) break;         /* See if end of linked list */
        if (p->id == End) break;    /* or 'End' also signifies end of work */
	switch(p->id) {             /* execute proper yapp id code */
            case Swap:   plswap();
			 break;
            case Xscale: plxscale(p->pars[0].r,p->pars[1].r);
			 break;
            case Yscale: plyscale(p->pars[0].r,p->pars[1].r);
			 break;
            case Ltype:  plltype(p->pars[0].i,p->pars[1].i);
                         break;
	    case Line:   plline(p->pars[0].r,p->pars[1].r);
			 break;
	    case Move:   plmove(p->pars[0].r,p->pars[1].r);
			 break;
	    case Color:  plcolor(p->pars[0].i);
			 break;
	    case Point:	 plpoint(p->pars[0].r,p->pars[1].r);
			 break;
            case Circle: plcircle(p->pars[0].r,p->pars[1].r,p->pars[2].r);
			 break;
            case Cross:  plcross(p->pars[0].r,p->pars[1].r,p->pars[2].r);
			 break;
            case Box:    plbox(p->pars[0].r,p->pars[1].r,p->pars[2].r);
			 break;
	    case Just:   pljust(p->pars[0].i);
			 break;
            case Text:   pltext(p->pars[0].s,p->pars[1].r,p->pars[2].r,
                                             p->pars[3].r,p->pars[4].r);
			 break;
	    case Flush:  plflush();
			 break;
            case Frame:  plframe();
			 break;
            case Init:   plinit(p->pars[0].s, p->pars[1].r, p->pars[2].r,
                                              p->pars[3].r, p->pars[4].r);
                         break;
            case Stop:   plstop();
			 break;
	    case Help:   layout_help();
                         break;
            case NOP:    break;
	    default: 	 warning("pl_exec: unknown id code %d",p->id);
		 	 break;
	}
        p = p->next;                /* point to the next plcommand */
    } /* for(;;) */
}

#define MAXWORD 256

#define skipspace(c)    while(isspace(*c)) {c++;}
#define skipword(c)     while(*c && !isspace(*c)) {c++;}
#define skipstring(c)   while(*c && *c!='"') {c++;}

local string *split(string line)
{
    char *cp=line;
    int n=0;
    permanent string sp[MAXWORD+1];

    for(;;) {
        skipspace(cp);
        if (*cp==0) {
            break;
        } else if (*cp=='"') {
            cp++;
            sp[n++] = cp;
            skipstring(cp);
            if (*cp=='"') {
                *cp=0;
                cp++;
            } else {        /* bad end of line: forgot to close " */
                warning("Badly formed string: %s",sp[n-1]);
                break;
            }
        } else {
            sp[n++] = cp;
            skipword(cp);
            if (*cp==0) break;   /* legal end of line */
            *cp=0;
            cp++;
        }
    }
    sp[n] = NULL;
    dprintf(1,"Command %s: argc=%d\n",sp[0],n);
    return sp;
}

local void layout_help(void)
{
  int i;

  printf("One of these days there will be help here, for now , these are the commands\n");
  for (i=0; www[i].command; i++) {
    if (strlen(www[i].args) > 0)
      printf("%-10s %s\n",www[i].command, www[i].args);
    else
      printf("%-10s\n",www[i].command);
  }
}



#ifdef TOOLBOX

#include <getparam.h>

string defv[] = {
    "in=???\n       Input yapp layout file",
    "init=f\n       Add plinit() before calling pl_exec()?",
    "stop=f\n       Add plstop() after calling pl_exec()?",
    "VERSION=0.1\n  16-sep-96 PJT",
    NULL,
};

string usage="YAPP interpreter";

nemo_main()
{
    if (getbparam("init"))
        plinit("***",0.0,20.0,0.0,20.0);
    pl_exec(pl_fread(getparam("in")));
    if (getbparam("stop"))
        plstop();
    
}

#endif
