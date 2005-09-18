/*
 * layout.h:  definitions for the YAPP layout interpreter 
 */

#ifndef _layout_h
#define _layout_h

/* 
 * pl_id:  enumeration of all valid YAPP commands 
 */

typedef enum { End, Swap, Xscale, Yscale, Ltype, Line, Move,    /* Valid Yapp */
	       Color,
               Point, Circle, Cross, Box, Just, Text, Flush,
               Frame, Init, Stop, Help, NOP} pl_id;


/*
 * pl_par:     union to hold any type of parameter we support
 */

typedef union pl_par {      /* hold a function parameters */
    char   c;
    int    i;
    float  f;
    double d;
    real   r;
    char  *s;
} pl_par;

#define MAXPAR  5

/*
 * plcommand:   a structure holding the pl_id of the calling function,
 *              and a list of par's, up to MAXPAR.
 */

typedef struct plcommand {  /* a YAPP command: name and parameters */
	pl_id  id;
        pl_par pars[MAXPAR];
        struct plcommand *next;
} plcommand;


plcommand *pl_fread(string file);
int        pl_lread(string line, plcommand *p);
void       pl_readlines(void);
void       pl_exec (plcommand *p);

#endif
