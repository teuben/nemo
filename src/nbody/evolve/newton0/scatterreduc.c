/* scatterreduc.c - */

/*
 *  scatterreduc.c: data reduction of the screen output of multiscatter
 *
 *      Feb 1989  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */

#include  <stdinc.h>

#define  MAX_LINE_LENGTH  1000
#define  MAX_BLOCK_LENGTH  10000

/*-----------------------------------------------------------------------------
 *  main  --  the following new blocks are deemed interesting enough to print:
 *            the ones containing:
 *                                      # previous blocks:  # following blocks:
 *               normal form                   1                  1
 *               starting analytic             2                  0
 *               starting a hierarchical       0                  0
 *               halt_criterion                0                  0
 *               =::=                          0                  0
 *               too much energy drift         0                  1
 *               report of                     0                  1
 *
 *            after each type of block is indicated how many blocks before
 *            and after should be printed; however, no block should be 
 *            printed twice, even if they have more than one independent 
 *            reason for being printed.
 *            
 *            if one or more blocks are not printed, ..... is printed instead.
 *            
 *            NOTE: this is written specifically for the following default
 *                  parameters of multiscatter: 
 *                     dt_minor=10.0
 *                     dt_major=100.0
 *                     de_rel_max=0.001
 *                     nstep_de=100
 *-----------------------------------------------------------------------------
 */
main()
    {
    int  how_many_before;
    int  how_many_after;
    bool  print_next;
    bool  print_new;
    bool  print_old;
    bool  print_anc;
    bool  skipping;
    char  new[ MAX_BLOCK_LENGTH ];
    char  old[ MAX_BLOCK_LENGTH ];
    char  anc[ MAX_BLOCK_LENGTH ];                 /* ancient */
    int  getblock();
    bool  isstartblock();
    bool  unexpected();
    bool  interesting();

    printf("\n");

    getblock(new, MAX_LINE_LENGTH, MAX_BLOCK_LENGTH);
    if (isstartblock(new) == FALSE)
	error_msg("\nmain: false start: first line not starting line\n");
    else
	printblock(new);

    getblock(new, MAX_LINE_LENGTH, MAX_BLOCK_LENGTH);
    printblock(new);
    getblock(new, MAX_LINE_LENGTH, MAX_BLOCK_LENGTH);
    printblock(new);
    getblock(new, MAX_LINE_LENGTH, MAX_BLOCK_LENGTH);
    printblock(new);
    
    getblock(anc, MAX_LINE_LENGTH, MAX_BLOCK_LENGTH);
    getblock(old, MAX_LINE_LENGTH, MAX_BLOCK_LENGTH);
    print_new = TRUE;
    print_old = TRUE;
    print_anc = TRUE;
    skipping = FALSE;

    while (getblock(new, MAX_LINE_LENGTH, MAX_BLOCK_LENGTH) > 0)
	{
	if (unexpected(new) == TRUE)
	    {
	    printblock(new);
	    error_msg("\nmain: last block printed was unexpected\n\n");
	    }
	if (interesting(new, &how_many_before, &how_many_after) == TRUE)
	    {
	    print_new = TRUE;
	    if (how_many_before > 1)
		print_anc = TRUE;
	    if (how_many_before > 0)
		print_old = TRUE;
	    if (how_many_after > 0)
		print_next = TRUE;
	    else
		print_next = FALSE;
	    }
	else
	    print_next = FALSE;

	if (print_anc == TRUE)
	    {
	    printblock(anc);
	    skipping = FALSE;
	    }
	else if (skipping == FALSE)
	    {
	    print_skip_message();
	    skipping = TRUE;
	    }	    

	copyblock(anc, old);
	copyblock(old, new);
	print_anc = print_old;
	print_old = print_new;
	print_new = print_next;
	}

    if (print_anc == TRUE)
	printblock(anc);
    if (print_old == TRUE)
	printblock(old);
    }

/*-----------------------------------------------------------------------------
 *  isstartblock  --  
 *-----------------------------------------------------------------------------
 */
local bool  isstartblock(b)
char *b;
    {
    bool  infirstline();

    return(infirstline(b, "Starting a series"));
    }

/*-----------------------------------------------------------------------------
 *  unexpected  --
 *                  NOTE: very inefficient as it stands; potentially costly.
 *-----------------------------------------------------------------------------
 */
local bool  unexpected(b)
char *b;
    {
    bool  infirstline();

    if (infirstline(b, "t = ") == TRUE)
	return(FALSE);
    else if (infirstline(b, "t_now") == TRUE)
	return(FALSE);
    else if (infirstline(b, "tarting") == TRUE)
	return(FALSE);
    else if (infirstline(b, "normal form") == TRUE)
	return(FALSE);
    else if (infirstline(b, "halt_criterion") == TRUE)
	return(FALSE);
    else if (infirstline(b, "=::=") == TRUE)
	return(FALSE);
    else if (infirstline(b, "the ROOT") == TRUE)
	return(FALSE);
    else if (infirstline(b, "Multiscatter") == TRUE)
	return(FALSE);
    else if (infirstline(b, "warning: no binary output") == TRUE)
	return(FALSE);
    else if (infirstline(b, "nbody = ") == TRUE)
	return(FALSE);
    else if (infirstline(b, "Setting up") == TRUE)
	return(FALSE);
    else if (infirstline(b, "too much energy drift") == TRUE)
	return(FALSE);
    else if (infirstline(b, "report of qualitative") == TRUE)
	return(FALSE);
    else
	return(TRUE);
    }

/*-----------------------------------------------------------------------------
 *  interesting  --  
 *-----------------------------------------------------------------------------
 */
local bool  interesting(b, beforeptr, afterptr)
char *b;
int  *beforeptr;
int  *afterptr;
    {
    if (infirstline(b, "normal form") == TRUE)
	{
	*beforeptr = 1;
	*afterptr = 1;
	return(TRUE);
	}
/**********************************************************************
    else if (infirstline(b, "starting analytic") == TRUE)
	{
	*beforeptr = 2;
	*afterptr = 0;
	return(TRUE);
	}
**********************************************************************/
    else if (infirstline(b, "starting a hierarchical") == TRUE)
	{
	*beforeptr = 0;
	*afterptr = 0;
	return(TRUE);
	}
    else if (infirstline(b, "=::=") == TRUE)
	{
	*beforeptr = 0;
	*afterptr = 0;
	return(TRUE);
	}
    else if (infirstline(b, "too much energy drift") == TRUE)
	{
	*beforeptr = 0;
	*afterptr = 1;
	return(TRUE);
	}
    else if (infirstline(b, "report of qualitative") == TRUE)
	{
	*beforeptr = 0;
	*afterptr = 1;
	return(TRUE);
	}
    else if (infirstline(b, "halt_criterion") == TRUE)
	{
	*beforeptr = 0;
	*afterptr = 0;
	return(TRUE);
	}
    else
	return(FALSE);    
    }

/*-----------------------------------------------------------------------------
 *  print_skip_message  --  
 *-----------------------------------------------------------------------------
 */
local void  print_skip_message()
    {
    printf("  .....\n\n");
    }

/*-----------------------------------------------------------------------------
 *  printblock  --  
 *-----------------------------------------------------------------------------
 */
local void  printblock(b)
    {
    printf("%s\n", b);
    }

/*-----------------------------------------------------------------------------
 *  copyblock  --  
 *-----------------------------------------------------------------------------
 */
local void  copyblock(b1, b2)
char *b1;
char *b2;
    {
    strcpy(b1, b2);
    }

/*-----------------------------------------------------------------------------
 *  getblock  --  
 *-----------------------------------------------------------------------------
 */
local int  getblock(s, line_room, block_room)
char *s;
int  line_room;
int  block_room;
    {
    int  i, j;
    int  getline();

    i = 0;

    if ((j = getline(s, line_room)) <= 1)
	;                             /* silly trick to stomach leading '\n' */
    else
	{
	s += j;
	*s++ = '\n';
	block_room -= j + 1;
	i++;
	}

    while ((j = getline(s, line_room)) > 1 && block_room > line_room)
	{
	s += j;
	*s++ = '\n';
	block_room -= j + 1;
	i++;
	}

    return(i);
    }

/*-----------------------------------------------------------------------------
 *  getline  --  
 *-----------------------------------------------------------------------------
 */
local int  getline(s, room)
char *s;
int  room;
    {
    int  i, c;

    i = 0;
    while (--room > 0 && (c = getchar()) != EOF  && c != '\n')
	s[i++] = c;

    if (c == 'n')
	s[i++] = c;

    s[i] = NULL;

    return(i);
    }

/*-----------------------------------------------------------------------------
 *  infirstline  --  
 *-----------------------------------------------------------------------------
 */
local bool  infirstline(s, sub)
char *s;
char *sub;
    {
    int  sub_length;

    sub_length = strlen(sub);

    while(*s != NULL && *s != '\n')
	{
	if (*s == *sub)
	    if (strncmp(s, sub, sub_length) == 0)
		return(TRUE);
	s++;
	}

    return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *  error_msg  --  
 *-----------------------------------------------------------------------------
 */
local void  error_msg(s)
char *s;
    {
    printf("%s", s);
    exit();
    }

/* endof: scatterreduc.c */
