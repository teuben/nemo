/*
 *  story.c: operations on stories and their lines
 *.............................................................................
 *    version 0:  July 1990   Piet Hut               email: piet@iassns.bitnet
 *                            Institute for Advanced Study, Princeton, NJ, USA
 *  9-oct-90    NEMO interface test                                         PJT
 *.............................................................................
 *  non-local functions: 
 **  operating on whole stories:
 *      cpstory, mkstory, rmstory
 **  modifying parts of stories:
 *      addchapter, addtext, checktitle, clean_chapter, givetitle, 
 *      rmbinding, take_out_last_chapter
 **  enquiring about stories:
 *      find_first_chapter, find_next_chapter, findchapter, is_empty_story
 **  I/O on stories:
 *      fgetstory, fputstory, getstory, putstory
 **  operating on (physical) quantities:
 *      addiq, addrq, addbq, addsq, _addivq, _addrvq
 *      findiq, findrq, findbq, findsq, _findivq, _findrvq, rmq
 *.............................................................................
 *     This file contains all functions which directly operate on stories.
 *  Although rather long, the file is not broken up into smaller files in order
 *  to keep the communication between all functions local.  A division in
 *  smaller files would require the communication between some individual
 *  functions to be external, and thereby accessible to user programs,
 *  introducing the possibility of accidental name conflicts and errors.  
 *  This problem is a consequence of the lack of hierarchical shielding of 
 *  environments in C.
 *     Organization: The globally accessible functions are first given in the
 *  order and categorizes as listed above. Then all local functions are given
 *  in alphabetic order.
 *.............................................................................
 *  see also: story.h, snode.h, storyio.h, storyq.h
 *.............................................................................
 */

#include "snode.h"      /* includes "story.h", which includes "stdinc.h" */

/*-----------------------------------------------------------------------------
 *  Local definitions of functions used in story.c
 *-----------------------------------------------------------------------------
 */
local void addstring(), clean_story(), extend_text(), fputchapter(),
           fputsnode(), give_text(), nl_detector(), pusnode(), rm_line_piece(),
           rmsnode(), write_iq(), write_rq(), write_bq(), write_sq(),
           write_ivq(), write_rvq();


/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*                                                                           */
/*                  Operations on whole stories:                             */
/*                                                                           */
/*      cpstory, mkstory, rmstory                                            */
/*                                                                           */
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

/*-----------------------------------------------------------------------------
 *  cpstory  --  ....
 *          note:
 *               the story tree structure is copied, but the same old line
 *               text strings are used. Beware!!!
 *-----------------------------------------------------------------------------
 */
snodeptr cpstory(old_story)
snodeptr old_story;
    {
    snodeptr  new_story;
    snodeptr  old_chapter;
    snodeptr  new_chapter;
    snodeptr  mksnode();

    if (old_story == NULL)
        return(NULL);

    new_story = mksnode();

    Text(new_story) = Text(old_story);

    if (isstory(old_story))
        {
        First_chapter(new_story) = cpstory(First_chapter(old_story));
	old_chapter = First_chapter(old_story);
	new_chapter = First_chapter(new_story);

	while ((old_chapter = Next_chapter(old_chapter)) != NULL)
	    {
	    Next_chapter(new_chapter) = cpstory(old_chapter);
	    new_chapter = Next_chapter(new_chapter);
	    }

	Last_chapter(new_story) = new_chapter;
	}
    
    return(new_story);
    }

/*-----------------------------------------------------------------------------
 *  mkstory  --  ....
 *               ....
 *          note:
 *               at least one empty line, to show that the title in the story
 *               is not just a line text.
 *-----------------------------------------------------------------------------
 */
snodeptr mkstory()
    {
    snodeptr new_story;
    snodeptr empty_chapter;
    snodeptr mksnode();

    new_story = mksnode();
    empty_chapter = mksnode();

    First_chapter(new_story) = Last_chapter(new_story) = empty_chapter;

    return(new_story);
    }

/*-----------------------------------------------------------------------------
 *  rmstory  --  ....
 *               ....
 *-----------------------------------------------------------------------------
 */
void  rmstory(old_story)
snodeptr  old_story;
    {
    rmsnode(old_story);
    }

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*                                                                           */
/*                  Modifying parts of stories:                              */
/*                                                                           */
/*      addchapter, addtext, checktitle, clean_chapter, givetitle,           */
/*      rmbinding, take_out_last_chapter                                     */
/*                                                                           */
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

/*-----------------------------------------------------------------------------
 *  addchapter  --  extends a story by adding a chapter at the end.
 *-----------------------------------------------------------------------------
 */
void  addchapter(old_story, new_chapter)
snodeptr  old_story;
snodeptr  new_chapter;
    {
    if (new_chapter == NULL)
        return;

    if (is_empty_story(old_story))
        {
	rmsnode(Last_chapter(old_story));                   /* rm dummy node */
	First_chapter(old_story) = Last_chapter(old_story) = new_chapter;
	}
    else
        {
        Next_chapter(Last_chapter(old_story)) = new_chapter;
        Last_chapter(old_story) = new_chapter;	
	}
    }

/*-----------------------------------------------------------------------------
 *  addtext  --  extends a story by adding text at the end of the last line of
 *               the story, if the last line does not yet end on a '\n', or
 *               otherwise to the beginning of a newly created last line of the
 *               story. If the text string contains more than one '\n', extra
 *               lines are added to the story, and the string is divided over
 *               the extra lines in such a way that every line contains at most
 *               one '\n', and only at the end, just before the terminating
 *               NULL of the string.
 *          note:
 *               The new string can have no, one, or more new-line characters;
 *               The old story can have a last chapter which is full, partly
 *               filled, or NULL.
 *-----------------------------------------------------------------------------
 */
void  addtext(old_story, new_string)
snodeptr  old_story;
char *new_string;
    {
    int  new_text_length;
    bool  contains_nl;                /* not used here */
    bool  past_nl;

    nl_detector(new_string, &new_text_length, &contains_nl, &past_nl);

    addstring(old_story, new_string, new_text_length);

    if (past_nl)
        addtext(old_story, new_string + new_text_length);
    }

/*-----------------------------------------------------------------------------
 *  checktitle  --  returns TRUE if the  old_story  has  new_string  as title;
 *                  or FALSE otherwise.
 *                     Typical applications: in I/O operations on stories, such
 *                  as in the function  fgetgnode()  in  gnode.c .
 *-----------------------------------------------------------------------------
 */
bool  checktitle(old_story, new_string)
snodeptr  old_story;
char *new_string;
    {
    if (streq(Text(old_story), new_string))
        return(TRUE);
    else
        return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *  clean_chapter  --  locate a chapter in a story which was already
 *                              read in, and remove the contents of that 
 *                              chapter while retaining its title.
 *                         note:
 *                              if more than one chapter is present with the
 *                              indicated title, then all these chapters are
 *                              cleaned.
 *                       beware:
 *                              the memory locations for the old contents are
 *                              freed.
 *-----------------------------------------------------------------------------
 */
bool  clean_chapter(a_story, title)
snodeptr  a_story;
char *title;
    {
    snodeptr  a_chapter;
    bool  found_such_chapter;

    found_such_chapter = FALSE;

    a_chapter = First_chapter(a_story);
    
    while (a_chapter != NULL)
        {
	if (streq(Text(a_chapter), title))
	    {
	    clean_story(a_chapter);
	    found_such_chapter = TRUE;
	    }
	a_chapter = Next_chapter(a_chapter);
	}

    return(found_such_chapter);
    }

/*-----------------------------------------------------------------------------
 *  givetitle  --  gives a title to a story. No '\n' characters are allowed 
 *                 in the title of a story, since that would not fit within 
 *                 the one-line story delimiter.
 *-----------------------------------------------------------------------------
 */
void  givetitle(old_story, new_string)
snodeptr  old_story;
char *new_string;
    {
    int  n;
    bool  contains_nl;
    bool  past_nl;                /* not used here */

    if (new_string == NULL)
        {
	if (Text(old_story) != NULL)
	    free(Text(old_story));
	Text(old_story) = NULL;
	}
    else
	{
        nl_detector(new_string, &n, &contains_nl, &past_nl);

        if (contains_nl)
            error("givetitle: new_string contains a new_line character");

        give_text(old_story, new_string, n);
	}
    }

/*-----------------------------------------------------------------------------
 *  rmbinding  --  remove the old_story, freeing up the links between its
 *                 chapters. The chapters remain intact.
 *            note:
 *                 It is the responsibility of the user to make sure that the
 *                 chapters are still pointed at by some other means, otherwise
 *                 these chapters can neither be read nor garbage collected!
 *-----------------------------------------------------------------------------
 */
void  rmbinding(old_story)
snodeptr  old_story;
    {
    snodeptr  a_chapter;
    snodeptr  new_chapter;

    a_chapter = First_chapter(old_story);
    while (a_chapter)
        {
	new_chapter = Next_chapter(a_chapter);
	Next_chapter(a_chapter) = NULL;
	a_chapter = new_chapter;
	}

    if (Text(old_story))
        free (Text(old_story));

    free(old_story);
    }

/*-----------------------------------------------------------------------------
 *  take_out_last_chapter  --  Return a pointer to the last chapter of a story,
 *                             after taking this chapter out of that story,
 *                             which now contains one less chapter.
 *-----------------------------------------------------------------------------
 */
snodeptr  take_out_last_chapter(st)
snodeptr  st;
    {
    snodeptr  ac;                   /* a chapter            */
    snodeptr  lc;                   /* last chapter         */
    snodeptr  nc;                   /* next-to-last chapter */

    if ((nc = First_chapter(st)) == NULL)
        error("take_out_last_chapter: story presented has no chapters");

    if ((lc = Next_chapter(nc)) == NULL)
        {
	First_chapter(st) = Last_chapter(st) = NULL;
	return(nc);
	}

    while (ac = Next_chapter(lc))
        {
	nc = lc;
	lc = ac;
	}

    Last_chapter(st) = nc;
    Next_chapter(nc) = NULL;

    return(lc);
    }

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*                                                                           */
/*                      Enquiring about stories:                             */
/*                                                                           */
/*      find_first_chapter, find_next_chapter, findchapter, is_empty_story   */
/*                                                                           */
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

/*-----------------------------------------------------------------------------
 *  find_first_chapter  --  an extra abstraction barrier: the real macros used
 *                          are defined in  snode.h , while the user
 *                          programs only include  story.h , which does not
 *                          contain the macro First_chapter().
 *                             Typical applications: in transformations between
 *                          intermediate representations (as stories) and
 *                          internal representations (as specific C data
 *                          structures); an example can be found in the
 *                          function  install_gnodes()  in the file  gnode.c .
 *-----------------------------------------------------------------------------
 */
snodeptr  find_first_chapter(a_story)
snodeptr  a_story;
    {
    return(First_chapter(a_story));
    }

/*-----------------------------------------------------------------------------
 *  find_next_chapter  --  an extra abstraction barrier: the real macros used
 *                         are defined in  snode.h , while the user
 *                         programs only include  story.h , which does not
 *                         contain the macro Next_chapter().
 *                             Typical applications: in transformations between
 *                          intermediate representations (as stories) and
 *                          internal representations (as specific C data
 *                          structures); an example can be found in the
 *                          function  install_gnodes()  in the file  gnode.c .
 *-----------------------------------------------------------------------------
 */
snodeptr  find_next_chapter(a_story)
snodeptr  a_story;
    {
    return(Next_chapter(a_story));
    }

/*-----------------------------------------------------------------------------
 *  findchapter  --  locate a chapter in a story which was already read in.
 *              note:
 *                   if more than one chapter is present with the indicated
 *                   title, then the first of these chapters is returned.
 *-----------------------------------------------------------------------------
 */
snodeptr  findchapter(a_story, title)
snodeptr  a_story;
char *title;
    {
    snodeptr  a_chapter;

    a_chapter = First_chapter(a_story);
    
    while (a_chapter != NULL)
	if (streq(Text(a_chapter), title))
	    return(a_chapter);
	else
	    a_chapter = Next_chapter(a_chapter);

    return(NULL);                       /* chapter not found */
    }

/*-----------------------------------------------------------------------------
 *  is_empty_story  --  A story is truely an empty story, if:
 *                      1) it has only one chapter, and:
 *                      2) this single chapter is itself not a story, and:
 *                      3) this single chapter contains no text
 *-----------------------------------------------------------------------------
 */
bool  is_empty_story(old_story)
snodeptr  old_story;
    {
    snodeptr  lc;

    if (! isstory(old_story))
        error("is_empty_story: snode presented is not a story");

    lc = Last_chapter(old_story);

    if (lc != First_chapter(old_story))
        return(FALSE);
    else if (isstory(lc))
        return(FALSE);
    else if (Text(lc) != NULL)
        return(FALSE);
    else
        return(TRUE);
    }

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*                                                                           */
/*                  Input/Output operations on stories:                      */
/*                                                                           */
/*      fgetstory, fputstory, getstory, putstory                             */
/*                                                                           */
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

#include "storyio.h"

/*-----------------------------------------------------------------------------
 *  fgetstory  --  read a whole story from a specified file
 *           NOTE:
 *                Current implementation simple, but not efficient: much
 *                character copying occurs.
 *-----------------------------------------------------------------------------
 */
snodeptr  fgetstory(fp)
FILE  *fp;
    {
    snodeptr  first_line;
    snodeptr  fgetline();
    snodeptr  fgetchapter();
    
    if ((first_line = fgetline(fp)) == NULL)
        return(NULL);

    if (! is_story_begin_line(first_line))
        {
	warning("fgetstory: first line not story_begin_line; returning NULL");
        return(NULL);
	}
    else
	return(fgetchapter(fp, first_line));
    }

/*-----------------------------------------------------------------------------
 *  fputstory  --  write a whole story onto a specified file
 *-----------------------------------------------------------------------------
 */
void  fputstory(fp, a_story)
FILE  *fp;
snodeptr  a_story;
    {
    fputchapter(fp, a_story);
    }

/*-----------------------------------------------------------------------------
 *  getstory  --  read a whole story from the standard input
 *-----------------------------------------------------------------------------
 */
snodeptr  getstory()
    {
    return( fgetstory(stdin) );
    }

/*-----------------------------------------------------------------------------
 *  putstory  --  write a whole story onto the standard output
 *-----------------------------------------------------------------------------
 */
void  putstory(a_story)
snodeptr  a_story;
    {
    fputsnode(stdout, a_story);
    }

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*                                                                           */
/*            Operations on (physical) quantities in stories:                */
/*                                                                           */
/*      addiq, addrq, addbq, addsq, _addivq, _addrvq,                        */
/*      findiq, findrq, findbq, findsq, _findivq, _findrvq, rmq              */
/*                                                                           */
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

/*-----------------------------------------------------------------------------
 *  addiq  --  write an integer quantity at the end of a story;
 *             or, if a previous quantity of the same name is found,
 *             overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void  addiq(a_story, name, iq)
snodeptr  a_story;
char *name;
int  iq;
    {
    snodeptr  new_line;
    snodeptr  mksnode();
    snodeptr  find_qmatch();

    if ((new_line = find_qmatch(a_story, name)) != NULL)
        free(Text(new_line));
    else
	{
	new_line = mksnode();
	addchapter(a_story, new_line);
	}

    write_iq(new_line, name, iq);
    }

/*-----------------------------------------------------------------------------
 *  addrq  --  write a real (i.e. double precision floating point) quantity at
 *             the end of a story;
 *             or, if a previous quantity of the same name is found,
 *             overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void  addrq(a_story, name, rq)
snodeptr  a_story;
char *name;
real  rq;
    {
    snodeptr  new_line;
    snodeptr  mksnode();
    snodeptr  find_qmatch();

    if ((new_line = find_qmatch(a_story, name)) != NULL)
        free(Text(new_line));
    else
	{
	new_line = mksnode();

	addchapter(a_story, new_line);
	}
    write_rq(new_line, name, rq);
    }

/*-----------------------------------------------------------------------------
 *  addbq  --  write an boolean quantity at the end of a story;
 *             or, if a previous quantity of the same name is found,
 *             overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void  addbq(a_story, name, bq)
snodeptr  a_story;
char *name;
bool  bq;
    {
    snodeptr  new_line;
    snodeptr  mksnode();
    snodeptr  find_qmatch();

    if ((new_line = find_qmatch(a_story, name)) != NULL)
        free(Text(new_line));
    else
	{
	new_line = mksnode();
	addchapter(a_story, new_line);
	}

    write_bq(new_line, name, bq);
    }

/*-----------------------------------------------------------------------------
 *  addsq  --  write a string quantity at the end of a story;
 *             or, if a previous quantity of the same name is found,
 *             overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void  addsq(a_story, name, sq)
snodeptr  a_story;
char *name;
char *sq;
    {
    snodeptr  new_line;
    snodeptr  mksnode();
    snodeptr  find_qmatch();

    if ((new_line = find_qmatch(a_story, name)) != NULL)
        free(Text(new_line));
    else
	{
	new_line = mksnode();
	addchapter(a_story, new_line);
	}

    write_sq(new_line, name, sq);
    }

/*-----------------------------------------------------------------------------
 *  _addivq  --  write an integer vector quantity at the end of a story;
 *               or, if a previous quantity of the same name is found,
 *               overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void  _addivq(a_story, name, ivq, n_dim)
snodeptr  a_story;
char *name;
int  *ivq;
int  n_dim;
    {
    snodeptr  new_line;
    snodeptr  mksnode();
    snodeptr  find_qmatch();

    if ((new_line = find_qmatch(a_story, name)) != NULL)
        free(Text(new_line));
    else
	{
	new_line = mksnode();
	addchapter(a_story, new_line);
	}

    write_ivq(new_line, name, ivq, n_dim);
    }

/*-----------------------------------------------------------------------------
 *  _addrvq  --  write a real (i.e. double precision floating point) vector
 *               quantity at the end of a story;
 *               or, if a previous quantity of the same name is found,
 *               overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void  _addrvq(a_story, name, rvq, n_dim)
snodeptr  a_story;
char *name;
real  *rvq;
int  n_dim;
    {
    snodeptr  new_line;
    snodeptr  mksnode();
    snodeptr  find_qmatch();

    if ((new_line = find_qmatch(a_story, name)) != NULL)
        free(Text(new_line));
    else
	{
	new_line = mksnode();
	addchapter(a_story, new_line);
	}

    write_rvq(new_line, name, rvq, n_dim);
    }

/*-----------------------------------------------------------------------------
 *  findiq  --  searches for an integer quantity from a line in a story
 *-----------------------------------------------------------------------------
 */
bool  findiq(a_story, name, iq_ptr)
snodeptr  a_story;
char *name;
int *iq_ptr;
    {
    snodeptr  a_line;
    bool  qmatch();
    char *get_qstring();

    a_line = First_chapter(a_story);

    while (a_line != NULL)
	if (qmatch(a_line, name))
	    {
	    *iq_ptr = atoi(get_qstring(a_line));
	    return(TRUE);
	    }
	else 
	    a_line = Next_chapter(a_line);

    return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *  findrq  --  searches for a real quantity from a line in a story
 *-----------------------------------------------------------------------------
 */
bool  findrq(a_story, name, rq_ptr)
snodeptr  a_story;
char *name;
real *rq_ptr;
    {
    snodeptr  a_line;
    bool  qmatch();
    char *get_qstring();

    a_line = First_chapter(a_story);

    while (a_line != NULL)
	if (qmatch(a_line, name))
	    {
	    *rq_ptr = atof(get_qstring(a_line));
	    return(TRUE);
	    }
	else
	    a_line = Next_chapter(a_line);

    return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *  findbq  --  searches for a boolean quantity from a line in a story
 *-----------------------------------------------------------------------------
 */
bool  findbq(a_story, name, bq_ptr)
snodeptr  a_story;
char *name;
bool *bq_ptr;
    {
    snodeptr  a_line;
    bool  qmatch();
    char *bq_string;
    char *get_qstring();

    a_line = First_chapter(a_story);

    while (a_line != NULL)
	if (qmatch(a_line, name))
	    {
	    bq_string = get_qstring(a_line);

	    if (streq(bq_string, "TRUE"))
		{
	        *bq_ptr = TRUE;
		return(TRUE);
		}
	    else if (streq(bq_string, "FALSE"))
		{
	        *bq_ptr = FALSE;
		return(TRUE);
		}
	    else
	        error("findbq: \"%s\" is neither TRUE nor FALSE", name);
	    }
	else 
	    a_line = Next_chapter(a_line);

    return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *  findsq  --  searches for an string quantity from a line in a story
 *         note:
 *              the found string is returned including the '\n' before the
 *              terminating NULL; if this is not desired, a new function should
 *              be written which first copies the string found, and then 
 *              removes the '\n' from the end.
 *-----------------------------------------------------------------------------
 */
bool  findsq(a_story, name, sq_ptr)
snodeptr  a_story;
char *name;
char **sq_ptr;
    {
    snodeptr  a_line;
    bool  qmatch();
    char *get_qstring();

    a_line = First_chapter(a_story);

    while (a_line != NULL)
	if (qmatch(a_line, name))
	    {
	    *sq_ptr = get_qstring(a_line);
	    return(TRUE);
	    }
	else 
	    a_line = Next_chapter(a_line);

    return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *     The following macros specify the buffer lengths for strings containing
 *  different types of quantities.
 *     The numbers given below are overestimates; if memory size is an 
 *  important issue, these numbers can be decreases somewhat, after careful
 *  analysis of how much space is actually needed in the worst case (beware of
 *  details such as minus signs, exponents, minus signs in exponents as well as
 *  in the mantissa, etc.).
 *     The macros below are not given in the include file  storyq.h , since
 *  they should not be generally visible to user application programs, which do
 *  include  storyq.h . The alternative, of creating yet another small include
 *  file, which should logically be called  snodeq.h  perhaps, seemed
 *  a bit unwieldy and unnecessary.
 *-----------------------------------------------------------------------------
 */
#define  BYTE_LENGTH          8
#define  SAFE_INT_LENGTH     (5 + (BYTE_LENGTH*sizeof(int))/3)
#define  SAFE_REAL_LENGTH    (10 + (BYTE_LENGTH*sizeof(real))/3)
#define  SAFE_BOOL_LENGTH     7    /* this should hold any boolean value     */
#define  SAFE_STRING_LENGTH   3    /* this should hold the string terminator */
#define  EXTRA_LENGTH  5           /* for initial "  " and " = " */

/*-----------------------------------------------------------------------------
 *  _findivq  --  searches for an integer vector quantity from a line in a
 *                story
 *-----------------------------------------------------------------------------
 */
bool  _findivq(a_story, name, ivq, n_dim)
snodeptr  a_story;
char *name;
int *ivq;
int  n_dim;
    {
    int  i;
    char *ivq_string;
    char  ivq_string0[SAFE_INT_LENGTH];
    char  ivq_string1[SAFE_INT_LENGTH];
    char  ivq_string2[SAFE_INT_LENGTH];
    snodeptr  a_line;
    bool  qmatch();
    char *get_qstring();

    a_line = First_chapter(a_story);

    while (a_line != NULL)
	if (qmatch(a_line, name))
	    {
	    ivq_string = get_qstring(a_line);

	    if (n_dim == 3)
	        {
	        sscanf(ivq_string, "%s %s %s", ivq_string0,
		       ivq_string1, ivq_string2);
	        ivq[0] = atoi(ivq_string0);
	        ivq[1] = atoi(ivq_string1);
	        ivq[2] = atoi(ivq_string2);
		return(TRUE);
		}
	    else if (n_dim == 2)
	        {
	        sscanf(ivq_string, "%s %s", ivq_string0,
		       ivq_string1);
	        ivq[0] = atoi(ivq_string0);
	        ivq[1] = atoi(ivq_string1);
		return(TRUE);
		}
	    else if (n_dim == 1)
	        {
	        sscanf(ivq_string, "%s", ivq_string0);
	        ivq[0] = atoi(ivq_string0);
		return(TRUE);
		}
            else
	        error("_find_ivq: NDIM = %d not yet implemented; sorry!",
		      n_dim);
	    }
	else 
	    a_line = Next_chapter(a_line);

    return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *  _findrvq  --  searches for an real vector quantity from a line in a
 *                story
 *-----------------------------------------------------------------------------
 */
bool  _findrvq(a_story, name, rvq, n_dim)
snodeptr  a_story;
char *name;
real *rvq;
int  n_dim;
    {
    int  i;
    char *rvq_string;
    char  rvq_string0[SAFE_REAL_LENGTH];
    char  rvq_string1[SAFE_REAL_LENGTH];
    char  rvq_string2[SAFE_REAL_LENGTH];
    snodeptr  a_line;
    bool  qmatch();
    char *get_qstring();

    a_line = First_chapter(a_story);

    while (a_line != NULL)
	if (qmatch(a_line, name))
	    {
	    rvq_string = get_qstring(a_line);

	    if (n_dim == 3)
	        {
	        sscanf(rvq_string, "%s %s %s", rvq_string0,
		       rvq_string1, rvq_string2);
	        rvq[0] = atof(rvq_string0);
	        rvq[1] = atof(rvq_string1);
	        rvq[2] = atof(rvq_string2);
		return(TRUE);
		}
	    else if (n_dim == 2)
	        {
	        sscanf(rvq_string, "%s %s", rvq_string0,
		       rvq_string1);
	        rvq[0] = atof(rvq_string0);
	        rvq[1] = atof(rvq_string1);
		return(TRUE);
		}
	    else if (n_dim == 1)
	        {
	        sscanf(rvq_string, "%s", rvq_string0);
	        rvq[0] = atof(rvq_string0);
		return(TRUE);
		}
            else
	        error("_find_rvq: NDIM = %d not yet implemented; sorry!",
		      n_dim);
	    }
	else 
	    a_line = Next_chapter(a_line);

    return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *  rmq  --  removes the story line containing a quantity from a story
 *-----------------------------------------------------------------------------
 */
bool  rmq(a_story, name)
snodeptr  a_story;
char *name;
    {
    snodeptr  a_line;
    bool  qmatch();

    a_line = First_chapter(a_story);

    while (a_line != NULL)
	if (qmatch(a_line, name))
	    {
	    pusnode(a_line, a_story);
	    return(TRUE);
	    }
	else 
	    a_line = Next_chapter(a_line);

    return(FALSE);
    }

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*                                                                           */
/*                  Local Functions (alphabetically):                        */
/*                                                                           */
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

/*-----------------------------------------------------------------------------
 *  addstring  --  .....
 *-----------------------------------------------------------------------------
 */
local void  addstring(old_story, new_string, n)
snodeptr  old_story;
char *new_string;
int  n;
    {
    int  old_text_length;
    bool  contains_nl;
    bool  past_nl;                   /* not used here */
    snodeptr  lc;
    snodeptr  mksnode();

    lc = Last_chapter(old_story);

    if (isstory(lc))
        {
	lc = mksnode();
	addchapter(old_story, lc);
	give_text(lc, new_string, n);
	}
    else if (Text(lc) == NULL)
	give_text(lc, new_string, n);
    else
        {
        nl_detector(Text(lc), &old_text_length, &contains_nl, &past_nl);
        if (contains_nl)
            {
	    lc = mksnode();
	    addchapter(old_story, lc);
	    give_text(lc, new_string, n);
	    }
	else
	    extend_text(lc, old_text_length, new_string, n);
	}
    }

/*-----------------------------------------------------------------------------
 *  clean_story  --  remove the contents of the story.
 *            beware:
 *                   the memory locations for the old contents are freed.
 *-----------------------------------------------------------------------------
 */
local void  clean_story(a_story)
snodeptr  a_story;
    {
    snodeptr  a_chapter;
    snodeptr  new_chapter;
    snodeptr  mksnode();

    a_chapter = First_chapter(a_story);
    
    while ((new_chapter = Next_chapter(a_chapter)) != NULL)
	{
	rmsnode(a_chapter);
	a_chapter = new_chapter;
	}
    rmsnode(a_chapter);

    First_chapter(a_story) = Last_chapter(a_story) = mksnode();
    }

/*-----------------------------------------------------------------------------
 *  extend_text  --  ......
 *-----------------------------------------------------------------------------
 */
local void  extend_text(old_snode, old_n, new_string, new_n)
snodeptr  old_snode;
int  old_n;
char *new_string;
int  new_n;
    {
    char *new_text;

    new_text = malloc((old_n + new_n + 1) * sizeof(char));
    if (new_text == NULL)
	error("extend_text: not enough memory left");

    strncpy(new_text, Text(old_snode), old_n);
    strncpy(new_text + old_n, new_string, new_n);
    new_text[old_n + new_n] = NULL;
    
    free(Text(old_snode));
    Text(old_snode) = new_text;
    }

/*-----------------------------------------------------------------------------
 *  fgetchapter  --  read a whole story chapter from a specified file
 *              NOTE:
 *                   Current implementation simple, but not efficient: much
 *                   character copying occurs.
 *-----------------------------------------------------------------------------
 */
local snodeptr  fgetchapter(fp, first_line)
FILE  *fp;
snodeptr  first_line;
    {
    char *title_of();
    snodeptr  new_chapter;
    snodeptr  new_line;
    snodeptr  fgetline();
    
    new_chapter = mkstory();

    if (! is_story_begin_line(first_line))
	error("fgetchapter: first line not a story_begin_line");

    givetitle(new_chapter, title_of(first_line));

    while ((new_line = fgetline(fp)) != NULL && !is_story_end_line(new_line))
        {
	if (is_story_begin_line(new_line))
	    addchapter(new_chapter, fgetchapter(fp, new_line));
	else
	    {
	    addtext(new_chapter, Text(new_line));
	    rmsnode(new_line);
	    }
	}
    if (new_line == NULL)
        error("fgetchapter: new_line == NULL before end of chapter");

    return(new_chapter);
    }

/*-----------------------------------------------------------------------------
 *  fgetline  --  read a whole line from a specified file.  fgetline()
 *                returns a '\n'-terminated line, or returns NULL if the input
 *                starts with EOF, or causes an error message to be printed if
 *                 EOF  is reached in mid-line, before '\n'.
 *-----------------------------------------------------------------------------
 */
local snodeptr  fgetline(fp)
FILE  *fp;
    {
    char *s;
    int  c;
    int  i;
    int  n;
    snodeptr  new_snode;
    snodeptr  mksnode();
    line_piece  first_buffer;
    line_piece_ptr  a_line_piece;
    line_piece_ptr  mk_line_piece();
    
    new_snode = mksnode();

    a_line_piece = &first_buffer;
    s = Text_piece(a_line_piece);
/*
 * initialization, as in mk_line_piece(), has to be done here too:
 */
    s[0] = NULL;
    Next_line_piece(a_line_piece) = NULL;

    i = 0;
    while ((c = getc(fp)) != EOF && c != '\n')
	{
	s[i % LINE_PIECE_LENGTH] = c;

	if (++i % LINE_PIECE_LENGTH == 0)
	    {
	    Next_line_piece(a_line_piece) = mk_line_piece();
	    a_line_piece = Next_line_piece(a_line_piece);
	    s = Text_piece(a_line_piece);
	    }
	}

    if (c == '\n')
	s[i % LINE_PIECE_LENGTH] = c;
    else                                              /* c == EOF */
        {
	if (i == 0)                                       /* empty line */
	    return(NULL);
	else                                              /* incomplete line */
            error("fgetline: EOF reached in mid-line");
	}

    n = length_of_line_text(&first_buffer);

    Text(new_snode) = malloc((n+1) * sizeof(char));
    if (Text(new_snode) == NULL)
	error("fgetline: not enough memory left");

    a_line_piece = &first_buffer;
    i = 0;
    while (n - i > LINE_PIECE_LENGTH)
        {
	strncpy(Text(new_snode) + i, Text_piece(a_line_piece),
	        LINE_PIECE_LENGTH);
	i += LINE_PIECE_LENGTH;
	a_line_piece = Next_line_piece(a_line_piece);
	}
    strncpy(Text(new_snode) + i, Text_piece(a_line_piece), n - i);
    Text(new_snode)[n] = NULL;

    if (Next_line_piece(&first_buffer))
	rm_line_piece(Next_line_piece(&first_buffer));	

    return(new_snode);
    }

/*-----------------------------------------------------------------------------
 *  find_qmatch  --  within a story, find a line which contains a physical
 *                   quantity with a matching name; or return NULL if such a
 *                   line is not found.
 *-----------------------------------------------------------------------------
 */
local snodeptr  find_qmatch(a_story, name)
snodeptr  a_story;
char *name;
    {
    snodeptr  a_chapter;
    bool  qmatch();

    a_chapter = First_chapter(a_story);

    while (a_chapter != NULL)
        {
	if (qmatch(a_chapter, name))
	    return(a_chapter);
	else
	    a_chapter = Next_chapter(a_chapter);
	}

    return(NULL);
    }

/*-----------------------------------------------------------------------------
 *  fput_headline  --  write the begin delimiter (including title, if present)
 *                    of a story
 *-----------------------------------------------------------------------------
 */
local void  fput_headline(fp, a_story)
FILE  *fp;
snodeptr  a_story;
    {
    if (Text(a_story))
	fprintf(fp, "%c%s\n", Story_begin_char, Text(a_story));
    else
	fprintf(fp, "%c\n", Story_begin_char);
    }

/*-----------------------------------------------------------------------------
 *  fput_tailline  --  write the end delimiter of a story
 *-----------------------------------------------------------------------------
 */
local void  fput_tailline(fp)
FILE  *fp;
    {
    fprintf(fp, "%c\n", Story_end_char);
    }

/*-----------------------------------------------------------------------------
 *  fputchapter  --  ...
 *-----------------------------------------------------------------------------
 */
local void  fputchapter(fp, a_snode)
FILE  *fp;
snodeptr  a_snode;
    {
    snodeptr  a_chapter;

    if (! isstory(a_snode))
        error("fputchapter: snode presented is not a story");

    fput_headline(fp, a_snode);

    a_chapter = First_chapter(a_snode);
    while (a_chapter != NULL)
	{
	fputsnode(fp, a_chapter);
	a_chapter = Next_chapter(a_chapter);
	}

    fput_tailline(fp);
    }

/*-----------------------------------------------------------------------------
 *  fputline  --  write a whole line to a file
 *-----------------------------------------------------------------------------
 */
local void  fputline(fp, a_line)
FILE  *fp;
snodeptr  a_line;
    {
    if (Text(a_line) != NULL)
	fprintf(fp, "%s", Text(a_line));
    }

/*-----------------------------------------------------------------------------
 *  fputsnode  --  write a snode and its descendents to a file
 *-----------------------------------------------------------------------------
 */
local void  fputsnode(fp, a_snode)
FILE  *fp;
snodeptr  a_snode;
    {
    if (isstory(a_snode))
	fputchapter(fp, a_snode);
    else
	fputline(fp, a_snode);
    }

/*-----------------------------------------------------------------------------
 *  get_qstring  --  
 *               note:
 *                    no check for presence of '=' ; BEWARE
 *-----------------------------------------------------------------------------
 */
local char *get_qstring(a_line)
snodeptr  a_line;
    {
    char *c;
    
    c = Text(a_line);
    
    while (*(++c) != '=')
	;
    while (*(++c) == ' ')
        ;

    return(c);
    }

/*-----------------------------------------------------------------------------
 *  give_text  --  ......
 *-----------------------------------------------------------------------------
 */
local void  give_text(a_snode, new_string, n)
snodeptr  a_snode;
char *new_string;
int  n;
    {

    Text(a_snode) = malloc((n + 1) * sizeof(char));
    if (Text(a_snode) == NULL)
	error("give_text: not enough memory left");

    strncpy(Text(a_snode), new_string, n);
    Text(a_snode)[n] = NULL;
    }

/*-----------------------------------------------------------------------------
 *  length_of_line_text  --  
 *-----------------------------------------------------------------------------
 */
local int  length_of_line_text(a_line_piece)
line_piece_ptr  a_line_piece;
    {
    int  n;
    char *s;

    n = 0;
    while(Next_line_piece(a_line_piece) != NULL)
        {
	n += LINE_PIECE_LENGTH;
	a_line_piece = Next_line_piece(a_line_piece);
	}

    s = Text_piece(a_line_piece);
    while (*(s++) != '\n')
        n++;
    n++;                           /* to count the final '\n' */

    return(n);
    }

/*-----------------------------------------------------------------------------
 *  mk_line_piece  --  
 *-----------------------------------------------------------------------------
 */
local line_piece_ptr  mk_line_piece()
    {
    line_piece_ptr lpp;

    lpp = (line_piece_ptr)malloc(sizeof(line_piece));
    if (lpp == NULL)
	error("mk_line_piece: not enough memory left");

    Text_piece(lpp)[0] = NULL;
    Next_line_piece(lpp) = NULL;

    return(lpp);
    }

/*-----------------------------------------------------------------------------
 *  mksnode  --  
 *-----------------------------------------------------------------------------
 */
local snodeptr  mksnode()
    {
    snodeptr new_snode;

    new_snode = (snodeptr)malloc(sizeof(snode));
    if (new_snode == NULL)
	error("mksnode: not enough memory left");

    Text(new_snode) = NULL;
    Next_chapter(new_snode) = NULL;
    First_chapter(new_snode) = Last_chapter(new_snode) = NULL;

    return(new_snode);
    }

/*-----------------------------------------------------------------------------
 *  nl_detector  --  .......
 *                   *nptr is the length of the string up to and including the
 *                   first new-line character, if present; if not, it is the
 *                   length of the string (not including the NULL terminator). 
 *-----------------------------------------------------------------------------
 */
local void  nl_detector(a_string, nptr, has_nl_ptr, past_nl_ptr)
char *a_string;
int *nptr;
bool *has_nl_ptr;
bool *past_nl_ptr;
    {
    int  i;

    i = 0;
    while (a_string[i] != NULL && a_string[i] != '\n')
        i++;

    if (a_string[i] == '\n')
        {
        i++;
        *has_nl_ptr = TRUE;
	if (a_string[i] != NULL)
	    *past_nl_ptr = TRUE;
	else
	    *past_nl_ptr = FALSE;
	}
    else
	*has_nl_ptr = *past_nl_ptr = FALSE;

    *nptr = i;
    }

/*-----------------------------------------------------------------------------
 *  pusnode  --  as rmsnode, but restores neighbor pointer links to pave over
 *              the gap created by removing the snode  old_snode .
 *-----------------------------------------------------------------------------
 */
local void  pusnode(old_snode, parent_snode)
snodeptr  old_snode;
snodeptr  parent_snode;
    {
    snodeptr  a_chapter, b_chapter;

    if (First_chapter(parent_snode) == old_snode)
        First_chapter(parent_snode) = Next_chapter(old_snode);
    else
        {
	a_chapter = First_chapter(parent_snode);
	while ((b_chapter = Next_chapter(a_chapter)) != old_snode)
	    a_chapter = b_chapter;
	Next_chapter(a_chapter) = Next_chapter(old_snode);
        }
	
    if (Last_chapter(parent_snode) == old_snode)
        Last_chapter(parent_snode) = NULL;

    rmsnode(old_snode);
    }

/*-----------------------------------------------------------------------------
 *  qmatch  --  checks whether a line contains a physical quantity with a
 *              matching name; if so, returns TRUE, otherwise FALSE.
 *         note: 
 *              leading blanks are discarded, both in the name and in the text
 *              of a_line ; trailing blanks in the text up to the equal sign
 *              are also discarded; if the equal sign is missing (i.e. if the
 *              line is not a proper quantity line) FALSE is returned.
 *-----------------------------------------------------------------------------
 */
local bool  qmatch(a_line, name)
snodeptr  a_line;
char *name;
    {
    int  i, j;
    char *s;

    if ((s = Text(a_line)) == NULL)
        return(FALSE);
    i = 0;
    while (s[i] == ' ')
        i++;

    j = 0;
    while (name[j] != NULL)
        {
	if (name[j] == ' ')
            j++;
	else 
	    break;
	}

    while (name[j] != NULL)
        {
	if (s[i] != name[j])
            return(FALSE);
	i++;
	j++;
	}

    while (s[i] == ' ')
        i++;

    if (s[i] != '=')
        return(FALSE);
    if (s[++i] != ' ')
        return(FALSE);

    return(TRUE);
    }

/*-----------------------------------------------------------------------------
 *  rm_line_piece  --  
 *-----------------------------------------------------------------------------
 */
local void  rm_line_piece(lpp)
line_piece_ptr  lpp;
    {
    if (Next_line_piece(lpp))
	rm_line_piece(Next_line_piece(lpp));

    free(Text_piece(lpp));
    free(lpp);
    }

/*-----------------------------------------------------------------------------
 *  rmsnode  --  ....
 *              ....
 *-----------------------------------------------------------------------------
 */
local void  rmsnode(old_snode)
snodeptr  old_snode;
    {
    snodeptr  a_chapter;
    snodeptr  new_chapter;

    if (isstory(old_snode))
        {
	a_chapter = First_chapter(old_snode);
	while ((new_chapter = Next_chapter(a_chapter)) != NULL)
	    {
	    rmsnode(a_chapter);
	    a_chapter = new_chapter;
	    }
	rmsnode(a_chapter);
	}
    else
        {
	if (Text(old_snode))
            free (Text(old_snode));
        free(old_snode);
        }
    }

/*-----------------------------------------------------------------------------
 *  title_of  --  ......
 *           note:
 *                this clobbers the string  a_line  by writing a NULL over the
 *                final  \n
 *-----------------------------------------------------------------------------
 */
local char *title_of(a_line)
snodeptr  a_line;
    {
    char *title_string;

    if (! is_story_begin_line(a_line))
	error("title_of: line presented is not a story_begin_line");

    if (Text(a_line)[1] == NULL)
        return(NULL);

    Text(a_line)[strlen(Text(a_line))-1] = NULL;
    
    return(Text(a_line) + 1);
    }

/*-----------------------------------------------------------------------------
 *  write_iq  --  write an integer quantity to a line.
 *-----------------------------------------------------------------------------
 */
local void  write_iq(a_line, name, iq)
snodeptr  a_line;
char *name;
int  iq;
    {
    char *new_string;
    int  new_string_length;
    
    new_string_length = EXTRA_LENGTH + strlen(name) + SAFE_INT_LENGTH;

    Text(a_line) = malloc(new_string_length * sizeof(char));
    if (Text(a_line) == NULL)
	error("write_iq: not enough memory left");

    sprintf(Text(a_line), "  %s = %d\n", name, iq);
    }

/*-----------------------------------------------------------------------------
 *  write_rq  --  write a real quantity to a line (see write_iq).
 *-----------------------------------------------------------------------------
 */
local void  write_rq(a_line, name, rq)
snodeptr  a_line;
char *name;
real  rq;
    {
    char *new_string;
    int  new_string_length;
    
    new_string_length = EXTRA_LENGTH + strlen(name) + SAFE_REAL_LENGTH;

    Text(a_line) = malloc(new_string_length * sizeof(char));
    if (Text(a_line) == NULL)
	error("write_rq: not enough memory left");

    sprintf(Text(a_line), "  %s = %23.16e\n", name, rq);
    }

/*-----------------------------------------------------------------------------
 *  write_bq  --  write a boolean quantity to a line (see write_iq).
 *-----------------------------------------------------------------------------
 */
local void  write_bq(a_line, name, bq)
snodeptr  a_line;
char *name;
bool  bq;
    {
    char *new_string;
    int  new_string_length;
    
    new_string_length = EXTRA_LENGTH + strlen(name) + SAFE_BOOL_LENGTH;

    Text(a_line) = malloc(new_string_length * sizeof(char));
    if (Text(a_line) == NULL)
	error("write_bq: not enough memory left");

    if (bq)
        sprintf(Text(a_line), "  %s = TRUE\n", name);
    else
        sprintf(Text(a_line), "  %s = FALSE\n", name);
    }

/*-----------------------------------------------------------------------------
 *  write_sq  --  write a string quantity to a line (see write_iq).
 *-----------------------------------------------------------------------------
 */
local void  write_sq(a_line, name, sq)
snodeptr  a_line;
char *name;
char *sq;
    {
    char *new_string;
    int  new_string_length;

    new_string_length = EXTRA_LENGTH + strlen(name) + strlen(sq)
    		        + SAFE_STRING_LENGTH;

    Text(a_line) = malloc(new_string_length * sizeof(char));
    if (Text(a_line) == NULL)
	error("write_sq: not enough memory left");

    sprintf(Text(a_line), "  %s = %s\n", name, sq);
    }

/*-----------------------------------------------------------------------------
 *  write_ivq  --  write an integer vector quantity to a line (see write_iq).
 *                 note: the SAFE_INT_LENGTH is extended by two for the 
 *                       separating space between the values in the vector.
 *-----------------------------------------------------------------------------
 */
local void  write_ivq(a_line, name, ivq, n_dim)
snodeptr  a_line;
char *name;
int  *ivq;
int  n_dim;
    {
    char *new_string;
    int  new_string_length;
    
    new_string_length = EXTRA_LENGTH + strlen(name)
		        + n_dim * (SAFE_INT_LENGTH + 2);

    Text(a_line) = malloc(new_string_length * sizeof(char));
    if (Text(a_line) == NULL)
	error("write_ivq: not enough memory left");

    if (n_dim == 3)
        sprintf(Text(a_line), "  %s = %d  %d  %d\n", name, ivq[0], ivq[1],
	        ivq[2]);
    else if (n_dim == 2)
        sprintf(Text(a_line), "  %s = %d  %d\n", name, ivq[0], ivq[1]);
    else if (n_dim == 1)
        sprintf(Text(a_line), "  %s = %d\n", name, ivq[0]);
    else
        error("write_ivq: NDIM = %d not yet implemented; sorry!", n_dim);
    }

/*-----------------------------------------------------------------------------
 *  write_rvq  --  write a real vector quantity to a line (see write_ivq).
 *-----------------------------------------------------------------------------
 */
local void  write_rvq(a_line, name, rvq, n_dim)
snodeptr  a_line;
char *name;
real  *rvq;
int  n_dim;
    {
    char *new_string;
    int  new_string_length;
    
    new_string_length = EXTRA_LENGTH + strlen(name) 
                        + n_dim * (SAFE_REAL_LENGTH + 2);

    Text(a_line) = malloc(new_string_length * sizeof(char));
    if (Text(a_line) == NULL)
	error("write_rvq: not enough memory left");

    if (n_dim == 3)
        sprintf(Text(a_line), "  %s = %23.16e %23.16e %23.16e\n", name,
	        rvq[0], rvq[1], rvq[2]);
    else if (n_dim == 2)
        sprintf(Text(a_line), "  %s = %23.16e  %23.16e\n", name, rvq[0],
		rvq[1]);
    else if (n_dim == 1)
        sprintf(Text(a_line), "  %s = %23.16e\n", name, rvq[0]);
    else
        error("write_rvq: NDIM = %d not yet implemented; sorry!", n_dim);
    }

/*===========================================================================*/

#ifdef TOOLBOX

/*-----------------------------------------------------------------------------
 *  main  --  driver to test addchapter() and addtext().
 *            Stories are read from the standard input. 
 *            If no story is provided, a new (short) story is created, with a
 *            title taken from the first command line argument (if present)
 *            and line contents taken from the following command line 
 *            arguments (if present).
 *            If one or more stories are read in, then the stories are written
 *            out as one (nested) story, again with a title taken from the 
 *            first command line argument (if present), to which are added as
 *            final lines the remaining command line arguments (if present).
 *
 *            Example:
 *                    story | tee tmp
 *                    a
 *                    story "my first line!" | tee -a tmp
 *                    a
 *                    story a short poem < tmp
 *            note:
 *                 the single a is typed above to indicate that no new story
 *                 appears on the standard input; anything else than a story
 *                 begin delimiter will do.  A more direct way would be to type
 *                 ^d, but that is more dangerous: typing it twice may log you
 *                 out.
 *         NOTE: 
 *              if a command line argument is "rq", the next command line
 *              arguments are interpreted as providing the name of an real
 *              quantity and its value, respectively; this information is then
 *              converted and written into the next line of the story at hand.
 *-----------------------------------------------------------------------------
 */

#include <getparam.h>

string defv[] = {
    "in=-\n        Input file [stdin]",
    "out=-\n       Output file [stdout]",
    "title=\n      Optional new title",
    "rq=\n         Optional new quantity",
    "VERSION=1.1\n 9-oct-90 PJT",
    NULL,
};

main(argc, argv)
int  argc;
char  *argv[];
{
    int  i;
    snodeptr  a_story;
    snodeptr  the_story;
    string title, rq;
    char *strchr(), *cp;
    stream instr, outstr;

    initparam(argv,defv);
    title = getparam("title");
    rq = getparam("rq");
    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");

    the_story = mkstory();

    while ((a_story = fgetstory(instr)) != NULL)
	addchapter(the_story, a_story);

    if (title)
        givetitle(the_story, title);

    if (*rq) {
        cp = strchr(rq,'=');
        if (cp==NULL)
            warning("No value supplied for quantity %s",rq);
        else {
            *cp = NULL;
            cp++;
	    addrq(the_story, rq, atof(cp));
	}
    } 
#if 0
    addtext(the_story, argv[i++]);
    addtext(the_story, "\n");
#endif

    fputstory(outstr,the_story);
}

#endif

/*===========================================================================*/

/* endof: story.c */
