/*
 *  story.h: data structure declarations & macros for story telling
 *.............................................................................
 *    version 0:  July 1990   Piet Hut               email: piet@iassns.bitnet
 *                            Institute for Advanced Study, Princeton, NJ, USA
 *	18-nov-91	included <malloc.h>		PJT
 *.............................................................................
 *     In this file, the basic data structure of a snode is defined. This is
 *  the minimum amount of information which should be included in a user
 *  application program which deals with stories.  However, the user should not
 *  address the components of a snode directly.  This is left to the standard
 *  story manipulation functions, given in the file  story.c .  These functions
 *  address the snode components only indirectly, through accessor macros,
 *  which are defined in the include file  snode.h , which should not be
 *  included in user programs.  In this way there are two levels of abstraction
 *  barriers separating the user from the underlying data structure, thereby
 *  allowing future changes in implementation as well as accessing which can
 *  be completely transparent for the user, as long as the names and results
 *  of the user accessible functions in  story.c  remain invariant.
 *     Externally available functions in  story.c  can only add to stories,
 *  not read arbitrary lines of them, except indirectly, through the reading
 *  and writing of physical quantities (see  storyq.h ).
 *     This file, together with  storyq.h , are the only headerfiles included
 *  in user application programs; the two headerfiles  snode.h  and 
 *   storyio.h  are reserved for inclusion in  story.c  only.
 *.............................................................................
 *  see also: story.c, snode.h, storyio.h, storyq.h
 *.............................................................................
 */

#include <stdinc.h>
#include <malloc.h>

/*-----------------------------------------------------------------------------
 *  snode     --  a node in a tree structure, containing one character pointer
 *                (pointing to a line of text) and three pointers to similar
 *                nodes.  Only the pointers to the next snode and to the first
 *                chapter (the first substory) are really necessary; the 
 *                pointer to the last chapter is redundant, but included for
 *                efficiency.
 *  snodeptr  --  a pointer to a snode
 *           note:
 *                the introduction of the additional name _snode is necessary
 *                because a more straightforward construction such as
 *
 *			 typedef  struct
 *			     {
 *			     snodeptr  next_snode;
 *			     snodeptr  first_chapter;
 *			     snodeptr  last_chapter;
 *			     char  *text;
 *			     } snode, *snodeptr;
 *
 *                is syntactically wrong, since the definition of 
 *                "snodeptr" here would be self-referential. 
 *                The name  _snode  is not used anywhere hereafter, once
 *                the names  snode  and  snodeptr  are defined.
 *-----------------------------------------------------------------------------
 */
struct _snode
    {
    struct _snode *next_chapter;
    struct _snode *first_chapter;
    struct _snode *last_chapter;
    char  *text;
    };

typedef  struct _snode  snode, *snodeptr;

/*-----------------------------------------------------------------------------
 *  story, storyptr  --  to provide compatibility with starlab versions before
 *                       July 1990
 *-----------------------------------------------------------------------------
 */
#define  story     snode
#define  storyptr  snodeptr

/*-----------------------------------------------------------------------------
 *  General declarations, freeing the user from the need to declare these
 *  functions (the remaining declarations can be found in  storyq.h ).
 *-----------------------------------------------------------------------------
 */
snodeptr  mkstory();
snodeptr  cpstory();
snodeptr  findchapter();
snodeptr  find_first_chapter();
snodeptr  find_next_chapter();
snodeptr  take_out_last_chapter();
snodeptr  fgetstory();
snodeptr  getstory();
bool  clean_chapter();
bool  checktitle();
bool  is_empty_story();

/* endof: story.h */
