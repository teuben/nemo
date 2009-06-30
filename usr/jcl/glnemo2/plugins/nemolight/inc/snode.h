/*
 *  snode.h: include file for story.c
 *.............................................................................
 *    version 1:  June 1989   Piet Hut               email: piet@iassns.bitnet
 *                            Institute for Advanced Study, Princeton, NJ, USA
 *.............................................................................
 *     This file should not be included in user application programs; see the
 *  corresponding heading in  story.h .
 *.............................................................................
 *  see also: story.c, story.h, storyio.h, storyq.h
 *.............................................................................
 */

#include "story.h"

/*-----------------------------------------------------------------------------
 *  The following macros provide access to the components of a story structure.
 *  They provide a barrier of data abstraction, which allows future changes in
 *  the actual implementation of such a structure without affecting the use of
 *  the structure components in existing programs, as long as the macros below
 *  are changed in accordance to the changes in the underlying structure.
 *
 *  Text  --  points to the character string containing the story title, or
 *            in case of a proper line, the text of the line.
 *  Prev_story  --  points to the previous story
 *  Next_story  --  points to the next story
 *  First_chapter  --  points to the first substory
 *  Last_chapter  --  points to the last substory
 *-----------------------------------------------------------------------------
 */
#define  Text(ptr)             ((ptr)->text)               /* type: charptr  */
#define  Prev_chapter(ptr)     ((ptr)->prev_chapter)       /* type: storyptr */
#define  Next_chapter(ptr)     ((ptr)->next_chapter)       /* type: storyptr */
#define  First_chapter(ptr)    ((ptr)->first_chapter)      /* type: storyptr */
#define  Last_chapter(ptr)     ((ptr)->last_chapter)       /* type: storyptr */

/*-----------------------------------------------------------------------------
 *  isstory  --  TRUE (i.e. non-zero) if a story, FALSE otherwise
 *-----------------------------------------------------------------------------
 */
#define  isstory(ptr)     First_chapter(ptr)

/* endof: snode.h */
