/*
 *  storyq.h: for internal translations of physical Quantities from stories
 *.............................................................................
 *    version 1:  May 1989   Piet Hut               email: piet@iassns.bitnet
 *                           Institute for Advanced Study, Princeton, NJ, USA
 *.............................................................................
 *     Normally stories can only be added to, not read or modified.  The only
 *  exception is through the use of reading and writing (physical) quantities,
 *  which can be of different types. In addition, they can be scalar (by 
 *  default) or vector (with a length NDIM). Examples of quantity types are
 *  integer, real, boolean and string type (the latter one allowing a very
 *  general type of "quantity": any text string is allowed).
 *     Quantities are represented in single story lines, in the form 
 *  "name = value" . Here "name" is an arbitrary string of characters, which
 *  can include letters, numbers, white spaces, and other symbols (except "=").
 *  The value associated with the name is given in "value".
 *     Example of a single integer quantitiy:
 *  	   particle number = 23
 *     Example of a real vector quantitiy:
 *  	   velocity = 0.123 -8.765 12.345
 *     Example of a string quantitiy:
 *         cluster name = NGC 1234
 *.............................................................................
 *  see also: story.c, story.h, snode.h, storyio.h
 *.............................................................................
 */

#include "vectmath.h"                   /* for vector operations             */

/*-----------------------------------------------------------------------------
 *  addivq, ...  --  names for functions contained in  story.c , but which
 *                   are depended on the value NDIM.
 *                      These macros are introduced here to avoid recompilation
 *                   of the file  story.c  when a different value of NDIM is
 *                   used in different user programs. 
 *                       This is accomplished by passing the value of NDIM to
 *                   the library in a way which is hidden for the user.
 *-----------------------------------------------------------------------------
 */
#define addivq(a, b, c)    (_addivq((a), (b), (c), NDIM))
#define addrvq(a, b, c)    (_addrvq((a), (b), (c), NDIM))
#define findivq(a, b, c)   (_findivq((a), (b), (c), NDIM))
#define findrvq(a, b, c)   (_findrvq((a), (b), (c), NDIM))

bool  findiq();
bool  findrq();
bool  findbq();
bool  findsq();
bool  _findivq();
bool  _findrvq();

/* endof: storyq.h */
