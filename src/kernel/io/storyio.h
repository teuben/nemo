/*
 *  storyio.h: data structure declarations & macros for story input/output
 *.............................................................................
 *    version 0:  Jul 1990   Piet Hut               email: piet@iassns.bitnet
 *                           Institute for Advanced Study, Princeton, NJ, USA
 *.............................................................................
 *     This file should be included in  story.c  but not in user application
 *  programs.  The latter should only access I/O functions contained in 
 *   story.c .  This provides an abstraction barrier, allowing user-transparent
 *  future changes in story I/O.
 *.............................................................................
 *  see also: story.c, story.h, snode.h, storyq.h
 *.............................................................................
 */

/*-----------------------------------------------------------------------------
 *  story_begin_string, story_end_string  --  story delimiters
 *-----------------------------------------------------------------------------
 */
#define  Story_begin_char	'('
#define  Story_end_char 	')'

/*-----------------------------------------------------------------------------
 *  ...
 *-----------------------------------------------------------------------------
 */
#define  is_story_begin_line(ptr)  (*(Text(ptr)) == Story_begin_char)

#define  is_story_end_line(ptr)	   (*(Text(ptr)) == Story_end_char)

/*-----------------------------------------------------------------------------
 *  LINE_PIECE_LENGTH  --  Size of the character array which holds the string
 *                         in a line piece structure (see below).
 *                         The actual value is arbitrary; a shorter value
 *                         would break up lines in more pieces, thereby 
 *                         incurring more function execution overhead; a longer
 *                         value would reserve more memory than is necessary.
 *                         In either case, the actual line length is completely
 *                         unrestricted, and will span as many line pieces as
 *                         needed.
 *-----------------------------------------------------------------------------
 */
#define  LINE_PIECE_LENGTH            80

/*-----------------------------------------------------------------------------
 *  line_piece  --  A data structure used as a buffer for reading story lines
 *             note:
 *                  the introduction of the additional name _line_piece is
 *                  necessary because a more straightforward construction such
 *                  as 
 *			 typedef  struct
 *			     {
 *			     char  text_piece[LINE_LENGTH];
 *			     line_piece_ptr  next_line_piece;
 *			     } line_piece, *line_piece_ptr;
 *
 *                  is syntactically wrong, since the definition of 
 *                  "line_piece_ptr" here would be self-referential.
 *-----------------------------------------------------------------------------
 */
struct _line_piece
    {
    char  text_piece[LINE_PIECE_LENGTH];
    struct  _line_piece *next_line_piece;
    };

typedef  struct _line_piece  line_piece, *line_piece_ptr;

/*-----------------------------------------------------------------------------
 *  The following macros provide access to the components of a line piece
 *  structure. They provide a barrier of data abstraction, which allows future
 *  changes in the actual implementation of such a structure without affecting 
 *  the use of the structure components in existing programs, as long as the 
 *  macros below are changed in accordance to the changes in the underlying 
 *  structure.
 *
 *  Text_piece   --  points to the character string contained in a line or line
 *                   piece
 *  Next_line_piece --  points to the next line in the line list.
 *-----------------------------------------------------------------------------
 */
#define  Text_piece(ptr)       ((ptr)->text_piece)         /* type: charptr  */
#define  Next_line_piece(ptr)  ((ptr)->next_line_piece)     
						     /* type: line_piece_ptr */
/* endof: storyio.h */
