#include "vogl.h"
#include <stdio.h>

#define LT_SIZE	11

typedef struct l {
	short		n;
	unsigned	ls;
	struct l	*next;
} LT;

static LT	*ls_table[LT_SIZE] = { (LT *)0, };

/*
 * deflinestyle
 *
 *	 Define a line style
 * 	 There can be 2*166 line styles!
 */
void
deflinestyle(n, ls)
	short		n;
	Linestyle	ls;
{
	LT	*l;

	/*
	 * If it's already there, then redefine it...
	 */
	for (l = ls_table[n % LT_SIZE]; l != (LT *)NULL; l = l->next)
		if (l->n == n) {
			if (n == 0) /* Not allowed  to redefine 0 (0xffff) */
				verror("vogl: deflinestyle can't redefine 0");

			l->ls = ls;
			return;
		}

	/*
	 * Othersize, we have to add it...
	 */

	l = (LT *)vallocate(sizeof(LT));
	l->n = n;
	l->ls = ls;
	l->next = ls_table[n % LT_SIZE];
	ls_table[n % LT_SIZE] = l;
}

/*
 * setlinestyle
 *
 *	Set the current linestyle...
 */
void
setlinestyle(n)
	short	n;
{
	LT	*l;
	Token	*tok;

	if (!vdevice.initialised)
		verror("setlinestyle: vogl not initialised");

	if (vdevice.inobject) {
		tok = newtokens(2);

		tok[0].i = LINESTYLE;
		tok[1].i = n;
		return;
        }

	/*
	 * Find it....
	 */
	for (l = ls_table[n % LT_SIZE]; l != (LT *)NULL; l = l->next) {
		if (l->n == n) {
			vdevice.attr->a.ls = l->ls;
			(*vdevice.dev.Vsetls)(l->ls);
			return;
		}
	}

	verror("vogl: undefined linestyle used.");
}

/*
 * linewidth
 *
 *	Set's the line width in pixels if we can....
 */
void
linewidth(w)
	short	w;
{
	Token	*tok;

	if (!vdevice.initialised)
		verror("linewidth: vogl not initialised");

	if (vdevice.inobject) {
		tok = newtokens(2);

		tok[0].i = LINEWIDTH;
		tok[1].i = w;
		return;
        }

	vdevice.attr->a.lw = w;
	(*vdevice.dev.Vsetlw)((int)vdevice.attr->a.lw);
}

/*
 * linewidthf
 *
 *	The 'float' version of the above. (Used for antialiased lines on
 *	a real SGI, here, we'll just round it to the nearest int and be
 *	the same as linewidth).
 */
void
linewidthf(w)
	float	w;
{
	linewidth((short)(w + 0.5));
}


