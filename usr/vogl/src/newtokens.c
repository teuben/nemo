#include <stdio.h>
#include "vogl.h"

static TokList		*current;

/*
 * newtokens
 *
 *	returns the space for num tokens
 */
Token *
newtokens(num)
	int	num;
{
	TokList	*tl;
	Token	*addr;
	int	size;

	if (vdevice.tokens == (TokList *)NULL || num >= MAXTOKS - current->count) {
		if ((tl = (TokList *)malloc(sizeof(TokList))) == (TokList *)NULL)
			verror("newtokens: malloc returns NULL");

		if (vdevice.tokens != (TokList *)NULL)
			current->next = tl;
		else 
			vdevice.tokens = tl;

		tl->count = 0;
		tl->next = (TokList *)NULL;
		if (num > MAXTOKS)
			size = num;
		else
			size = MAXTOKS;
		if ((tl->toks = (Token *)malloc(size * sizeof(Token))) == (Token *)NULL)
			verror("newtokens: malloc returns NULL");

		current = tl;
	}

	addr = &current->toks[current->count];
	current->count += num;

	return(addr);
}
