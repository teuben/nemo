#include "vogl.h"

static	int	sync = 1;
/*
 * backbuffer
 *
 *	swap drawing to backbuffer - returns -1 if no
 * backbuffer is available.
 */
void
backbuffer(yes)
	int	yes;
{
	Token	*tok;

	if (!vdevice.initialised)
		verror("backbuffer: vogl not initialised.");

	if (vdevice.inobject) {
		tok = newtokens(2);
		tok[0].i = BACKBUFFER;
		tok[1].i = yes;
		return;
	}

	if (vdevice.attr->a.mode == SINGLE)
		return;

	if (yes) {
		if ((*vdevice.dev.Vbackb)() < 0)
			verror("device doesn't support double buffering\n");
		
		vdevice.inbackbuffer = 1;
		vdevice.sync = 0;
	} else
		vdevice.inbackbuffer = 0;
	
}

/*
 * frontbuffer
 *
 *	start drawing in the front buffer again. This
 * will always work!
 */
void
frontbuffer(yes)
	int	yes;
{
	Token	*tok;

	if (!vdevice.initialised)
		verror("frontbuffer: vogl not initialised.");

	if (vdevice.inobject) {
		tok = newtokens(2);
		tok[0].i = FRONTBUFFER;
		tok[1].i = yes;
		return;
	}

	if (vdevice.attr->a.mode == SINGLE)
		return;

	if (yes) {
		(*vdevice.dev.Vfrontb)();
		vdevice.inbackbuffer = 0;
		vdevice.sync = sync;
	} else
		backbuffer(1);
}

/*
 * swapbuffers
 *
 *	swap the back and front buffers - returns -1 if
 * no backbuffer is available.
 */
void
swapbuffers()
{
	Token	*tok;

	if (vdevice.inobject) {
		tok = newtokens(1);
		tok[0].i = SWAPBUFFERS;
		return;
	}

	if (!vdevice.initialised)
		verror("swapbuffers: vogl not initialised.");

/* TEST HACK
	if (vdevice.inbackbuffer != 1)
		verror("swapbuffers: double buffering not initialised.\n");
*/

	if ((*vdevice.dev.Vswapb)() < 0)
		verror("swapbuffers device doesn't support double buffering\n");

}

/*
 * doublebuffer()
 *
 *	Flags our intention to do double buffering....
 *	tries to set it up etc etc...
 */
void
doublebuffer()
{
	if (!vdevice.initialised)
		verror("doublebuffer: vogl not initialised.");

	if ((*vdevice.dev.Vbackb)() < 0)
		verror("device doesn't support double buffering\n");

	vdevice.inbackbuffer = 1;
	sync = vdevice.sync;
	vdevice.sync = 0;
}

/*
 * singlebuffer()
 *
 *	Goes back to singlebuffer mode....(crock)
 */
void
singlebuffer()
{ 
	if (vdevice.attr->a.mode == SINGLE)
		return;

	(*vdevice.dev.Vfrontb)();

	vdevice.inbackbuffer = 0;

	if (sync)
		vdevice.sync = 1;
}
