#ifndef lint
static char sccsid[] = "@(#)get_view_surface.c 1.1 86/07/07 SMI";
#endif

/*
 *	get_view_surface  --  Determines from command-line arguments and
 *			      the environment a reasonable view surface
 *			      for a SunCore program to run on.
 *	Mods by Peter Wisnovsky --
 *		Will not process command line args;
 *		If the device is a color device, will allocate
 *		(but not initialize) a three bit color table.
 *	16-feb-90:	new SUN3COLOR option for SUN OS 4.0.3 ??
 *	14-mar-90:	make GCC happy [PJT]
 *      30-may-90:      added fbtype = FBTYPE_SUNFAST_COLOR	PJT
 *	28-may-92:	fixed bug (but not quite tested yet) on argv NULL	PJT
 */

#include <sunwindow/window_hs.h>
#include <sys/file.h>
#include <sys/ioctl.h>
#include <sun/fbio.h>
#include <stdio.h>

int bw1dd();		/* All device-independent/device-dependent   */
int bw2dd();		/* routines are referenced in this function.      */
int cg1dd();		/* This means the linker will pull in all of them */
int cg2dd();
int cg3dd();
int cg4dd();
int gp1dd();
int pixwindd();
int cgpixwindd();
int gp1pixwindd();

static struct vwsurf nullvs = NULL_VWSURF;

static char *devchk;
static int devhaswindows;
static int chkdevhaswindows();
int my_get_view_surface(vsptr, argv, cmssize)
struct vwsurf *vsptr;
char **argv;
int cmssize;
	{
	int devfnd, fd, fbtype;
	char *wptr, dev[DEVNAMESIZE], *getenv();
	struct screen screen;

	*vsptr = nullvs;
	devfnd = FALSE;
	if ((argv) && (*argv))
		/*
		If command-line arguments are passed, process them using
		win_initscreenfromargv (see the Programmer's Reference Manual
		for the Sun Window System).  The only option used by
		get_view_surface is the -d option, allowing the user to
		specify the display device on which to run.
		*/
		{
		win_initscreenfromargv(&screen, argv);
		if (screen.scr_fbname[0] != '\0')
			{
			/* -d option was found */
			devfnd = TRUE;
			strncpy(dev, screen.scr_fbname, DEVNAMESIZE);
			/*
			Check to see if this device has a window system
			running on it.  If so devhaswindows will be TRUE
			following the call to win_enumall.  win_enumall is
			a function in libsunwindow.a.  It takes a function
			as its argument, and applies this function to every
			window being displayed on any screen by the window
			system.  To do this it opens each window and passes
			the windowfd to the function.  The enumeration
			continues until all windows have been tried or the
			function returns TRUE.
			*/
			devchk = dev;
			devhaswindows = FALSE;
			win_enumall(chkdevhaswindows);
			}
		}
	if (!devfnd)
		/* No -d option was specified */
		if (wptr = getenv("WINDOW_ME"))
			{
			/*
			Running in the window system.  Find the device from
			which this program was started.
			*/
			devhaswindows = TRUE;
			if ((fd = open(wptr, O_RDWR, 0)) < 0)
			    {
			    fprintf(stderr, "get_view_surface: Can't open %s\n",
				wptr);
			    return(1);
			    }
			win_screenget(fd, &screen);
			close(fd);
			strncpy(dev, screen.scr_fbname, DEVNAMESIZE);
			}
		else
			{
			/*
			Not running in the window system.  Assume device is
			/dev/fb.
			*/
			devhaswindows = FALSE;
			strncpy(dev, "/dev/fb", DEVNAMESIZE);
			}
	/* Now have device name.  Find device type. */
	if ((fd = open(dev, O_RDWR, 0)) < 0)
		{
		fprintf(stderr, "get_view_surface: Can't open %s\n", dev);
		return(1);
		}
	if ((fbtype = pr_getfbtype_from_fd(fd)) == -1)
		{
		fprintf(stderr, 
			"get_view_surface: pr_getfbtype_from_fd() failed for %s\n",
			dev);
		close(fd);
		return(1);
		}
	close(fd);
	/* Now have device type and know if window system is running on it. */
	if (devhaswindows)
		switch(fbtype)
			{
		case FBTYPE_SUN1BW:
		case FBTYPE_SUN2BW:
			vsptr->dd = pixwindd;
			break;
		case FBTYPE_SUN1COLOR:
		case FBTYPE_SUN2COLOR:
		case FBTYPE_SUN3COLOR:		/* new for SUN OS 4.0.3 */
		case FBTYPE_SUN4COLOR:
		case FBTYPE_SUNFAST_COLOR:	/* added for SUN OS 4.0.3 */
			vsptr->dd = cgpixwindd;
			vsptr->cmapsize = cmssize;
			break;
		case FBTYPE_SUN2GP:
			vsptr->dd = gp1pixwindd;
			break;
		default:
			fprintf(stderr,
			"get_view_surface: %s is unknown fbtype %d\n", 
			dev,fbtype);
			return(1);
			}
	else
		switch(fbtype)
			{
		case FBTYPE_SUN1BW:
			vsptr->dd = bw1dd;
			break;
		case FBTYPE_SUN2BW:
			vsptr->dd = bw2dd;
			break;
		case FBTYPE_SUN1COLOR:
			vsptr->dd = cg1dd;
			break;
		case FBTYPE_SUN2COLOR:
			vsptr->dd = cg2dd;
			break;
		case FBTYPE_SUN3COLOR:	/* new for SUN OS 4.0.3 */
			vsptr->dd = cg3dd;
			break;
		case FBTYPE_SUN4COLOR:
			vsptr->dd = cg4dd;
			break;
		case FBTYPE_SUN2GP:
			vsptr->dd = gp1dd;
			break;
		default:
			fprintf(stderr,
			"get_view_surface: %s is unknown fbtype %d\n", 
			dev,fbtype);
			return(1);
			}
	/* Now SunCore device driver pointer is set up. */
	if (!devhaswindows || devfnd)
		/*
		If no window system on device or -d option was specified,
		tell SunCore which device.  Otherwise, let SunCore figure
		out the device itself from WINDOW_GFX so the default
		window will be used if desired.
		*/
		strncpy(vsptr->screenname, dev, DEVNAMESIZE);
	return(0);
	}

static int chkdevhaswindows(windowfd)
int windowfd;
	{
	struct screen windowscreen;

	win_screenget(windowfd, &windowscreen);
	if (strcmp(devchk, windowscreen.scr_fbname) == 0)
		{
		/*
		If this window is on the display device we are checking, set
		the flag TRUE.  Return TRUE to terminate the enumeration.
		*/
		devhaswindows = TRUE;
		return(TRUE);
		}
	return(FALSE);
	}
