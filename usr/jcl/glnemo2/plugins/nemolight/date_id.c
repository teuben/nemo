/* 
 *  DATE_ID:    return pointer to (static) string containing
 *              current date and time in the format
 *              "Sun Sep 16 01:03:52 1973" (see: ctime(3))
 *	See also:  on BSD:  	time(3)
 *		   on SYSV: 	time(2)
 *
 *  Written:    xx-xxx-89   but with SUN bug		PJT
 *  mods:	11-mar-90   consistent on BSD and SYSV	PJT
 *			    by defining, not using #include
 *               1-jul-90   turned sun back on		PJT
 *		 9-dec-90   the 'real time;' declaration caused trouble PJT
 *		19-feb-94   ANSI
 *              12-jan-99   kkkkon solaris and w/ ccmalloc ctime() crashed PJT
 *              24-nov-03   %ld instead of %d
 */

#include <stdinc.h>
#include <time.h>	/* time(), ctime() */

string date_id(void)
{
    permanent char did[32];
#if 0    
    char *ct;
    time_t clock, tloc;

    ct = ctime(&lt);
    strcpy(did,ct);
    did[strlen(did)-1] = '\0';  /* get rid of last \n; at position 24 */
#else
    long lt = time(0);

    sprintf(did,"DATE_ID=%ld",lt);
#endif
    return(did);
}

#if defined(TESTBED)

main()
{
    printf("date_id=%s\n",date_id());
}
#endif
