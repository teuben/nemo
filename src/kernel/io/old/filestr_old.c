/*=== This file contains pieces of old code that no longer should have to be
      maintained ===*/
/*
 * Old-style types and flags - if used at all
 *	This section comes from filesecret.c and is now taken
 *	out as we don't expect to support old filestruct
 *	anymore.
 */

#if defined(FILESTR_OLD)

#define OldTagFlag    001       /* item has a tag field */
#define OldVecFlag    002       /* item is a vector type */

#define OldAnyType    000       /* anything at all */
#define OldCharType   001       /* printable chars */
#define OldByteType   002       /* unprintable chars */
#define OldShortType  003       /* short integers */
#define OldIntType    004       /* standard integers */
#define OldLongType   005       /* long integers */
#define OldFloatType  006       /* short floating */
#define OldDoubleType 007       /* long floating */
#define OldSetType    010       /* begin compound item */
#define OldTesType    011       /* end of compound item */

#endif



	/* piece of code taken out of gethdr() in filesecret.c */
	/* it does a last resort check on header to see if it */
	/* could be old-style header */
#if 0
    else {					/* old-style (probably)     */
	oflg = num >> 8;			/*   extract old flag value */
	otyp = num & 0377;			/*   and old type code      */
	if ((oflg & ~(OldTagFlag | OldVecFlag)) ||
	      otyp < OldAnyType || OldTesType < otyp ||
	    					/*   catch bad flag or type */
	        (oflg & OldTagFlag) == 0 && otyp != OldTesType)
						/*   catch untagged data    */
	    error("gethdr: bad magic: %o\n", num);
	num = (oflg & OldVecFlag) ? PlurMagic : SingMagic;
						/*   deduce new-style magic */
	typ = convtype(otyp);			/*   convert to new type    */
	if (fread(&olen, sizeof(short), 1, str) != 1)
	    					/*   read trailing length   */
	    error("gethdr: EOF reading old length\n");
	if (olen != baselen(typ))		/*   check lengths match    */
	    error("gethdr: type %s old length %d inconsistent\n", typ, olen);
	if (firsttime)				/*   first time through?    */
	    fprintf(stderr, "[filestruct: reading old-style file]\n");
	firsttime = FALSE;
    }
#endif


#if defined(FILESTR_OLD)
/*
 * CONVTYPE: use tl_tab to get new type from old one.
 */

local string convtype(otyp)
int otyp;
{
    char *copxstr();

    return ((string) copxstr(tl_tab[otyp].tl_typ, sizeof(char)));
}
#endif

