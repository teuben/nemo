/*
 *  MATCH: Minimal match of sub-strings in 'STRING' with substrings
 *		in the compare string 'COMPAR'. The result is returned
 *		in 'MASK' (if a match is possible with the I-th substring
 *		in 'COMPAR', then the I-th bit in 'MASK' will be set.
 *		The comparison is only done for all printable ascii
 *		ASCII characters except for the blank and the comma,
 *		since they are used as substring delimiters.
 *
 *	Adapted by Peter Teuben from an original program by Kor Begeman 
 * 	written in SHELTRAN for GIPSY, (c) University of Groningen 1985.
 *	
 *	notes:  maximum sub-strings is 32 (depending on bit-length of integer)
 *		needs substring-comparison routine 'STREQ' or so 
 *		
 *	output: match returns:
 *			 1	no error
 *			 0      length of input string is zero
 *			-1	length of compare string is zero
 *			-2	a sub-string in STRING can be matched in more
 *				than one sub-string in COMPAR
 *			-3	no sub-string in STRING can be matched to any
 *				sub-string in COMPAR.
 *
 *  ===> Bug: multiple options doesn't work  <===
 *            i.e. with return code -2 it always seems mask=0x00
 *
 *	 2-jun-88 last time tinkered with it			    PJT
 *	15-may-91 finally included n NEMO - changed par order       PJT
 *	25-feb-92 happy gcc2.0					    PJT
 *       4-apr-92 fixed bug that it now returns 1 on success        PJT
 *      19-jun-92 renamed string not to clash with NEMO's type      PJT
 *      27-feb-94 ansi + <stdinc.h> now
 *	12-apr-95 no more ARGS, full prototypes
 *	 7-apr-01 gcc warning         
 */

#include <stdinc.h>

int match(string, string, int *);
int partialstreq (char *, int ,int, char *, int, int);

#define INTLEN 32	/* length of an integer in bits, i.e. max ss's */

#define C_SPACE 0x20
#define C_DEL   0x7F
#define C_COMMA 0x2C

int match (char *istring, char *compar, int *mask)
{
    int nsub, isub, ichr,
	ncom, 			/* length of COMPAR */
	nmat, 			/* length of STRING */
	nummat, 		/* number of substrings in STRING */
	numcom,			/* number of substrings in COMPAR */
	i,j,i1,i2,j1,j2, imatch;
    char c;
    int  lcom[INTLEN+1], ecom[INTLEN+1], lmat[INTLEN+1], emat[INTLEN+1];	
	
    ncom=strlen(compar);
    nmat=strlen(istring);
    if (ncom==0) 
        return(-1);
    else if (nmat==0)
	return(0);

    nsub=0;
    isub=0;
    for (ichr=0; ichr<ncom; ichr++) {		/*   parse COMPAR  */
        c=compar[ichr];
	if ( (c>C_SPACE) && (c<C_DEL) && (c!=C_COMMA) ) /* check legal char */
	    isub++;
	else if (isub>0) {			/* separator ? */
	    nsub++;			    /* number of substrings */
	    lcom[nsub]=isub;	    /* length of substring  */
	    ecom[nsub]=ichr;	    /* end position of ss   */
	    isub=0;
	}
    }
    if (isub>0) {
	nsub++;
	lcom[nsub]=isub;
	ecom[nsub]=ichr;
    }
    numcom=nsub;			/* number of sub-strings in COMPAR */
#ifdef DEBUG	
    printf ("Number of substrings in COMPAR: %d\n",numcom);
    for (i=1; i<=numcom; i++) 
	printf ("lcom ecom = %d %d\n",lcom[i],ecom[i]);
#endif
    nsub=0;
    isub=0;
    for (ichr=0; ichr<nmat; ichr++) {		/* parse STRING */
	c=istring[ichr];
	if ( (c>C_SPACE) && (c<C_DEL) && (c!=C_COMMA) )
	    isub++;
	else if (isub>0) {
	    nsub++;
	    lmat[nsub]=isub;
	    emat[nsub]=ichr;
	    isub=0;
	}
    }
    if (isub>0) {
        nsub++;
	lmat[nsub]=isub;
	emat[nsub]=ichr;
    }
    nummat=nsub;
#ifdef DEBUG
    printf ("Number of substrings in STRING: %d\n",nummat);
    for (i=1; i<=nummat; i++) 
        printf ("lmat emat = %d %d\n",lmat[i],emat[i]);
#endif
		
    *mask = 0;
    for (i=1; i<=nummat; i++) {		/* do for every ss in STRING */
	i1=emat[i]-lmat[i];
	i2=emat[i]-1;
	imatch= -1;			    /* set for no match */
	for (j=1; j<=numcom; j++) {	 /* check every ss in COMPAR */
	    j1=ecom[j]-lcom[j];
	    j2=j1+lmat[i]-1;
	    if (lmat[i] <= lcom[j]) {    /* match string short enough */
	        if (partialstreq(compar,j1,j2,istring,i1,i2)) {
		    if (imatch<0)        /* is it a first match? */
		        imatch=j-1;       /* set substring bit code */
		    else
			return(-2);
		}
	    }
	}
	if (imatch<0) 
	    return (-3);
	else {
#ifdef DEBUG			
	    printf ("Masking 0x%x with 0x%x\n",*mask,imatch);
#endif
	    *mask |= (1<<imatch);	/* set appropriate bit */
	}
    }
    return(1);		/* successfull return */
}


/*
 * PARTIALSTREQ: compares two strings partially
 *	 
 *	input  a,i1,i2:  string and its starting end ending position to test
 *	       b,j1,j2:  test-string, its starting and ending position
 *		
 *	returns   1  (true)  : if strings have minimum equality
 *		  0  (false) : not equal
 *	It stops comparing when either one of 'a' or 'b' strings is
 *	exhausted according to the specified bounds. The caller is
 *	is responsible for len(a) {<|=|>} len(b).
 */
 
int partialstreq (char *a, int i1,int i2, char *b,int j1,int j2)
{
    while ( (i1<=i2) && (j1<=j2))
        if (a[i1++]!=b[j1++])
            return 0;			/*  not equal */
    return 1;				/*    equal   */
}


#ifdef TESTBED

#include <stdinc.h>
#include <getparam.h>

string defv[] = {
    "istring=foo\n                          Test options",
    "compar=foo bar and barred galaxies\n   Compare string",
    "VERSION=1\n                            27-feb-94 PJT",   
    NULL,
};
string usage="testbed for match.c";

nemo_main()
{
	string compare=getparam("compar");
        string word=getparam("istring");
	int i,j,mask;

	printf ("ISTRING :%s: with length %d\n",word,strlen(word));
	printf ("COMPAR  :%s: with length %d\n",compare,strlen(compare));
	printf ("Match returns: %d\n",match(word,compare,&mask));	
	printf ("mask=0x%x\n",mask);
}

#endif
