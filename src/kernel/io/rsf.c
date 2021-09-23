/*
 * RSF.C: read a structured file.
 *
 * V1.0  Panayotis Skordos 6/87	first implementation
 * V2.0  Joshua Barnes 4/88	new filestruct package
 * V2.1  Peter Teuben 12-sep-90	helpvec
 * V2.2  Peter Teuben  8-dec-90	bug reading multiline character strings
 * V2.3  PJT          12-jun-91 nemo_main
 *     a              19-feb-92 usage
 *     b	      19-may-92 malloc() -> allocate()
 *     c              19-aug-92 removed some lint complaints
 *     d              20-feb-94 ansi
 *     e               5-feb-95 no more ARGS (alpha complained)
 *     f              12-mar-97 NULL -> 0
 */

#include <stdinc.h>
#include <getparam.h>
#include <ctype.h>
#include <filestruct.h>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			input file name (ascii structured file)",
    "out=???\n			output file name (binary structured file)",
    "VERSION=2.3f\n		15-mar-05 PJT",
    NULL,
};

string usage = "Read an ascii structured file into binary structured file";

/* local functions */

bool parse_line (stream, string, string, int *, byte **);
void parse_data_item (stream, string, string, int *, byte **);
bool read_token (stream , string);
bool is_any    (string);
bool is_char   (string);
bool is_byte   (string);
bool is_short  (string);
bool is_int    (string);
bool is_long   (string);
bool is_float  (string);
bool is_double (string);
bool is_set    (string);
bool is_tes    (string); 
void safe_read_string     (stream, char *, int, char);
void safe_read_name       (stream, char *);
void safe_read_dimensions (stream, int *);
void safe_read_byte       (stream, void *);
void safe_read_char       (stream, char *);
void safe_read_short      (stream, short *);
void safe_read_int        (stream, int *);
void safe_read_long       (stream, long *);
void safe_read_float      (stream, float *);
void safe_read_double     (stream, double *); 

void nemo_main()
{
    stream istr, ostr;
    char typ[MaxTagLen], tag[MaxTagLen];
    int dims[MaxVecDim];
    byte *data;

    dprintf(2,"RSF: MaxTagLen = %d  MaxVecDim = %d\n",MaxTagLen, MaxVecDim);
    istr = stropen(getparam("in"), "r");
    ostr = stropen(getparam("out"), "w");
    while (parse_line(istr, typ, tag, dims, &data)) {
	if (streq(typ, SetType))
	    put_set(ostr, tag);
	else if (streq(typ, TesType))
	    put_tes(ostr, tag);
	else
	    put_data_sub(ostr, tag, typ, data, dims, FALSE);
	if (data)
	    free((char *)data);
    }
}

/*
 * PARSE_LINE: scan input file for one item worth of info.
 * Returns: TRUE for successful completion, FALSE at E.O.F..
 * The pointer returned should be de-allocated after use.
 */

bool parse_line(
    stream istr,
    string type,
    string name,
    int dims[],
    byte **data)          /* a pointer to a char pointer */
{
    char keyw[MaxTagLen];
  
    if (! read_token(istr, keyw))	/* read in type keyword             */
	return (FALSE);
    if (is_any(keyw)) {              /* classify and dispatch on keyword */
	parse_data_item(istr, keyw, name, dims, data);
	strcpy(type, AnyType);
    } else if (is_char(keyw)) {
	parse_data_item(istr, keyw, name, dims, data);
	strcpy(type, CharType);
    } else if (is_byte(keyw)) {
	parse_data_item(istr, keyw, name, dims, data);
	strcpy(type, ByteType);
    } else if (is_short(keyw)) {
	parse_data_item(istr, keyw, name, dims, data);
	strcpy(type, ShortType);
    } else if (is_int(keyw)) {
	parse_data_item(istr, keyw, name, dims, data);
	strcpy(type, IntType);
    } else if (is_long(keyw)) {
	parse_data_item(istr, keyw, name, dims, data);
	strcpy(type, LongType);
    } else if (is_float(keyw)) {
	parse_data_item(istr, keyw, name, dims, data);
	strcpy(type, FloatType);
    } else if (is_double(keyw)) {
	parse_data_item(istr, keyw, name, dims, data);
	strcpy(type, DoubleType);
    } else if (is_set(keyw)) {
	safe_read_name(istr, name);	/*   that's all the parsing needed  */
	*data = NULL;
	strcpy(type, SetType);
    } else if (is_tes(keyw)) {
	*data = NULL;
	strcpy(type, TesType);
    } else if (streq(keyw,"."))    /* . . . is used by tsf for incomplete data */
	error("parse_line: possibly incomplete data in or after %s",name);
    else 
	error("parse_line: unknown keyword in tsf file %s", keyw);
    return TRUE;
}

/*
 * PARSE_DATA_ITEM: handle scan and allocation of data item.
 */

void parse_data_item(
    stream istr,
    string keyw,
    string name,
    int dims[],
    byte **data)
{
    int nbytes, i;

    safe_read_name(istr, name);			/* read the name first */
    safe_read_dimensions(istr, dims);		/* read dimensions if any   */
    nbytes = 1;					/* compute nbytes from dims */
    for (i = 0; dims[i] !=  0; i++)
	nbytes *= dims[i];
    if (is_any(keyw) || is_byte(keyw)) {	/* dispatch on type of data */
	*data = (byte *) allocate(nbytes * sizeof(byte));
	for (i = 0; i < nbytes; i++) 
	    safe_read_byte(istr, &((byte *) *data)[i]);
    } else if (is_char(keyw)) {
	*data = (byte *) allocate(nbytes * sizeof(char));
	safe_read_string(istr, (char *)*data, nbytes, dims[0] == 0 ? '\'' : '\"');
					/* use s.q. for scalar, else d.q.   */
    } else if (is_short(keyw)) {
	*data = (byte *) allocate(nbytes * sizeof(short));
	for (i = 0; i < nbytes; i++)
	    safe_read_short(istr, &((short *) *data)[i]);   /* PPAP */
    } else if (is_int(keyw)) {
	*data = (byte *) allocate(nbytes * sizeof(int));
	for (i = 0; i < nbytes; i++)
	  safe_read_int(istr, &((int *) *data)[i]);         /* PPAP */
    } else if (is_long(keyw)) {
	*data = (byte *) allocate(nbytes * sizeof(long));
	for (i = 0; i < nbytes; i++)
	  safe_read_long(istr, &((long *) *data)[i]);       /* PPAP */
    } else if (is_float(keyw)) {
	*data = (byte *) allocate(nbytes * sizeof(float));
	for (i = 0; i < nbytes; i++)
	  safe_read_float(istr, &((float *) *data)[i]);     /* PPAP */
    } else if (is_double(keyw)) {
	*data = (byte *) allocate(nbytes * sizeof(double));
	for (i = 0; i < nbytes; i++)
	    safe_read_double(istr, &((double *) *data)[i]); /* PPAP */
    } else
	error("parse_data_item: illegal type");
}

string str_uppercase(string str)
{
    int i, c;
    permanent char buf[MaxTagLen];

    for(i = 0; str[i] != 0; i++) {
	c = str[i];
	buf[i] = (islower(c) ? toupper(c) : c);
    }
    buf[i] = 0;
    return buf;	    		/* note volatile storage used!      */
}

/*
 * Predicates for keyword types; case insensitive.
 */

bool is_any(char *keyw)
{
    return (streq(str_uppercase(keyw), "ANY"));
}

bool is_char(char *keyw)
{
    return (streq(str_uppercase(keyw), "CHAR"));
}

bool is_byte(char *keyw)
{
    return (streq(str_uppercase(keyw), "BYTE"));
}

bool is_short(char *keyw)
{
    return (streq(str_uppercase(keyw), "SHORT"));
}

bool is_int(char *keyw)
{
    return (streq(str_uppercase(keyw), "INT"));
}

bool is_long(char *keyw)
{
    return (streq(str_uppercase(keyw), "LONG"));
}

bool is_float(char *keyw)
{
    return (streq(str_uppercase(keyw), "FLOAT"));
}

bool is_double(char *keyw)
{
    return (streq(str_uppercase(keyw), "DOUBLE"));
}

bool is_set(char *keyw)
{
    return (streq(str_uppercase(keyw), "SET"));
}

bool is_tes(char *keyw)
{
    return (streq(str_uppercase(keyw), "TES"));
}

bool get_char_ok(stream fp, char thechar)
{
    char temp;

    fscanf(fp, " %c", &temp);
    return (temp == thechar); 
}

/*
 * SAFE_READ_STRING: rather smart routine understands backslashes
 * (skips them and reads in the next) and fills with NULLS if there
 * are not enough chars in the double-quote-delimeted space.
 *	dec-90   bug (?) removed; can now read multiline strings   PJT
 *	mar-94   ansi (bug in sun cc)
 */

void safe_read_string(stream fp, char *word, int len, char delim)
{
    bool more_data=TRUE, just_newline = FALSE;
    int  i;
    char char_read;
  
    if (! get_char_ok(fp, delim))
	error("safe_read_string: opening quote not read");
    for (i = 0; i < len; i++) {
	if (more_data) {
	    safe_read_char(fp, &char_read);             /* always read char */
            while (char_read == ' ' && just_newline) {/* if at start of line */
                safe_read_char(fp, &char_read);    /* read until non-space */
            }
            just_newline = FALSE;        /* always flag we're in the middle */
	}
	if (char_read == '\n') {                /* if newline read in       */
            i--;                                /*  decrease counter and    */
            just_newline = TRUE;                /*  set flag for reader and */
            continue;                           /*  read next char again    */
        }
	if (char_read == '\\')			/* if backslash read in     */
	    safe_read_char(fp, &char_read);	/*   read whatevers next    */
	else if (char_read == '\"') {		/* if closing quote seen    */
	    more_data = FALSE;			/*   read no more data      */
	    char_read = 0;			/*   store NULL instead     */
	}
	word[i] = char_read;
    }
    if (i != len)  /* due to terminating zero, len is 1 longer than string is */
        warning("safe_read_string: expected %d characters, read %d:\n%s",
                    len,i,word);
    i = strlen(word)+1;
    if (i != len)
        warning("safe_read_string: expected %d characters, saw %d:\n%s",
                    len,i,word);
    if (more_data && ! get_char_ok(fp, delim))
	error("safe_read_string: closing quote not read:\n %s",word);
}

/*
 * standard read abstractions for parse.c
 */

bool read_token(stream fp, char *word)
{
    return (fscanf(fp, "%s", word) != EOF);
}

void safe_read_name(stream fp, char *word)
{
    if (fscanf(fp, " %[^[ \n]", word) == EOF)
				/* the ^[ stops at the dimensions, if any */
	error("safe_read_name: unexpected EOF");
}

void safe_read_dimensions(stream fp, int dims[MaxVecDim])
{
    int i;

    for (i = 0; i < MaxVecDim; i++)
	dims[i] = 0;			/* zero to mark end of dimensions;  *
					 * note scalars have dims[0] = 0    */
    i = 0;
    fscanf(fp, "[%d]", &dims[i]);   /* try to read the first dimension */
    while ((i < MaxVecDim) && (dims[i] != 0)) {
                                     /* if you haven't failed, try again */
	i += 1;
	fscanf(fp, "[%d]", &dims[i]);
    }
}

void safe_read_byte(stream fp, void *vword)
{
    int temp;
    char *word = (char *) vword;

    if (fscanf(fp, "%o", &temp) == EOF)
	error("read_byte: EOF found");
    *word = ((char) temp);
}

void safe_read_char(stream fp, char *word)
{
    char temp;

    if (fscanf(fp, "%c", &temp) == EOF)
	error("read_char: EOF found");
    *word = temp;
}

void safe_read_short(stream fp, short *word)
{
    short temp;

    if (fscanf(fp, "%ho", &temp) == EOF)
	error("read_short: EOF found");
    *word = temp;
}

void safe_read_int(stream fp, int *word)
{
    int temp;

    if (fscanf(fp, "%o", &temp) == EOF)
	error("read_int: EOF found");
    *word = temp;
}

void safe_read_long(stream fp, long *word)
{
    long temp;

    if (fscanf(fp, "%lo", &temp) == EOF)
	error("read_long: EOF found");
    *word = temp;
}

void safe_read_float(stream fp, float *word)
{
    float temp;

    if (fscanf(fp, "%f", &temp) == EOF)
	error("read_float: EOF found");
    *word = temp;
}

void safe_read_double(stream fp, double *word)
{
    double temp;

    if (fscanf(fp, "%lf", &temp) == EOF)
	error("read_double: EOF found");
    *word = temp;
}
