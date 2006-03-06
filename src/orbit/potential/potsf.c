/*
 * POTSF: extract potential (potname/pars/file) from a dataset
 *        *** cloned from tsf ***
 *
 *   0.1  13-sep-05        Example for Peyaud 
 *   0.2                   added first=
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <history.h>

string defv[] = {
  "in=???\n                     input file name ",
  "first=t\n                    Only print the first one found",
  "VERSION=0.2\n		13-sep-05 PJT ",
  NULL,
};

string usage="find potential parameters in a structured file";

string cvsid="$Id$";




local stream instr;			/* input stream from struct. file   */

local bool Qpot   = FALSE;
local bool Qfirst = TRUE;
local int npot = 0;

/* local functions */
void   print_item   (string);
void   print_set    (string);
void   print_tes    (string);
void   print_header (string, string, int *);
void   print_data   (string, string, int *);
bool   outstr       (string);
void   end_line     (void);

#define BUFLEN   256


void nemo_main()
{
  string *tags, *hist, *sp;

  instr = stropen(getparam("in"), "r");
  Qfirst = getbparam("first");
  
  while ((tags = list_tags(instr)) != NULL) {
    print_item(*tags);
    free(*tags);
    free((char *)tags);
  }
  if (npot==0) {
    warning("No potential found, trying history....");
    warning("eventually the program will parse these....");
    hist = ask_history();
    for (sp=hist; *sp; sp++) {
      dprintf(0,"%s\n",*sp);
    }
  }
}

void print_item(string tag)
{
  string type, *tags, *tp;
  int *dims;

  type = get_type(instr, tag);
  if (streq(type, SetType)) {
    get_set(instr, tag);
    print_set(tag);
    if (streq(tag,"Potential")) {
      Qpot = TRUE;
      npot++;
    }
    tags = list_tags(instr);
    for (tp = tags; *tp != NULL; tp++)
      print_item(*tp);
    get_tes(instr, tag);
    if (streq(tag,"Potential")) {
      if (Qfirst) stop(0);
      Qpot = FALSE;
    }
    print_tes(tag);
    for (tp = tags; *tp != NULL; tp++)
      free(*tp);
    free((char *)tags);
  } else {
    dims = get_dims(instr, tag);
    print_header(tag, type, dims);
    print_data(tag, type, dims);
    end_line();
    if (dims != NULL)
      free((char *)dims);
  }
  free(type);
}

void print_set(string tag)
{
  end_line();
}
    
void print_tes(string tag)
{
  end_line();
}

void print_header(string tag, string type, int *dims)
{
  if (!Qpot) return;
    
  if (streq(tag,"Name")) printf("potname=");
  if (streq(tag,"Pars")) printf("potpars=");
  if (streq(tag,"File")) printf("potfile=");
}

void print_data(string tag, string type, int *dims)
{
  int dlen;
  char buf[BUFLEN];  /* danger: buffer overflow */
  byte *dat, *dp;
  
  dlen = get_dlen(instr, tag);
  dat = (byte*) allocate(dlen);
  get_data_sub(instr, tag, type, dat, dims, FALSE);
  
  for (dp = dat; dp < dat + dlen; ) {		/* loop over data array     */
    if (streq(type, AnyType)) {		/*   output generic data?   */
      dp += sizeof(byte);
    } else if (streq(type, CharType)) {	/*   output readable chars? */
      if (Qpot) sprintf(buf, "%c", *((char *) dp));
      dp += sizeof(char);
    } else if (streq(type, ByteType)) {	/*   output bytes of data?  */
      dp += sizeof(byte);
    } else if (streq(type, ShortType)) {	/*   output short integers? */
      dp += sizeof(short);
    } else if (streq(type, IntType)) {	/*   output standard ints?  */
      dp += sizeof(int);
    } else if (streq(type, LongType)) {	/*   output long integers?  */
      dp += sizeof(long);
    } else if (streq(type, FloatType)) {	/*   output floating point? */
      dp += sizeof(float);
    } else if (streq(type, DoubleType)) {	/*   output double numbers? */
      dp += sizeof(double);
    } else
      error("print_data: type %s unknown\n", type);
    if (Qpot)
      outstr(buf);
  }
  free((char *)dat);
}


bool outstr(string str)
{
  if (Qpot) printf("%s", str);                          /* output string */
}

void end_line()
{
  if (Qpot) printf("\n");
}

