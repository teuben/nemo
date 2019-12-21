/*
 * BSF: benchmark regression of the floating point data in a 
 *      structured file
 *
 * 24-nov-2019   PJT     Created
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#ifdef unix
# include <sys/ioctl.h>
# include <stdio.h>
# include <unistd.h>
#endif
#include <moment.h>

string defv[] = {
  "in=???\n              input file name",
  "test=\n               If given, this is the regression test",
  "fmt=%g\n              output format for floating point values",
  "eps=\n                Accuracy comparison (not implemented)",
  "label=\n              Override the in= filename in reporting",
  "ignore=cputime\n      Items to ignore in checksum",
  "VERSION=1.1\n         21-dec-2019 PJT ",
  NULL,
};

string usage="benchmark stats on all floating point values in a (binary) structured file";

string cvsid="$Id$";

local stream instr;	
local Moment m;

local string ignore = "cputime";

/* local functions */
void   accum_item   (string);
void   accum_data   (string, string, int *);
string find_name    (string);

#define BUFLEN   512


void nemo_main()
{
  string *tags;
  string infile = getparam("in");
  string fmt = getparam("fmt");
  string test = getparam("test");
  string ignore = getparam("ignore");
  string label;
  char fmt4[32];
  char current[128];

  dprintf(2,"TSF: MaxSetLen = %d\n",MaxSetLen);
  ini_moment(&m,2,0);
  sprintf(fmt4,"%s %s %s %s %%d",fmt,fmt,fmt,fmt);
  dprintf(1,"%s\n",fmt4);
  if (hasvalue("label"))
    label = getparam("label");
  else
    label = infile;
    
  instr = stropen(infile, "r");
  while ((tags = list_tags(instr)) != NULL) {
    accum_item(*tags);
    free(*tags);
    free((char *)tags);
  }
  
  sprintf(current,fmt4,
	 mean_moment(&m), sigma_moment(&m),
	 min_moment(&m), max_moment(&m),
	 n_moment(&m));

  if (hasvalue("test")) {
    if (streq(current,test))
      printf("BSF %s: %s OK\n",label,current);
    else {
      printf("BSF %s: %s FAIL\n",label,current);
      printf("BSF %s: %s expected\n",label,test);
    }
  } else {
    printf("BSF %s: %s\n",label,current);
  }
}

void accum_item(string tag)
{
  string type, *tags, *tp;
  int *dims;

  type = get_type(instr, tag);
  if (streq(type, SetType)) {
    get_set(instr, tag);

    tags = list_tags(instr);
    for (tp = tags; *tp != NULL; tp++)
      accum_item(*tp);
    get_tes(instr, tag);

    for (tp = tags; *tp != NULL; tp++)
      free(*tp);
    free((char *)tags);
  } else {
    dims = get_dims(instr, tag);

    accum_data(tag, type, dims);
    if (dims != NULL)
      free((char *)dims);
  }
  free(type);
}

void accum_data(string tag, string type, int *dims)
{
  size_t dlen;
  byte *dat, *dp;
  real rv;

  dlen = get_dlen(instr, tag);
  dat = (byte*) allocate(dlen+1);
  get_data_sub(instr, tag, type, dat, dims, FALSE);
  
  for (dp = dat; dp < dat + dlen; ) {	/* loop over data array     */
    if (streq(type, AnyType)) {             /*   generic data?   */
      dp += sizeof(byte);
    } else if (streq(type, CharType)) {     /*   readable chars? */
      dp += sizeof(char);
    } else if (streq(type, ByteType)) {     /*   bytes of data?  */
      dp += sizeof(byte);
    } else if (streq(type, ShortType)) {    /*   short integers? */
      dp += sizeof(short);
    } else if (streq(type, IntType)) {      /*   standard ints?  */
      dp += sizeof(int);
    } else if (streq(type, LongType)) {     /*   long integers?  */
      dp += sizeof(long);
    } else if (streq(type, HalfpType)) {    /*   short int? */
      dp += sizeof(short);
    } else if (streq(type, FloatType)) {    /*   floating point? */
      rv = (real) *((float *) dp);
      dp += sizeof(float);
      if (streq(tag,ignore))
	  dprintf(1,"Ignoring %s = %g\n",tag,rv);
      else
	  accum_moment(&m,rv,1.0);
    } else if (streq(type, DoubleType)) {   /*   double numbers? */
      rv = (real) *((double *) dp);      
      if (streq(tag,ignore))
	  dprintf(1,"Ignoring %s = %g\n",tag,rv);
      else
	  accum_moment(&m,rv,1.0);
      dp += sizeof(double);
    } else
      error("accum_data: type %s unknown", type);
  }
  free((char *)dat); 
}


