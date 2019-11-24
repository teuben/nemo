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
  "in=???\n              input file name ",
  "VERSION=0.1\n         24-nov-2019 PJT ",
  NULL,
};

string usage="benchmark stats on all floating point values in a (binary) structured file";

string cvsid="$Id$";

local stream instr;	
local Moment m;

/* local functions */
void   accum_item   (string);
void   accum_data   (string, string, int *);
string find_name    (string);

#define BUFLEN   512


void nemo_main()
{
  string *tags;

  dprintf(2,"TSF: MaxSetLen = %d\n",MaxSetLen);
  ini_moment(&m,2,0);
  
  instr = stropen(getparam("in"), "r");
  while ((tags = list_tags(instr)) != NULL) {
    accum_item(*tags);
    free(*tags);
    free((char *)tags);
  }
  printf("%g %g %g %g %d\n",
	 mean_moment(&m), sigma_moment(&m),
	 min_moment(&m), max_moment(&m),
	 n_moment(&m));
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

  dlen = get_dlen(instr, tag);
  dat = (byte*) allocate(dlen+1);
  get_data_sub(instr, tag, type, dat, dims, FALSE);
  
  for (dp = dat; dp < dat + dlen; ) {		/* loop over data array     */
    if (streq(type, AnyType)) {             /*   output generic data?   */
      dp += sizeof(byte);
    } else if (streq(type, CharType)) {     /*   output readable chars? */
      dp += sizeof(char);
    } else if (streq(type, ByteType)) {     /*   output bytes of data?  */
      dp += sizeof(byte);
    } else if (streq(type, ShortType)) {    /*   output short integers? */
      dp += sizeof(short);
    } else if (streq(type, IntType)) {      /*   output standard ints?  */
      dp += sizeof(int);
    } else if (streq(type, LongType)) {     /*   output long integers?  */
      dp += sizeof(long);
    } else if (streq(type, HalfpType)) {	/*   output short int? */
      dp += sizeof(short);
    } else if (streq(type, FloatType)) {	/*   output floating point? */
      dp += sizeof(float);
    } else if (streq(type, DoubleType)) {	/*   output double numbers? */
      accum_moment(&m, *((double *) dp),1.0);
      dp += sizeof(double);
    } else
      error("accum_data: type %s unknown", type);
  }
  free((char *)dat); 
}


