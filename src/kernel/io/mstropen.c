/*
 * MSTROPEN:   aid opening multiple files serially using a template
 *             useful if you want to split simulation files in
 *             pieces under your own control
 *
 *    8-may-2002  created on Amtrak,  Amar suggesting this for the 2GB problem   PJT
 *   21-may-2002  oops, template is a reserved word in C++ :-)                   pjt
 *   23-may-2002  added mode to mstr_open; though currently only write allowed
 *   14-sep-2002  added query what counter we're in                              pjt
 *    2-jul-2004  g++ safer
 */

#include <stdinc.h>

mstr *mstr_init(string tmplate)
{
  mstr *mp = (mstr *) allocate(sizeof(mstr));

  mp->tmplate = (char *) allocate(strlen(tmplate)+1);
  mp->filename = (char *) allocate(strlen(tmplate)+20);
  strcpy(mp->tmplate,tmplate);
  mp->status = 0;
  mp->count = 0;
  mp->mode = (strchr(tmplate,'%') != NULL);
  if (mp->mode) {
    int n;
    n = sprintf(mp->filename,mp->tmplate,mp->count);
    if (n < 0) error("Problem using %s as tmplate filename",mp->tmplate);
  } else
    strcpy(mp->filename,mp->tmplate);
  dprintf(1,"mstr_init: %s\n",mp->filename);

  mp->ostr = stropen(mp->filename,"w");
  return mp;
}

stream mstr_open(mstr *mp, string mode)
{
  if (mp->count == 0) {
    mp->count++;
    return mp->ostr;
  }
  if (mp->mode) {
    sprintf(mp->filename,mp->tmplate,mp->count);
    strclose(mp->ostr);
    mp->ostr = stropen(mp->filename,"w");
  }
  dprintf(1,"mstr_open: %s count=%d\n",mp->filename,mp->count);
  mp->count++;
  return mp->ostr;
}

void mstr_close(mstr *mp)
{
  dprintf(1,"mstr_close: %s count=%d\n",mp->tmplate,mp->count);
  if (mp->ostr) strclose(mp->ostr);
  free(mp->filename);
  free(mp->tmplate);
  free(mp);
}

int mstr_count(mstr *mp)
{
  return mp->count;
}

int mstr_multi(mstr *mp)
{
  return mp->mode;
}


#if defined(TESTBED)

#include <getparam.h>

string defv[] = {
  "out=test%03d.dat\n     template name",
  "n=10\n                 number to write",
  "VERSION=1.0a\n         21-may-2002 PJT",
  NULL,
};

string usage="testing mstropen";

nemo_main()
{
  mstr *mp;
  int i, n = getiparam("n");


  mp = mstr_init(getparam("out"));
  for (i=0; i<n; i++) {
    stream ostr = mstr_open(mp,"w");
    fprintf(ostr,"Hello world, i=%d\n",i);
  }
  mstr_close(mp);
}
#endif
