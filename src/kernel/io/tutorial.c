/*
 * tutorial.c :  level dependent showing a file.
 *
 *    10-jul-2025   trying out some ideas
 */

#include <stdinc.h>
#include <stdarg.h>
#include <getparam.h>

int tutorial_level=0;	/* needs to be global; see also getparam.c */


void nemo_set_tutorial(int tutorial)
{
  tutorial_level = tutorial;
}

local string red   = "\033[0;31m";
local string green = "\033[0;32m";
local string nc    = "\033[0m";

void tutorial(int level, string name)
{
  char cmd[128];
  
  fprintf(stderr,"%sTUTORIAL(%s) [%s]%s\n", red,name, getargv0(),nc);
  sprintf(cmd,"cat $NEMO/docs/tutorial/%s.tut  1>&2", name);
  system(cmd); 
}

