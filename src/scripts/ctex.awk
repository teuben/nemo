#! /usr/bin/awk -f
#
#  extract CTEX comments from C source code...
#  Each paragraph should be structured as follows:
#
#	/*CTEX
#	 *  bla bla bla (tex source)
#	 */
#
#  Usage:
#	ctex.awk file.c > file.tex
#
#  See also the 'ctex' shell scripts, which combines various ctex.awk
#  calls into one (la)tex file

BEGIN {
  on = 0;
  printf("%% File: %s\n",ARGV[1]);
}

{
  if ($1 == "*" && on) {
    for (i=2; i<=NF; i++)
        printf("%s ",$i);
    printf("\n");
  } else if ($1 == "/*CTEX") {
    printf("%% CTEX Line: %d\n",NR);
    on = 1;
  } else if ($1 == "*/") {
    printf("\n");
    on=0;
  }
}  
