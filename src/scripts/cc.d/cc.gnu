#! /bin/csh -f
#  cc script for generic GNU cc compiler
#  unfortunately the -static switch causes the compiler to
#  ignore sharable libraries, and hence create larger executables -
#  but is needed for loadobj. If you don't use loadobj, -static
#  can be taken out.
#  Note: for sun4 or sun3 an extra -Dsun3 or -Dsun4 is recommended
#	 you also need -Dbsd if symbols appear prepended with _ (as is on sun3/4)
#	 But e.g. on solaris and linux, this is not needed

# June 1998:
#  need -DSOLARIS on sun5's
#  took out -static, added -ldl, since that appears to work now

exec gcc -I$NEMOINC -L$NEMOLIB -ansi -fwritable-strings -Wall $*
