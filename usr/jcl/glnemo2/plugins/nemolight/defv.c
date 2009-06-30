/*
 * defv.c: dummy definition of the required "extern string defv[]"
 *         This is to ensure that programmers who forgot to declare it,
 *         deliberate or not, will have a linkable program.
 *
 */

#include <stdinc.h>

string defv[] = {
    "VERSION=0.0\n 26-feb-92 PJT  ### NEMO program without keywords (defv)",
    NULL,
};
