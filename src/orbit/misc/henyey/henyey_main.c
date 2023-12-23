/* THIS FILE HAS BEEN CREATED BY ftoc - do not edit */

#include <stdinc.h>

string defv[] = {
    "x0=1\n                  Initial x coordinate",
    "y0=0\n                  Initial y coordinate",
    "u0=0\n                  Initial u velocity",
    "v0=1\n                  Initial v velocity",
    "per=3\n                 Period",
    "type=1\n                Type of orbit {1,2,3}",
    "ome=0.0\n               Pattern speed of potential",
    "norbit=1\n              Number of orbits to compute",
    "step=0.1\n              Step to increase non-zero position",
    "VERSION=1.0\n           18-sep-91 PJT",
    NULL,
};

string usage="NEMO program with unknown intent";

nemo_main()
{
    nemomain_();
}
