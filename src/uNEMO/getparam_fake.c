/*
 *   getparam_fake:     two functions needed by uNEMO when the real
 *			getparam is not available.
 *
 */
#include <stdinc.h>

static string _fake = "WARNING: getparam_fake used";

string getparam(string key)
{
    return key;
}


void finiparam(void)
{
}
