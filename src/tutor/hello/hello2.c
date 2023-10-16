/*
 *  HELLO2:  Hello World, probably the most simple NEMO main(),
 *           without any extraneous nemo definitions.
 *
 */

#include <nemo.h>

void nemo_main() 
{
    // hello world as we always do
    printf("Hello2.\n");
    // NEMO's debug level printf(); use debug=1 or higher to see this (and more)
    dprintf(1,"Hello debug=1 World\n");
}
