/*
 *  YAPPMAIN:  TOOLKIT version for yapp
 *
 *       
 *
 *  V1.1 6-jan-2011     created
 */

#include <stdinc.h>
#include <getparam.h>
#include <yapp.h>
#include <layout.h>

string defv[] = {
  "name=***\n     Override for yapp= system keyword (*** = default)",
  "layout=\n      Optional layout file",
  "VERSION=0.1\n  6-dec-2011 PJT",
  NULL,
};

local plcommand *layout;

nemo_main()
{
    int i, j, ip, np, nc;
    string name, dumpfile, headline;
    char label[80];

    name = getparam("name");
    if (hasvalue("layout"))
      layout = pl_fread(getparam("layout"));
    else
      layout = NULL;

    plinit(name, 0.0, 20.0, 0.0, 20.0);     /* open device */
    if (layout) pl_exec(layout);
    plstop();
}


