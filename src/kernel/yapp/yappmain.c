/*
 *  YAPPMAIN:  TOOLKIT version for yapp, mainly to allow plotting a layout file
 *             or just to test if it works.
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
  "pgplot=\n      Query PGPLOT devices",
  "VERSION=0.3a\n 27-oct-2022 PJT",
  NULL,
};

local plcommand *layout;

void nemo_main()
{
    string name;

    if (hasvalue("pgplot")) {
      warning("pgplot= not implemented yet");
    }

    name = getparam("name");
    if (hasvalue("layout"))
      layout = pl_fread(getparam("layout"));
    else
      layout = NULL;

    plinit(name, 0.0, 20.0, 0.0, 20.0);     /* open device */
    if (layout)
      pl_exec(layout);
    else {
      warning("No layout= given, just yapping");
      pljust(0);
      pltext("yapp yapp", 10.0, 10.0, 0.5, 0.0);
    }
    plstop();
}


