/*
 *    JD2ISO:   convert a julian date to ISO civil time
 */

#include <nemo.h>
#include <time.h>

string defv[] = {
  "jd=\n          Julian Day",
  "mjd=\n         Modified Julian Day (MDJ=JD-2,400,000.5)",
  "tjd=\n         Truncated Julian Day (TjD=JD-2,440,000.5)",
  "VERSION=0.1\n  8-dec-2022 PJT",
  NULL,
};

string usage="Convert (Modified) Julian Day to a (ISO) civil time";


local real jd1970 = 2440587.50000;
local real mjd0   = 2400000.5;
local real tjd0   = 2440000.5;
  

local void jd2iso(real jd)
{
  time_t tnow;
  time(&tnow);
  dprintf(1,"now=%ld\n", tnow);
  dprintf(1,"jd=%f\n",jd);
  dprintf(1,"mjd=%f\n",jd-mjd0);

  tnow = (time_t) ((jd-jd1970)*24*3600);
  dprintf(1,"sec-1970=%ld\n",tnow);
  if (tnow < 0) error("cannot parse negative seconds: %ld ~ %g years",tnow, tnow/(365.35*24*3600));
  printf("%s", ctime(&tnow));
  // @todo timezone is local, not UT
}

void nemo_main()
{
  real jd = getrparam("jd");
  real mjd = getrparam("mjd");
  real tjd = getrparam("tjd");

  if (hasvalue("jd"))
    jd2iso(jd);
  if (hasvalue("mjd"))
    jd2iso(mjd + mjd0);
  if (hasvalue("tjd"))
    jd2iso(tjd + tjd0);
  
}
