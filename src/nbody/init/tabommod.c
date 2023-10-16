/*
 * TABOMMOD: tabulate an Osipkov-Merritt model
 *
 *     5-jun-2023   created - PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

string defv[] = {
  "in=???\n			  input OM model",
  "VERSION=0.1\n		  5-jun-2023 PJT",
  NULL,
};

string usage =  "tabulate an Osipkov-Merritt model";

/*
 * The structure of the model to generate is specified by the following.
 */

#define MXTB 1024		/* max number of values per table             */

local real anisorad;		/* anisotropy radius (if le 0, isotropic)     */
local int ntab;			/* count of tabulated values                  */

local real mtab[MXTB]; /* table of interior masses                   */
local real rtab[MXTB]; /* table of radii,          */
local real ptab[MXTB]; /* table of potentials,   */
local real ftab[MXTB]; /* table of distrib. func,  */


local void readmodel(string);

void nemo_main()
{
    readmodel(getparam("in"));			/* read tabulated model     */
}

/*
 * READMODEL: initalize tables for radius, enclosed mass, potential,
 * anisotropic "energy" and distribution function from a file.
 */

local void readmodel(string name)
{
    stream instr;
    int i;

    instr = stropen(name, "r");
    get_history(instr);
    get_set(instr, "OsipkovMerrittModel");
    get_data_coerced(instr, "AnisoRadius", RealType, &anisorad, 0);
    get_data(instr, "Ntab", IntType, &ntab, 0);
    if (ntab > MXTB) error("Not enough space for tables, MXTB=%d",MXTB);
    get_data_coerced(instr, "Radius", RealType, rtab, ntab, 0);
    get_data_coerced(instr, "Mass", RealType, mtab, ntab, 0);
    get_data_coerced(instr, "Potential", RealType, ptab, ntab, 0);
    get_data_coerced(instr, "DistFunc", RealType, ftab, ntab, 0);
    get_tes(instr, "OsipkovMerrittModel");
    strclose(instr);

    dprintf(1, "[readmode: ntab = %d]\n", ntab);
    printf("# Radius Mass Potential DistFunc\n");
    for (i=0; i<ntab; i++) {
      printf("%g %g %g %g\n", rtab[i], mtab[i], ptab[i], ftab[i]);
      
      if (i==0) continue;
      if (rtab[i-1] > rtab[i]) warning("Radius %d not increasing?",i+1);
      if (mtab[i-1] > mtab[i]) warning("Mass %d not increasing?",i+1);
      if (ptab[i-1] > ptab[i]) warning("Potential %d not increasing?",i+1);
      if (ftab[i-1] < ftab[i]) warning("DistFunc %d not decreasing?",i+1);
    }
}

