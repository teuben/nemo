/*
 * yancpars:	add YANC parameters to a snapshot, so YANC has (some
 *		better) values for its integration parameters
 *              Not sure if this will be the official way to setup YANC runs
 * 
 *	26-may-2001	Written in Mexico, for Walter Dehnen's YANC code  PJT
 */

#if 0
/* need to fix getrparam problem !!! */
#define SINGLEPREC
#endif

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>

string defv[] = {
    "in=???\n		Input file",
    "out=???\n		Output file",
    "eps=0.05\n		Softening length",
    "kernel=21\n        Softening kernel #",
    "theta=0.5\n        Treecode opening angle",
    "hmin=5\n           Timestep dt=2^-hmin",
    "convert=d2f\n      Conversion options",
    "VERSION=1.0\n	26-may-01 PJT",
    NULL,
};

string usage = "Copy a binary structured file";


extern string *burststring(string,string);

nemo_main()
{
    stream outstr, instr = stropen(getparam("in"),"r");
    int    i;
    string *tags;
    real eps = getrparam("eps");
    real theta = getrparam("theta");
    real hmin = getrparam("hmin");
    int kernel = getiparam("kernel");
    string *cvt = burststring(getparam("convert"),", ");

    dprintf(0,"eps=%g theta=%g hmin=%g kernel=%d\n",eps,theta,hmin,kernel);
    
    get_history(instr);
    if (!get_tag_ok(instr,SnapShotTag))
        error("Input file is not a valid snapshot");

    outstr = stropen(getparam("out"),"w");
    put_history(outstr);

    get_set(instr,SnapShotTag);
    put_set(outstr,SnapShotTag);

    tags = list_tags(instr);

    for (i=0; tags[i]; i++)  {                 /* copy old stuff */
      if (streq(tags[i],YancParametersTag)) continue;
      dprintf(1,"Copying %s\n",tags[i]);
#if 0
      copy_item_cvt(outstr,instr,tags[i],cvt);
#else
      copy_item(outstr,instr,tags[i]);
#endif
    }
    dprintf(1,"Done\n");
    free((char *)*tags);
    free((char *)tags);
    
    /* add new YANC parameters */
    put_set(outstr,YancParametersTag);
    put_data(outstr, YancEpsTag, RealType, &eps, 0);
    put_data(outstr, YancThetaTag, RealType, &theta, 0);
    put_data(outstr, YancHminTag, RealType, &hmin, 0);
    put_data(outstr, YancKernelTag, IntType, &kernel, 0);
    put_tes(outstr, YancParametersTag);
    
    get_tes(instr,SnapShotTag);
    put_tes(outstr,SnapShotTag);
    
    strclose(instr);
    strclose(outstr);
}
