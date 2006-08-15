/*
 * yancpars:	add YANC parameters to a snapshot, so YANC has (some
 *		better) values for its integration parameters
 *
 *		now deprecated, since YancNemo (now called gyrfalcON) is now the official way
 *		but we keep this code since the method is rather unique
 * 
 *	26-may-2001	Written in Mexico, for Walter Dehnen's YANC code  PJT
 *      30-jun-2001	1.0a now declared it deprecated
 *      15-aug-2006     1.0b document gyrfalcON and the fact that it does not compile anymore
 *                           due to missing things like YancEpsTag etc.
 */

#if 0
/* need to fix getrparam problem !!! */
#define SINGLEPREC
#endif

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>

string defv[] = {
    "in=???\n		Input file",
    "out=???\n		Output file",
    "eps=0.05\n		Softening length",
    "kernel=21\n        Softening kernel #",
    "theta=0.5\n        Treecode opening angle",
    "hmin=5\n           Timestep dt=2^-hmin",
    "convert=d2f\n      Conversion options",
    "VERSION=1.0b\n	30-jun-01 PJT",
    NULL,
};

string usage = "Add (yanc) parameter to a snapshot";


extern string *burststring(string,string);

void nemo_main()
{
    stream outstr, instr = stropen(getparam("in"),"r");
    int    i;
    string *tags;
    real eps = getrparam("eps");
    real theta = getrparam("theta");
    real hmin = getrparam("hmin");
    int kernel = getiparam("kernel");
    string *cvt = burststring(getparam("convert"),", ");

    dprintf(1,"eps=%g theta=%g hmin=%g kernel=%d\n",eps,theta,hmin,kernel);
    dprintf(1,"using conversion %s\n",cvt[0]);
    
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
