/* 
 * CCDRGB: read an R, G and B image, and create composite
 *		** unfinished **
 *
 *	 8-jun-98 	written	- very quick and dirty      pjt
 *                      
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
        "in1=???\n       Input image(s)",
        "in2=???\n       Input image(s)",
	"out=???\n      Output image file",
	"c=0.5\n	color composition index RGB=(A, B+c*A, B)",
	"range1=\n      Range in image 1",
	"range2=\n      Range in image 1",
	"VERSION=1.0\n  8-jun-98 PJT",
	NULL,
};

string usage = "create composite RGB color images from multiple input images";


void nemo_main(void)
{
    stream   instr1, instr2, outstr;
    int      nx, ny, n1, n2, i, j;
    imageptr iptr1=NULL, iptr2=NULL, optr=NULL;      /* pointer to images */
    real     a, b, o, c = getdparam("c");
    real     range1[2], range2[2], scale1, offset1, scale2, offset2;

    instr1 = stropen(getparam("in1"), "r");
    instr2 = stropen(getparam("in2"), "r");
    n1 = nemoinpd(getparam("range1"),range1,2);
    n2 = nemoinpd(getparam("range1"),range2,2);
    if (n1==0) {
        scale1 = 1.0;
        offset1 = 0.0;
    } else if (n1==2) {
        scale1 = 1.0/(range1[1]-range1[0]);
        offset1 = range1[0];
    } else
        error("range1= needs two values");

    if (n2==0) {
        scale2 = 1.0;
        offset2 = 0.0;
    } else if (n2==2) {
        scale2 = 1.0/(range2[1]-range2[0]);
        offset2 = range2[0];
    } else
        error("range2= needs two values");


    read_image( instr1, &iptr1);
    read_image( instr2, &iptr2);

    outstr = stropen(getparam("out"), "w");

    for (j=0; j<ny; j++) {                  /* make copy of image */
    	for (i=0; i<nx; i++) {
            a = MapValue(iptr1,i,j);
            b = MapValue(iptr2,i,j);

            a = scale1*a + offset1;
            b = scale2*b + offset2;

	}
    }
    write_image(outstr, iptr1);
}


