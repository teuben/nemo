/*
 *
 * format: little TESTBED program with which you can figure
 *	   out your endianism and if you're float/double is
 *	   in IEEE format
 * big endian: (sun)
 *	double    : 3F F0 00 00 00 00 00 00
 *	float     : 3F 80 00 00
 *	int       : 00 00 00 01
 *
 *	mar-94	created			peter teuben
 *      jun-94  default file=           pjt
 */



#include <stdinc.h>
#include <getparam.h>

string defv[]= {
    "file=/tmp/junk\n   File where I/O will take place",
    "double=1.0\n  	Double value to write",
    "float=1.0\n        Float value to write",
    "int=1\n            Integer value to write",
    "VERSION=1.0a\n     10-jun-94 PJT",
    NULL,
};

string usage="bit pattern tester for FITS format conversion";

void nemo_main(void)
{
    string fname = getparam("file");
    stream str = stropen(fname,"w!");
    float fval = getdparam("float");
    double dval = getdparam("double");
    int ival = getiparam("int");
    unsigned char dpat[32], fpat[32], ipat[32];
    int nd, nf, ni;

    nd = sizeof(double);
    nf = sizeof(float);
    ni = sizeof(int);

    printf("double=%lf float=%f int=%d valued\n",dval,fval,ival);
    printf("double=%d float=%d int=%d bytes\n",nd,nf,ni);

    if (fwrite(&dval,nd,1,str) != 1) error("writing double");
    if (fwrite(&fval,nf,1,str) != 1) error("writing float");
    if (fwrite(&ival,ni,1,str) != 1) error("writing int");
    strclose(str);

    str = stropen(fname,"r");
    if (fread(dpat,1,nd,str) != nd) error("reading double");
    if (fread(fpat,1,nf,str) != nf) error("reading float");
    if (fread(ipat,1,ni,str) != ni) error("reading int");
    strclose(str);

    show("double",dpat,nd);
    show("float",fpat,nf);
    show("int",ipat,ni);
}


show(name,pat,n)
string name;
unsigned char *pat;
int n;
{
    int i;

    printf("%-10s: ",name);
    for (i=0; i<n; i++)
        printf("%02X ",pat[i]);
    printf("\n");
}    
