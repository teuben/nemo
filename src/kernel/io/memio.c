/*
 * MEMIO: safe parsing of misaligned data
 *
 *	10-may-94	Created			Peter Teuben
 *       3-jul-94       added memcread & memsread        PJT
 *       7-mar-96       added optional binary output     PJT
 */

#include <stdinc.h>
#include <memio.h>


mem *memopen(char *buf, int buflen)
{
    mem *mp = (mem *) allocate(sizeof(mem));
    
    dprintf(2,"memopen: %d bytes\n",buflen);
    mp->buf = mp->bp =  buf;
    mp->buflen = buflen;
    return mp;
}

void memclose(mem *mp)
{
    dprintf(2,"memclose: read %d bytes\n",mp->bp-mp->buf);
    free((char *)mp);
}

float memfread(mem *mp)
{
    float x;
    char *cp = (char *) &x;
    int i;

    if (mp->bp - mp->buf + sizeof(float) > mp->buflen)
        error("memfread: not enuf memory to read float; buflen=%d",mp->buflen);
    for (i=0; i<sizeof(float); i++)
        *cp++ = *mp->bp++;
    return x;
}

double memdread(mem *mp)
{
    double x;
    char *cp = (char *) &x;
    int i;

    if (mp->bp - mp->buf + sizeof(double) > mp->buflen)
        error("memdread: not enuf memory to read double; buflen=%d",mp->buflen);
    for (i=0; i<sizeof(double); i++)
        *cp++ = *mp->bp++;
    return x;
}

int memiread(mem *mp)
{
    int x;
    char *cp = (char *) &x;
    int i;

    if (mp->bp - mp->buf + sizeof(int) > mp->buflen)
        error("memiread: not enuf memory to read int; buflen=%d",mp->buflen);
    for (i=0; i<sizeof(int); i++)
        *cp++ = *mp->bp++;
    return x;
}

short memsread(mem *mp)
{
    short x;
    char *cp = (char *) &x;
    int i;

    if (mp->bp - mp->buf + sizeof(short) > mp->buflen)
        error("memsread: not enuf memory to read int; buflen=%d",mp->buflen);
    for (i=0; i<sizeof(short); i++)
        *cp++ = *mp->bp++;
    return x;
}

char memcread(mem *mp)
{
    char x;

    if (mp->bp - mp->buf + sizeof(char) > mp->buflen)
        error("memcread: not enuf memory to read int; buflen=%d",mp->buflen);
    x = *mp->bp++;
    return x;
}

/* memseek:
 * 	 mode:  0=from start  1=from current  2=from end
 */	 

void memseek(mem *mp, int loc, int mode)
{
    if (mode==0)
        mp->bp = mp->buf + loc;
    else if (mode==1)
        mp->bp += loc;
    else if (mode==2)
        mp->bp = mp->buf + mp->buflen - loc;
}


#ifdef TOOLBOX
   

#include <getparam.h>

string defv[] = {
    "value=\n       enter value if to parse",
    "itype=i\n      type of input value (i,f,d)",
    "otype=i\n      type of output value (i,f,d)",
    "fmt=\n         format used in output value",
    "out=\n         Optional output",
    "VERSION=1.1\n  7-mar-96 PJT", 
    NULL,
};

string usage = "misaligned data/memory routines with optional disk I/O";

nemo_main()
{
    stream ostr;
    string itype = getparam("itype");
    string otype = getparam("otype");
    string fmt = getparam("fmt");
    char buf[64], format[32];
    float fval;
    double dval;
    int ival;
    mem *mp;

    if (hasvalue("value")) {
        switch(*itype) {
        case 'f':
            fval = getdparam("value");
            memcpy(buf,(char *)&fval,sizeof(float));
            break;
        case 'd':
            dval = getdparam("value");
            memcpy(buf,(char *)&dval,sizeof(double));
            break;
        case 'i':
            ival = getiparam("value");
            memcpy(buf,(char *)&ival,sizeof(int));
            break;
        default:
            error("invalid input type %s",itype);
        } 
    }


    if (hasvalue("fmt")) sprintf(format,"%s\n",getparam("fmt"));

    if (hasvalue("out"))
        ostr = stropen(getparam("out"),"w");
    else
        ostr = NULL;


    mp = memopen(buf,32);
    switch(*otype) {
        case 'f':
            if (ostr) fwrite(buf,sizeof(float),1,ostr);
            fval = memfread(mp);
            if (!hasvalue("fmt")) sprintf(format,"%s\n","%g");
            printf(format,fval);
            break;
        case 'd':
            if (ostr) fwrite(buf,sizeof(double),1,ostr);
            dval = memdread(mp);
            if (!hasvalue("fmt")) sprintf(format,"%s\n","%g");
            printf(format,dval);
            break;
        case 'i':
            if (ostr) fwrite(buf,sizeof(int),1,ostr);
            ival = memiread(mp);
            if (!hasvalue("fmt")) sprintf(format,"%s\n","%d");
            printf(format,ival);
            break;
        default:
            error("invalid output type %s",otype);
    }

    if (ostr) strclose(ostr);



}


#endif
