/*
 * Full test of all of Micro Nemo's capabilities
 *	4-oct-90	Created for new 'NEMO'		PJT
 *	27-feb-94	ansi
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>

string defv[] = {
    "in=\n         Example input file",
    "out=\n        Example output file",
    "nesting=1\n   Nesting of i/o example",
    "nstack=1\n    Stacking data items",
    "VERSION=1.0a\n  27-feb-94 PJT",
    NULL,
};
string usage="testing filestruct routines";

nemo_main()
{
    string name;
    char fname[32];
    int n,i;
    stream str;

    name = getparam("in");
    if (*name)
        get_file(name);
    else
        dprintf(0,"No input file specified\n");
    
    name = getparam("out");
    if (*name)
        put_file(name,getiparam("nesting"),getiparam("nstack"));
    else
        dprintf(0,"No output file specified");

    do_stropen();
}

/*  MAXSTREAM is max. number of open files to try */
#define MAXSTREAM 512

local int break1=0;
local stream str[MAXSTREAM];

do_stropen()
{
    int i, catch1();
    char fname[32];

    dprintf(0,"Going to try how many open files\n");
    recover(catch1);        /* since we know it will fail - catch it */
    for(i=0;;i++) {
        sprintf(fname,"junk_%d",i);
        dprintf(1,"Opening file %s\n",fname);
        str[i] = stropen(fname,"w");
        if (str[i]==NULL || break1) {
            str[i] = NULL;
            break;
        }
        put_history(str[i]);
    }
    printf("%d files could be opened simultaneously\n",i);
    recover(NULL);          /* reset error recover */
    for(i=0;;i++) {         /* close all open files */
        if (str[i]==NULL) break;
        sprintf(fname,"junk_%d",i);
        dprintf(1,"Closing file %s\n",fname);
        strclose(str[i]);
    }
    for(i=0;;i++) {         /* delete all those junk files */
        sprintf(fname,"junk_%d",i);
        if (file_size(fname) < 0) break;
        dprintf(1,"Deleting file %s\n",fname);
        if (unlink(fname)<0) error("Error deleting file %s\n",fname);
    }
}

catch1()
{
    printf("### Recovered from fatal file table overflow: ###\n");
    break1 = 1;
}

get_file(name)
string name;
{
    stream str;
    string *names, nm;
    char tag[32];
    int i;

    str = stropen(name,"r");
    get_history(str);
    if (!get_tag_ok(str,"FStest")) {
        names = list_tags(str);
        printf("ITEMS: ");
        for (i=0; names[i]!=NULL; i++)
                printf("%s ",names[i]);
        printf("\n");
        error("%s: Not a valid input file - must be FStest",name);
    }
    get_set(str,"FStest");
    get_tes(str,"FStest");
    if (!get_tag_ok(str,"STACK")) {
        get_set(str,"STACK");
        for(;;){
            sprintf(tag,"I_%d",i);
#if 1
            if (!get_tag_ok(str,tag)) break;        /* search */
#else
            get_data(str,tag,IntType,&i,0);
            printf("Reading tag %s with data %d\n",tag,i);
#endif
        }
        get_tes(str,"STACK");
    }
    strclose(str);
}

put_file(name,n,m)
string name;
int n;
int m;
{
    stream str;
    int i, catch1();
    char setname[16];
    

    str = stropen(name,"w");
    put_history(str);
    put_set(str,"FStest");
    for (i=0; i<n; i++) {
        sprintf(setname,"ID_%d",i);
        put_set(str,setname);
        put_data(str,"i",IntType,&i,0);
    }
    for (i=n-1; i>=0; i--) {
        sprintf(setname,"ID_%d",i);
        put_tes(str,setname);
    }
    put_tes(str,"FStest");
    put_set(str,"STACK");
    for (i=0;i<m;i++) {
        sprintf(setname,"I_%d",i);
        put_data(str,setname,IntType,&i,0);
    }
    put_tes(str,"STACK");
    strclose(str);
}

