#include <nemo.h>

string defv[] = {
    "saveparam=\n       List of parameters to save values"
    "VERSION=1\n        18-dec-99 pjt",
    NULL
};

string usage = "test";

nemo_main()
{
    string *wp = burststring(getparam("saveparam"),", ");
    string *sp;

    for (sp=wp; *sp; sp++) {
        dprintf(0,"Saving %s\n",*sp);
    }
        
}

void writeparam(string key, string value)
{

}
