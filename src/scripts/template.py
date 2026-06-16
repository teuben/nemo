#! /usr/bin/env python3
import os
from datetime import datetime
import sys


if len(sys.argv) == 1:
    print(f"Usage: {sys.argv[0]} progname key1=def1 ...")
    print("version 24-jan-2012")
    print("Write template NEMO program based on keywords and default values")
    print("Output will be in a file progname.c, which cannot be overwritten")
    exit(0)

prog = sys.argv[1]
if ("." in prog):
    name = prog
else:
    name = f"{prog}.c"

date = datetime.now().strftime("%d-%b-%Y")

if os.path.exists(name):
    print(f"File {name} already exists, will not overwrite it.")
    exit(0)

with open(name, "w") as f:
    f.write("/*\n")
    f.write(" *  TEMPLATE: created by NEMO \n")
    f.write(" */\n\n")
    f.write("#include <nemo.h>\n\n")
    f.write("string defv[] = {\n")

    for arg in sys.argv[2:]:
        key, value = arg.split("=")
        f.write(f'    "{arg}\\n    some help",\n')

    f.write(f'    "VERSION=0.0\\n       {date} XYZ",\n')
    f.write("    NULL,\n")
    f.write("};\n\n")
    f.write("string usage=\"TEMPLATE program usage -- fill in yourself\";\n\n")
    f.write("void nemo_main()\n")
    f.write("{\n")
    f.write("    /* your code goes here */\n")
    f.write("    /* e.g.   string s = getparam(\"s\");  */\n")
    f.write("    /* e.g.   int    i = getiparam(\"i\"); */\n")
    f.write("    /* e.g.   real   r = getrparam(\"r\"); */\n")
    f.write("    /* e.g.   double d = getdparam(\"d\"); */\n")
    f.write("    /* followed by lots of lovely calculations */\n")
    f.write("    warning(\"New template NEMO program\");\n")
    f.write("    /*\n")

    for arg in sys.argv[2:]:
        key, _ = arg.split("=")
        f.write(f"    string {key} = getparam(\"{key}\");\n")

    f.write("     */\n")
    f.write("    /*\n")
    f.write("     *  other useful things:\n")
    f.write("     *    if (hasvalue(\"xxx\") && !hasvalue(\"yyy\"))\n")
    f.write("     *        error(\"need value for yyy\");\n")
    f.write("     */\n")
    f.write("}\n")

print("All done, {name} created. You can compile it using the command")
print("        bake {prog}")
print("  or:   mknemo {prog}")
print("When done, don't forget to create a template man page:")
print(f"        $NEMO/src/scripts/mkman {prog} > $NEMO/man/man1/{prog}.1")