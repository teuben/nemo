#! /usr/bin/env python
#
#  Extract the TLDR markdown information from a man page
#  Typical usage:
#       man2tldr.py tsf.1 > ../tldr/tsf.md
#  where tsf.md will be 0-length if there is no TLDR section
#
#  See also discussion in https://github.com/tldr-pages/tldr/issues/6154
#
# 21-jun-2021  quick hack - PJT 
#


import os, sys

def my_match(word1, word2):
    """   when is word2 equal or inside of word1
          (e.g. word1 could contain quotes)
    """
    if word1.find(word2) < 0:
        return False
    return True


def checkMan(file):
    """   extract the TLDR section of the man page
          formatted in markdown style
    """
    have_name = False
    have_tldr = False
    lines = open(file).readlines()
    for line in lines:
        line = line.strip()
        if line[:3] == ".SH":
            words = line.split()
            if my_match(words[1],"NAME"):
                have_name = True
                # print(line)
            if have_tldr:
                have_tldr = False
            if my_match(words[1],"TLDR"):
                have_tldr = True                
                # print(line)
        elif have_name:
            words = line.split()
            name = words[0]
            desc = words[2:]
            have_name = False
            print_name = False
        elif have_tldr:
            # first write the header
            if not print_name:
                print("# %s" % name)
                print(" ")
                print("> %s" % ' '.join(desc))
                print(" ")
                print_name = True
            # then the TLDR section
            # this needs an alternate re-formatting, TBD. Just for show now
            # alternatively, the TLDR is already in the expected format 
            print(line)



if __name__ == "__main__":
    if len(sys.argv) > 1:
        checkMan(sys.argv[1])    
