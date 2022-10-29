#
# command line parameter parsing tools
#       28-oct-2022   finally a good commandline parser
#
#
# things to think about:
#   - more pythonic ?     p.get("a")   vs.     p["a"]
#                         p.has("a")   vs.     "a" in p

import os
import sys

class Param(object):
    """
    Simple NEMO style parameter parsing,
    - still limited error checking
    - only program keywords, though -h and --help give reminders

    """
    def __init__(self, keyval, usage, debug=0):
        self.keyval = keyval
        self.k   = []
        self.key = {}
        self.help = {}
        for a in keyval:
            w=a.split('\n')
            kv=w[0].split('=')
            self.k.append(kv[0])
            self.key[kv[0]] = kv[1]
            self.help[kv[0]] = w[1].lstrip()
        if debug>0: print(self.key)
        self.name = sys.argv[0]
        self.usage = usage
        equals = False
        na = 0
        for a in sys.argv[1:]:
            if a=='-h':
                for i in range(len(self.k)):
                    print("%s=%s" % (self.k[i], self.key[self.k[i]]))
                sys.exit(0)
            if a=='--help':
                for i in range(len(self.k)):
                    print("%-10s: %s [%s]" % (self.k[i], self.help[self.k[i]], self.key[self.k[i]]))
                print(self.usage)
                sys.exit(0)
            if a=='--version':
                print("%s  %s  %s"  % (self.name, self.key["VERSION"], self.help["VERSION"]))
                sys.exit(0)            
            na = na + 1
            if debug>0: print("processing ",a)
            if a.find('=') < 0:
                if na > len(self.k):
                    print("too many args")                
                    sys.exit(0)
                if equals == True:
                    print("bad order, no = found")
                    sys.exit(0)                    
                keyb = self.k[na-1]
                self.key[keyb] = a
            else:
                equals = True
                kv = a.split('=')
                self.key[kv[0]] = kv[1]
            
    def get(self, keyword):
        return self.key[keyword]

    def has(self, keyword):
        if len(self.key[keyword]) > 0:
            return True
        return False

    def listf(self, keyword):
        return [float(x) for x in self.get(keyword).split(',')]

    def listi(self, keyword):
        return [int(x) for x in self.get(keyword).split(',')]



if __name__ == "__main__":
    keyval = [
        "a=1\n       help for a",
        "b=2\n       help for b",
        "c=3\n       help for c",
    ]

    p = Param(keyval)
    print('a:',p.get("a"))
    print('b:',p.get("b"))
    print('c:',p.get("c"))
