#! /usr/bin/env python
#
#   align columns in a file.   doesn't work well if number of columns
#   not the same through out the file
#

import os, sys

def readtable(file=sys.stdin):
    """ gobble the file
    """
    if file is not sys.stdin:
        fp = open(file,'r')
    else:
        fp = sys.stdin
    lines = fp.readlines()
    if file is not sys.stdin:
        fp.close()
    return lines


def process(lines):
    """ measure the file, and output 
    """
    # loop over all lines and find first line that has non-zero words
    wlen = 0
    for l in lines:
        w = l.split()
        if len(w) > 0:
            wlen = len(w)
            break

    # initialize a width[] list for each column
    width = list(range(wlen))
    for i in range(wlen):
        width[i] = 0

    # compute the max width for each column
    for l in lines:
        words = l.strip().split()
        for i in range(wlen):
            if len(words[i]) > width[i]: width[i] = len(words[i])

    # prepare a list with the proper format statement
    fmt = []
    for i in range(wlen):
        fmt.append('%%-%ds' % width[i])

    # finally loop over the table again, now printing out each
    # column with the correct width
    for l in lines:
        words = l.strip().split()
        s = ""
        for i in range(wlen):
            # print fmt[i], words[i]
            s = s + fmt[i] % words[i] + " "
        print(s.strip())

#
def process_dot(lines):
    """ measure the file, and output 
    """
    # loop over all lines and find first line that has non-zero words
    wlen = 0
    for l in lines:
        words = l.strip().split()
        if len(words) > 0:
            wlen = len(words)
            # initialize a width[] list for each column
            width = list(range(wlen))
            dot   = list(range(wlen))
            dot1  = list(range(wlen))
            for i in range(wlen):
                width[i] = 0
                dot[i] = words[i].find('.')
                dot1[i] = 0
            break

    #                    dot width
    #       aaaa.        4   5
    #         aa.aa      2   5
    #           .aaaa    0   5

    # compute the max width for each column
    for l in lines:
        words = l.strip().split()
        for i in range(wlen):
            if dot[i] < 0: 
                if len(words[i]) > width[i]: width[i] = len(words[i])
            else:
                d = words[i].find('.')
                if d<0:
                    print("Cannot deal with column %d having no dot" % (i+1))
                    print("LINE: ",l)
                    sys.exit(0)
                d1 = len(words[i]) - d - 1
                if d  > dot[i]:  dot[i]  = d
                if d1 > dot1[i]: dot1[i] = d1
                if dot[i] + dot1[i] + 1 > width[i]: width[i] = dot[i] + dot1[i] + 1
                # print width[i]

    # prepare a list with the proper format statement
    fmt = []
    for i in range(wlen):
        fmt.append('%%-%ds' % width[i])

    # finally loop over the table again, now printing out each
    # column with the correct width
    # @todo the dot ones need some spaces on either side...
    for l in lines:
        words = l.strip().split()
        s = ""
        for i in range(wlen):
            # print fmt[i], words[i]
            s = s + fmt[i] % words[i] + " "
        print(s.strip())

#

if __name__ == '__main__':
    if len(sys.argv) == 1:
        lines = readtable()
    else:
        lines = []
        for f in sys.argv[1:]:
            lines = lines + readtable(f)
    process_dot(lines)








