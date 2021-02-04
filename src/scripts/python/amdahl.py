#! /usr/bin/env python
#
#  Amdahl's law:   a fraction of the code can be parallized
#  Output are, for each N
#       r_u,  r_e, f_p
#  where both r_u (user_ratio) and r_e (elapsed_ratio) should ideally
#  be 1, else there is overhead of running in parallel, and
#  f_p (parallel fraction) is the derived fraction of code that's parallized.
#  If e.g. f_p depends a lot on N, something is not ideal.
#
#  Code that does parallize has r_u=1 r_e=N and f_p=0
#  Perfect parallel code        r_u=1 r_e=1 and f_p=1
#  
#
#  nov 22, 2020:   first version
#  jan 31, 2021:   cleanup up, generalized for files that contain
#                  N <time_output>
#                  1 15.16user 0.00system 0:15.18elapsed   99%CPU
import sys, os

def par(u1,e1,np,un,en):
    """   amdahl's law
    returns:
    r1 = user ratio (r_u)
    r2 = elapsed ratio (r_e)
    fp = parallel fraction (f_p)
    """
    r1 = un/u1
    r2 = en/(e1/np)
    if np > 1:
        fp = (np-r2)/(np-1)
    else:
        fp = 1.0
    return (r1,r2,fp)

def hms2s(hms):
    w = hms.split(':')
    if len(w) == 1:
        return float(w[0])
    elif len(w) == 2:
        return float(w[0])*60 + float(w[1])
    elif len(w) == 3:
        return float(w[0])*3600 + float(w[1])*60 + float(w[0])
    else:
        return -1.0

def readf1(filename):
    """  example line:
         2  72.24 10.58 1:47.36
         1  19.72user 0.13system 0:19.95elapsed   99%CPU
    """
    lines = open(filename).readlines()
    for line in lines:
        if line[0] == '#':
            print(line.strip())
            continue
        pline = line[:len(line)-1]
        line = line.replace('user','').replace('system','').replace('elapsed','').replace('%CPU','')
        w = line.strip().split()
        n = int(w[0])
        u0 = float(w[1])
        s0 = float(w[2])
        e0 = hms2s(w[3])
        if n == 1:
            un0 = u0
            sn0 = s0
            en0 = e0
        (r10,r20,fp0) = par(un0,en0,n,u0,e0)
        # print("%d %g %g %g  %g %g %g" % (n,u0,s0,e0, r10,r20,fp0))
        print("%s %.2f %.2f %.2f" % (pline,r10,r20,fp0))
    

def readf2(filename):
    """  example line:
         2  72.24 10.58 1:47.36   882.77  33.23 14:37.16
    """
    lines = open(filename).readlines()
    for line in lines:
        if line[0] == '#': continue
        w = line.strip().split()
        n = int(w[0])
        u0 = float(w[1])
        s0 = float(w[2])
        e0 = hms2s(w[3])
        u1 = float(w[4])
        s1 = float(w[5])
        e1 = hms2s(w[6])
        if n == 1:
            un0 = u0
            sn0 = s0
            en0 = e0
            un1 = u1
            sn1 = s1
            en1 = e1
        (r10,r20,fp0) = par(un0,en0,n,u0,e0)
        (r11,r21,fp1) = par(un1,en1,n,u1,e1)
        print(n,u0,s0,e0,u1,s1,e1, r10,r20,fp0, r11,r21,fp1)

        

readf1(sys.argv[1])        
