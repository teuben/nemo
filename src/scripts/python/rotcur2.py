#! /usr/bin/env python
#
# Generalize rotation curve plotter in order to compare rotation curves
#
# -u   rotcur format 
# -i   ringfit format
# -s   velfitss07 format  [not activated yet]
#
# -r   keep in radio format
# -o   convert from radio to optical convention (needs vsys)
#
# -z   (not implemented yet) convert everything to relativistic format
#
# key=value    allowed to override
#       
# Example of use:
#       rotcur2.py NGC2347 NGC2347.co.ringfit -o NGC2347.ha.ringfit -u -r try2.rotcur 
#
#
# Example rotcur2.txt file: (first non-comment line needs to be the column names!!!)
#<<
##  name             rmax   vsys  inc   w50  w90
##  NGC2347          30.0   4421  50.2  416  445
#>>
#

from __future__ import print_function

from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
import math
import os,sys

version = "25-apr-2017 PJT"

degrad = 57.2957795
c = 299792.458

parfile = 'rotcur2.txt'

def properties(name, file='rotcur2.txt'):
    """
    return a dictionary of properties for a galaxy/object by name
    """
    data = ascii.read(file)
    names = data['name'].data.tolist()
    idx = names.index(name)
    if True:
        print(data[idx])
        return data[idx]
    else:
        p = {}
        p['rmax'] = 30.0
        p['vsys'] = 4400.0
        p['inc'] = 50.0
        p['w50'] = 416.0
        p['w90'] = 445.0
        return p

def rotcurtab(file):
    """ reads a outputfile from tab= in rotcur
    if no comments in it, perhaps use popen("tabcomment $file -")
    """
    if True:
        os.system('tabcomment %s - > %s.tmp' % (file,file))
        data = ascii.read('%s.tmp' % file)
    else:
        data = ascii.read(file)
    data.sort('col1')
    return data

def get_ringplot(data, efactor=1):
    """
    radius from the center in arcsec, vrot, error, vrad, error, vsys, error. 
    """
    r     = data['col1']
    vt    = data['col2'] 
    d_vt  = data['col3']*efactor
    vr    = data['col4']
    d_vr  = data['col5']*efactor
    vsys  = data['col6']
    d_vsys= data['col7']*efactor
    return (r,vt)

def get_rotcur(data, efactor=1):
    """
    radius from the center in arcsec, vrot, error, vrad, error, vsys, error. 
    """
    r     = data['col1']
    vt    = data['col4'] 
    d_vt  = data['col5']*efactor
    return (r,vt)

def get_velfit(data, efactor=1):
    """
    runvelfitss07
    r   npt  vt   eVt Vr eVr Vm,t  eVm,t   Vm,r   eVm,r
    """
    r     = data['col1']
    vt    = data['col3'] 
    return (r,vt)

def junk():
    ax1.scatter(r,inc)
    ax1.errorbar(r,inc,yerr=d_inc,fmt='ro')
    ax1.set_title('Inclination (arcsec)')
    ax1.xaxis.label.set_size(10)
    ax1.yaxis.label.set_size(10)
    ax1.set_xlim([0,rmax])
    ax1.set_ylim([0,90])
    ##
    ax2 = fig.add_subplot(2,4,8)
    #ax2.set_title('RMS velocities in ring')
    ax2.scatter(r,vrms)
    ax2.set_title('RMS (km/s)')
    ax2.xaxis.label.set_size(10)
    ax2.yaxis.label.set_size(10)
    ax2.set_xlim([0,rmax])
    ax2.set_ylim([0,10])
    ##
    ax3 = fig.add_subplot(2,4,1)
    ax3.scatter(r,v)
    ax3.errorbar(r,v,yerr=d_v,fmt='ro')
    ax3.set_title('Vrot (km/s)')
    ax3.xaxis.label.set_size(10)
    ax3.yaxis.label.set_size(10)
    ax3.set_xlim([0,rmax])
    ax3.set_ylim([0,110])
    ##
    ax4 = fig.add_subplot(2,4,5)
    #ax2.set_title('V.sin(INC)')
    ax4.scatter(r,vsini)
    ax4.set_xlabel('Radius (arcsec)')
    ax4.set_title('V.sin(i) (km/s)')
    ax4.xaxis.label.set_size(10)
    ax4.yaxis.label.set_size(10)
    ax4.set_xlim([0,rmax])
    ax4.set_ylim([0,110])
    ##
    ax5 = fig.add_subplot(2,4,3)
    ax5.set_title('X-center')
    ax5.scatter(r,xpos)
    ax5.errorbar(r,xpos,yerr=d_xpos,fmt='ro')
    ax5.xaxis.label.set_size(10)
    ax5.yaxis.label.set_size(10)
    ax5.set_xlim([0,rmax])
    ax5.set_ylim([150,400])
    #
    ax6 = fig.add_subplot(2,4,7)
    ax6.set_title('Y-center')
    ax6.scatter(r,ypos)
    ax6.errorbar(r,ypos,yerr=d_ypos,fmt='ro')
    ax6.xaxis.label.set_size(10)
    ax6.yaxis.label.set_size(10)
    ax6.set_xlim([0,rmax])
    ax6.set_ylim([250,500])
    #
    ax7 = fig.add_subplot(2,4,4)
    ax7.set_title('VSYS')
    ax7.scatter(r,vsys)
    ax7.errorbar(r,vsys,yerr=d_vsys,fmt='ro')
    ax7.set_xlim([0,rmax])
    ax7.set_ylim([140,190])
    ax7.xaxis.label.set_size(10)
    ax7.yaxis.label.set_size(10)
    #
    ax8 = fig.add_subplot(2,4,6)
    ax8.set_title('PA')
    ax8.scatter(r,pa)
    ax8.errorbar(r,pa,yerr=d_pa,fmt='ro')
    ax8.set_xlim([0,rmax])
    ax8.set_ylim([0,90])
    ax8.xaxis.label.set_size(10)
    ax8.yaxis.label.set_size(10)
    #
    plt.show()
    fig.savefig('junk.pdf')

def region_ds9(data,ds9,scale=0.0083333):
    """ create a ds9 region of the ellipses found in rotcur solution 
    Also needs the scale to convert to pixels in a map
    """
    (xpos,ypos)  = data['col10'],data['col12']
    inc = data['col8']
    pa  = data['col6']
    r   = data['col1']/3600.0
    maj = r / scale
    min = r*np.cos(inc/degrad) / scale
    if False:
        print("Center: ",xpos,ypos)
        print("PA,INC: ",pa,inc)
        print("Radius: ",r)
    r1='global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
    r2='image\n'
    #
    fp = open(ds9,'w')
    fp.write('# Region file format: DS9 version 4.1\n')
    fp.write('# Filename: smc.mom1f/fits\n')
    fp.write(r1)
    fp.write(r2)
    for i in range(len(data)):
        r3='ellipse(%f,%f,%f,%f,%f) # color=black width=2\n' % (xpos[i],ypos[i],min[i],maj[i],pa[i])
        print(r3)
        fp.write(r3)
    fp.close()
    print(data)

def region_cgdisp(data,ds9,scale=0.0083333):
    """ create an overlay file for miriad::cgdisp 
    """

def plabel(umode,scale):
    """
    """
    if umode:
        lab = 'u'
    else:
        lab = 'i'
    if scale:
        lab = lab + '-o'
    else:
        lab = lab + '-r'
    return lab

def print_usage(argv):
    print("Multi-table rotation curve plotter and comparisons - version %s " % version)
    print("Usage:")
    print("%s name  [key=val] [-u] [-i] [-r] [-o] curve1  [-u] [-i] [-r] [-o] curve2 ..." % argv[0])
    print("   name     Required name, an ID grab default values from %s" % parfile)
    print("   key=val  Optional re-assigned keyword")
    print("   -u       rotcur type table  (for this and all following tables until reset) ")
    print("   -i       ringfit type table  (for this and all following tables until reset) ")
    print("   -r       radio convention  (for this and all following tables until reset) ")
    print("   -o       optical convention  (for this and all following tables until reset) ")
    print("Currently all curves are *plotted* in the radio convention")
    print("")
    print("In addition, for a limited number of keywords, a new value can be given:")
    print("   rmax")
    print("   vsys")
    print("   inc")
    print("   w50")
    print("   w90")

    sys.exit(0)

    

if __name__ == "__main__":
    print("LEN:",len(sys.argv))
    if len(sys.argv) < 2: print_usage(sys.argv)
    gal = sys.argv[1]
    p = properties(gal)
    rmax = p['rmax']
    vsys = p['vsys']
    inc = p['inc']
    w50 = p['w50']
    w90 = p['w90']
    #
    fig = plt.figure()
    plt.title('%s   :   VSYS=%g    INC=%g' % (gal,vsys,inc))
    ax = fig.add_subplot(1,1,1)
    scale = False      # scale from optical to radio convention?   (-o and -r)
    umode = False      # -u: rotcur format       -i: ringfit format (default)
    for name in sys.argv[2:]:
        if name.find('=') > 0:
            print("EXEC: ",name)
            exec(name)
            continue
        if name=='-o':
            scale = True
            continue
        if name=='-r':
            scale = False
            continue
        if name=='-u':
            umode = True
            continue
        if name=='-i':
            umode = False
            continue
        data = rotcurtab(name)
        if umode:
            (r1,v1) = get_rotcur(data)      # 'u'
        else:
            (r1,v1) = get_ringplot(data)    # 'i'
            #(r1,v1) = get_velfit(data)     # 's'
        if scale:
            o2r = 1.0-2.0*vsys/c
            v1 = v1 * o2r
        n = len(data)
        print("Found %d radii for %s" % (n,name))
        ax.plot(r1,v1,label="%s[%s]" % (name,plabel(umode,scale)))
    (rmin,rmax) = ax.get_xlim()
    (vmin,vmax) = ax.get_ylim()
    ax.set_xlim([0.0,rmax])
    ax.set_ylim([0.0,vmax])
    sini = math.sin(inc*math.pi/180.0)
    v50 = 0.5*w50/sini
    v90 = 0.5*w90/sini
    print(v50,v90)
    ax.plot([0.9*rmax,rmax],[v50,v50],'k-',label='HI-W50',linestyle='-',linewidth=2)
    ax.plot([0.8*rmax,rmax],[v90,v90],'k-',label='HI-W90',linestyle='-',linewidth=2)
    ax.plot([0.0,0.1*rmax],[v50,v50],'k-',linestyle='-',linewidth=2)
    ax.plot([0.0,0.2*rmax],[v90,v90],'k-',linestyle='-',linewidth=2)
    ax.legend(loc='best',prop={'size':8})
    ax.set_xlabel('Radius (arcsec)')
    ax.set_ylabel('Velocity (km/s)')
    plt.show()
