#! /usr/bin/env python
#

from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
import math
import sys

degrad = 57.2957795

def rotcurtab(file):
    """ reads a outputfile from tab= in rotcur
    if no comments in it, perhaps use popen("tabcomment $file -")
    """
    data = ascii.read(file)
    data.sort('col1')
    return data

def plot1(data, label='',efactor=1,rmax=2):
    r = data['col1']
    v   = data['col4'] 
    d_v = data['col5']*efactor
    inc   = data['col8']
    d_inc = data['col9']*efactor
    pa    = data['col6']
    d_pa  = data['col7']*efactor
    xpos  = data['col10']
    d_xpos= data['col11']*efactor
    ypos  = data['col12']
    d_ypos= data['col13']*efactor
    vsys  = data['col2']
    d_vsys= data['col3']*efactor
    vrms  = data['col15']
    vsini = v * np.sin(inc/degrad)
    fig = plt.figure()
    plt.suptitle('Rotation Curve Analysis :  (%s)' % label)
    #plt.subplots_adjust(bottom=0.14) 
    #plt.title('Tilted Ring Rotation Curve analysis')
    #ax0 = fig.add_subplot(1,1,1)
    #ax0.set_title('Tilted  Ring Rotation Curve Summary')
    ## 
    ax1 = fig.add_subplot(2,4,2)
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
	print "Center: ",xpos,ypos
	print "PA,INC: ",pa,inc
	print "Radius: ",r
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
	print r3
	fp.write(r3)
    fp.close()
    print data

def region_cgdisp(data,ds9,scale=0.0083333):
    """ create an overlay file for miriad::cgdisp 
    """


if __name__ == "__main__":
    name = sys.argv[1]
    rmax = float(sys.argv[2])
    data = rotcurtab(name)
    plot1(data,name,efactor=3.0,rmax=rmax)
    n = len(data)
    print "Found %d rings" % n
    m = 1
    region_ds9(data[0:n:m],'ds9.reg')
