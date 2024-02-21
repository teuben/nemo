#! /usr/bin/env python
#
#  Average spectrum in a GBT data set at given (ra,dec) size and (dra,ddec)
#
#  Either a single FITS cube can be given, or one or more calibrated SDFITS files
#
#  Peter Teuben - 17-mar-2022 - Created
#                 17-oct-2022 - various updates
#                 26-oct-2022 - write out spectrum for plotsp3.py
#                 28-oct-2022 - using new nemopy.getparam
#
# Make a comma separated list of files in unix:
#     ls NGC0001/NGC0001_scan_* | sed -z 's/\n/,/g;s/,$/\n/'

import os
import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from astropy.io import fits
from spectral_cube import SpectralCube
import astropy.units as u
from nemopy import getparam
#import pdb
#pdb.set_trace()


keyval = [
    "in=\n          input fits cube",
    "ra=\n          RA   hh:mm:ss.s (or decimal degrees) [ref pixel]",
    "dec=\n         DEC  dd:mm:ss.s (or decimal degrees) [ref pixel]",
    "size=20\n      Size of area in arcsec (FWHM)",
    "dra=0\n        RA offset in arcsec",
    "ddec=0\n       DEC offset in arcsec",
    "reg=\n         Alternate 'ds9' region file",
    "sdfits=\n      optional sdfits file(s) (not re-implemented yet)",
    "tab=\n         If given, write out spectrum here",
    "vrange=\n      Plotting range in velocity",
    "irange=\n      Plotting range in intensity",
    "iscale=1\n     Intensity scaling factor",
    "blorder=0\n    Order of baseline fits (if >= 0)",
    "blregion=\n    Pairs of sections along velocity axis for baseline fit",
    "vshift=0\n     Apply velocity shift to the SDFITS data",
    "vlsr=\n        Place VLSR marker if this is set",
    "winc=0\n       Window smoothing applied to FITS data cube",
    "wins=11\n      Window smoothing applied to SDFITS spectral data",
    "edge=5\n       Number of edge channels to exclude",
    "plot=raw\n     Plot 'r(aw)' or 's(ubtracted)' spectrum - 's' will show flux",
    "savefig=\n     Save figure instead of interactive pllot",
    "VERSION=0.7\n  14-dec-2023 PJT",
]

usage = """
plot a spectrum from a cube and/or set of sdfits files

this script can take a FITS cube and plot a spectrum in a region defined
by an offset from an (RA,DEC) position and a  size.  Optionally it can
do a baseline fit in a blregion= and plot the residuals

It will also do SDFITS files once this is fully enabled again
"""

p = getparam.Param(keyval,usage)

def stats3(x):
    n = len(x)
    n1 = n//3
    n2 = 2*n1
    x1 = x[:n1]
    x2 = x[n1:n2]
    x3 = x[n2:]
    print("stats  : Showing mean/std/min/max in 3 sections of spectrum")
    print("stats-1: ",x1.mean(), x1.std(), x1.min(), x1.max())
    print("stats-2: ",x2.mean(), x2.std(), x2.min(), x2.max())
    print("stats-3: ",x3.mean(), x3.std(), x3.min(), x3.max())
    print("stats  : Warning; units are mK really")

def sexa2deci(s, scale=1.0):
    """ s is a hms/dms string   12.34 or 12:34
        returns the float value
        if it's already a float, returns "as is", no
        factor 15 conversion  (this is the FITS convention)
    """
    if s.find(':') > 0:
        dms = [float(x) for x in s.split(':')]
        if len(dms) == 2:
            dms.append(0)
        if s[0] == '-':
            sign = -1
        else:
            sign = +1
        r = abs(dms[0]) + (dms[1] + dms[2]/60.0)/60.0
        return r*scale*sign
    else:
        return float(dms)

def read_ds9_region(reg):
    """ read a region file
    """
    lines = open(reg,"r").readlines()
    ds9_region = ""
    for line in lines:
        if line[0] == '#': continue
        if line.find('global') == 0: continue
        if len(ds9_region) == 0:
            ds9_region = line.strip()
        else:
            ds9_region = ds9_region + '; ' +  line.strip()
    return ds9_region

def nr_smooth(x,window_len=11,window='hanning'):

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='same')
    return y

def my_smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def fit_poly(x, y, p_order=1, bl = []):
    """ from array X between Xmin and Xmax fit a polynomial
    """

    if len(bl) == 0:
        p = np.poly1d(np.polyfit(x,y,p_order))
        t = x
        r = y - p(x)
    else:
        first = True
        for b in bl:
            # print('B',b)
            if first:
                m = ((x>b[0]) & (x<b[1]))
                first = False
            else:
                m = m | ((x>b[0]) & (x<b[1]))
                
        p = np.poly1d(np.polyfit(x[m],y[m],p_order))
        t = x[m]
        r = y[m] - p(x[m])
    return (p,t,r)

def diff_rms(y):
    """ take the differences between neighboring signals
    and compute their rms. this should be sqrt(2)*sigma
    if there is no  trend in the input signal, and if
    the input signal is not correlated (e.g. hanning)
    """
    #y1 = y[1:]
    #y2 = y[:-1]
    return (y[1:]-y[:-1]).std() / 1.414

def add_spectrum(filename):
    (v,t) = np.loadtxt(filename).T
    plt.plot(v,t,label=filename)


def flux_spectrum(x,y,xstart,xend):
    """  measure the flux in the (x,y) spectrum between xrange[0] and xrange[1]
    """
    dx = x[1]-x[0]
    if dx<0:
        print("flux_spectrum channel width < 0 should not happen",dx)
    m = ma.masked_outside(x,xstart,xend)
    sm = ma.MaskedArray(y,m.mask)
    flux = dx * sm.sum()
    return flux
    # print("Flux between %g %g = %g\n" % (xstart,xend,flux))

# --------------------------------------------------------------------------

if p.has("ra") and p.has("dec"):
    needr  = False
    ra0    = sexa2deci(p.get("ra"), 15.0)
    dec0   = sexa2deci(p.get("dec"))
else:
    needr  = True
    print("Using reference pixel as ra=,dec=")

size   = float(p.get("size"))
dra    = float(p.get("dra"))
ddec   = float(p.get("ddec"))


blorder = int(p.get("blorder"))
if p.has("blregion"):
    blregion = p.listf("blregion")
    nbl = len(blregion)
    if nbl%2 != 0:
        print("Need even number of baseline sections")
        sys.exit(0)
    bl = []
    nbl = nbl // 2
    for i in range(nbl):
        vmin = blregion[2*i]
        vmax = blregion[2*i+1]
        bl.append((vmin,vmax))
    print("bl=",bl)
else:
    nbl = 0
print("nbl=",nbl)

iscale = float(p.get('iscale'))
winc = int(p.get("winc"))
wins = int(p.get("wins"))
    
radmax = size / 2 /  3600.0
edge = int(p.get("edge"))
c = 299792.458
restfr = -1

Qplot = True
if p.get("plot")[0] == 'r':
    Qraw = True
elif p.get("plot")[0] == 's':
    Qraw = False
else:
    Qraw = True
    print("Warning: invalid plot=, assuming raw mode")



#  test if a cube
ff = p.get("in")
hdu = fits.open(ff)
if hdu[0].header['NAXIS'] > 2:      
    cube = SpectralCube.read(ff).with_spectral_unit(u.km/u.s)
    if needr:
        needr = False
        
        ra0 = hdu[0].header['CRVAL1']
        dec0 = hdu[0].header['CRVAL2']
        cosd0 = np.cos(dec0*np.pi/180)
        dec0  = dec0 + ddec/3600.0
        ra0   = ra0 + dra/3600.0/cosd0
        print("PJT",ra0,dec0,restfr)

        bmaj = hdu[0].header['BMAJ']*3600
        bmin = hdu[0].header['BMIN']*3600
        print("BEAM: %g x %g arcsec" % (bmaj,bmin))
            

    if len(p.get("reg")) > 0:
        ds9 = read_ds9_region(p.get("reg"))
    else:
        ds9 = 'fk5; circle(%.10f, %.10f, %g")'  % (ra0,dec0,size/2)
    print('DS9',ds9)
    spectrum = cube.subcube_from_ds9region(ds9).mean(axis=(1, 2))
    spec1 = iscale * 1000 * spectrum 

    crpix3 = hdu[0].header['CRPIX3']
    crval3 = hdu[0].header['CRVAL3']
    cdelt3 = hdu[0].header['CDELT3']
    if 'RESTFRQ' in hdu[0].header:
        restfr = hdu[0].header['RESTFRQ']
    elif 'RESTFREQ' in hdu[0].header:
        restfr = hdu[0].header['RESTFREQ']        
    nchan = len(spectrum)
    chan  = np.arange(0,nchan)+1
    vrad1 = (chan-crpix3) * cdelt3 + crval3
    gal   = ff.split('_')[0]    
    gal1  = ff
    dv = vrad1[1]-vrad1[0]
    if dv < 0:
        vrad1 = np.flip(vrad1)
        spec1 = np.flip(spec1)
    print("Found FITS",nchan," channels ",vrad1[1]-vrad1[0], " vel ", vrad1[0],vrad1[nchan-1])    
else:
    spec1 = None
sdfits = p.list('sdfits')


# Spectra, units are mK
# spec1 :   spectrum from a FITS cube (w/ nchan and chan)
# spec2 :   from SDFITS file(s)

spec2 = None
nspec = 0
nchan = 0
vshift = float(p.get("vshift"))

# if no cube, assume they are SDFITS file(s)

for ff in sdfits:
    hdu = fits.open(ff)
    d2 = hdu[1].data
    crpix1 = d2['CRPIX1'][0]
    crval1 = d2['CRVAL1'][0]
    cdelt1 = d2['CDELT1'][0]
    restfr = d2['RESTFREQ'][0]
    restfr = 115.2712018 * 1e9
    tsys   = d2['TSYS'][0]
    if nchan == 0:
        # first time around, find the number of channels
        # we are just adding all spectra, no WCS check; make sure crval1 doesn't vary more than cdelt1
        nchan = len(d2['DATA'][0])
        chan  = np.arange(0,nchan) + 1
        spec2 = np.zeros(nchan)
        freq  = (chan - crpix1)  * cdelt1 + crval1
        vrad2 = (1-freq/restfr)*c - vshift
        nspec = 0
        gal   = ff.split('_')[0]
        gal2  = ff
        print("Found SDFITS",nchan," channels ",vrad2[1]-vrad2[0], " vel ", vrad2[0],vrad2[nchan-1])
    ra = d2['CRVAL2']
    dec = d2['CRVAL3']
    idx2 = ra[ra==0]
    idx3 = dec[dec==0]
    rad = np.sqrt( (ra-ra0)**2 + (dec-dec0)**2)
    r2c = rad[rad<radmax]
    idx = np.where(rad<radmax)[0]
    print(ff,len(idx),tsys,crval1,crpix1,cdelt1)
    nspec = nspec + len(idx)
    spec2 = spec2 + d2['DATA'][idx].sum(axis=0)

if nspec > 0:
    print("Found %d spectra from SDFITS file(s)" % nspec)
    spec2 = 1000 * spec2 / nspec
    dv = vrad2[1]-vrad2[0]
    if dv < 0:
        vrad2 = np.flip(vrad2)
        spec2 = np.flip(spec2)
    if False:
        tab = "plot_spectrum.txt"
        fp = open(tab,"w")
        # fp.write("# generated by ",sys.argv)
        for (v,s) in zip(vrad2[edge:-edge], spec2[edge:-edge]):
            fp.write("%g %g\n" % (v,s))
        fp.close()
        print("Wrote %s average of %d spectra from SDFITS file(s)" % (tab,nspec))

f1 = -1
f2 = -1
smooth = nr_smooth


if Qplot:
    plt.figure()
    xlim = []

    if False:
        # disable the cube 
        spec1 = None

    # spectrum from FITS cube
    if spec1 != None:
        if winc > 0:
            spec1s = smooth(spec1,winc)
            vrad1s = smooth(vrad1,winc)
            print('SMOOTH: %d %d' % (len(spec1),len(spec1s)))
            nchan1 = len(vrad1)
            vcen1 = vrad1[nchan1//2]
            icen1  = np.where(vrad1s==vcen1)[0][0]
            print("PJT-10",nchan1//2,icen1)
            spec1  = spec1s[icen1-nchan1//2+edge:icen1+nchan1//2-edge]
            vrad1  = vrad1s[icen1-nchan1//2+edge:icen1+nchan1//2-edge]
            
        if Qraw:
            plt.plot(vrad1,spec1,label=gal1)
        if p.has("vrange"):
            xlim = p.listf("vrange")
        else:
            xlim = [vrad1[0],vrad1[-1]]
        plt.plot(xlim,[0,0],c='black')
        stats3(spec1)
        # subtract baseline?
        if nbl>0 and blorder>=0:
            (p2,t2,r2) = fit_poly(vrad1, spec1.value, blorder, bl)
            rms2 = r2.std()
            rms3 = diff_rms(r2)
            if Qraw:
                plt.plot(vrad1,p2(vrad1),'-',label='POLY %d SMTH %d' % (blorder,-1))
                plt.plot(t2, r2, '-', label='RMS %.3g %.3g' % (rms2, rms3))
            else:
                resid = spec1.value-p2(vrad1)
                f1 = flux_spectrum(vrad1,resid,bl[0][1],bl[1][0])/1000
                plt.plot(vrad1,resid,'-',label='RESID fits; flux=%g' % f1)
            # plt.plot([v2[0],v2[-1]], [0.0, 0.0], c='black', linewidth=2, label='baseline BAND %d' % do_band)
        
    # spectra from SDFITS
    if nspec > 0:
        if wins > 0:
            spec2s = smooth(spec2,wins)
            vrad2s = smooth(vrad2,wins)
            print('SMOOTH: %d %d' % (len(spec2),len(spec2s)))
            if False:
                y = spec2s[edge:-edge]
                x = vrad2[:len(y)]
                if Qraw:
                    plt.plot(x,y,label=gal2+'-s')
            else:
                nchan2 = len(vrad2)
                vcen = vrad2[nchan2//2]
                icen2  = np.where(vrad2s==vcen)[0][0]
                print("PJT-10",nchan2//2,icen2)
                spec2  = spec2s[icen2-nchan2//2+edge:icen2+nchan2//2-edge]
                vrad2  = vrad2s[icen2-nchan2//2+edge:icen2+nchan2//2-edge]
                plt.plot(vrad2,spec2,label=gal2+'-s')                
        else:
            if Qraw:
                plt.plot(vrad2[edge:-edge], spec2[edge:-edge],label=gal2)

        # draw black baseline at T=0
        plt.plot([vrad2[0],vrad2[-1]], [0,0], c='black')
        if len(xlim) == 0:
            xlim = [vrad2[0],vrad2[-1]]

        stats3(spec2)

        # subtract baseline?
        if nbl>0 and blorder>=0 and wins > 0:
            (p3,t3,r3) = fit_poly(vrad2, spec2, blorder, bl)
            rms2 = r3.std()
            rms3 = diff_rms(r3)
            if Qraw:
                plt.plot(vrad2,p3(vrad2),'-',label='POLY %d SMTH %d SDFITS' % (blorder,-1))
                plt.plot(t3, r3, '-', label='RMS %.3g %.3g SDFITS' % (rms2, rms3))
                plt.plot([bl[0][1],bl[1][0]],[-4*rms2,-4*rms2], c='black')                
            else:
                resid2 = spec2-p3(vrad2)
                f2 = flux_spectrum(vrad2,resid2,bl[0][1],bl[1][0])/1000
                plt.plot(vrad2,resid2,'-',label='RESID sdfits; flux=%g' % f2)
                plt.plot(vrad2,p3(vrad2),'-',label='POLY %d SMTH %d SDFITS' % (blorder,-1))
            # plt.plot([v2[0],v2[-1]], [0.0, 0.0], c='black', linewidth=2, label='baseline BAND %d' % do_band)

    #plt.plot([bl[0][1],bl[1][0]],[-4.0*rms2,-4.0*rms2], c='black')
    #plt.plot([bl[0][1],bl[1][0]],[-4.5*rms2,-4.5*rms2], c='black')
    
    if f1 > 0:
        print("Flux   FITS = %g K.km/s in %g arcsec beam" % (f1,size))
    if f2 > 0:
        print("Flux SDFITS = %g K.km/s in %g arcsec beam" % (f2,size))
        
    plt.xlabel('Vrad (km/s)')
    plt.ylabel('T (mK)')
    if xlim[0] < xlim[1]:
        plt.xlim(xlim)
    else:
        plt.xlim([xlim[-1],xlim[0]])
    if p.has("irange"):
        ylim = p.listf("irange")
        plt.ylim(ylim)
    if p.has("vlsr"):
        vlsr = float(p.get("vlsr"))
        ax = plt.gca()
        fmax = ax.get_ylim()[1]
        plt.arrow(vlsr,fmax,0.0,-0.5*fmax,
                  head_width=20, head_length=0.1*fmax,
                  length_includes_head=True, facecolor='red')
        plt.annotate('VLSR=%g' % vlsr, xy=(vlsr, fmax), multialignment='center')
    if p.has("reg"):
        plt.title('%s @ %f %f region %s"' % (gal,ra0,dec0,p.get("reg")))
    else:
        plt.title('%s @ %f %f size %g"' % (gal,ra0,dec0,size))
    plt.legend(loc='best')
    if p.has("savefig"):
        plt.savefig(p.get("savefig"))
    else:
        plt.show()
    #

if p.has('tab'):
    sp_out = p.get('tab')
    sp_data = np.squeeze(np.dstack((vrad1,spec1.value)))
    np.savetxt(sp_out,sp_data,fmt='%.4f')
    print("Written ",sp_out)

# @todo:   weight by Tsys
