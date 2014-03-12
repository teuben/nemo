#!/usr/bin/python

import matplotlib
import numpy as np
import matplotlib.colors
import matplotlib.pyplot as plt
import math
import time
import scipy.ndimage as ndi
from py_unsio import *            # import py_unsio package

# cnd line
import sys, getopt,os.path


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# main
def main(argv):
    prog = os.path.basename(argv.pop(0)) # program name
    file=''
    component='all'
    times='all'
    out=''
    xrange=20
    sigma=10.0
    
    try:
        opts,args=getopt.getopt(argv,"hi:o:c:t:r:s:",["in=","out","comp=","times=","range=","sigma="])

    except getopt.GetoptError:
        printHelp(prog)
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            printHelp(prog)
            sys.exit()
        elif opt in ("-i", "--in"):
            file = arg
        elif opt in ("-o", "--out"):
            out = arg        
        elif opt in ("-c", "--comp"):
            component = arg
        elif opt in ("-t", "--times"):
            times = arg
        elif opt in ("-r", "--range"):
            xrange = float(arg)
        elif opt in ("-s", "--sigma"):
            sigma = float(arg)


    if (file==''):
        print "\n\nyou must specify an input file name !!!\n"
        printHelp(prog)
    else:
        compute(file,out,component,times,xrange,sigma)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# printHelp
def printHelp(prog):
    help= """\

    ----------------------------------------------
    plot UNS simulation
    ----------------------------------------------

    Syntaxe : %s  -i <inputfile> -c <components> -t <times> -r <range> -s <sigma> -o <out_img>

    Notes :
        inputfile  : UNS snapshot with gas component
        components : component
        times      : selected time
        range      : -range<plot <range
        sigma      : gaussin sigma
	out_img    : if blank (default) display on screen, else on file

    """
    print help % (prog) 


def grid_density_gaussian_filter(x0, y0, x1, y1, w, h, x,y,r,rho):
    to = time.clock()	
    kx = (w - 1) / (x1 - x0)
    ky = (h - 1) / (y1 - y0)
    print "sigma r=",r
    #r=10
    border = 0#r
    imgw = (w + 2 * border)
    imgh = (h + 2 * border)
    img = np.zeros((imgh,imgw))
    img += 0.5
    #img+= (rho.min()/2.)
    #ix = 1
    #iy = 1
    #for x, y in data:
    #    ix = int((x - x0) * kx) + border
    #    iy = int((y - y0) * ky) + border
    #    if 0 <= ix < imgw and 0 <= iy < imgh:
    #        img[iy][ix] += 1

    myix= (x - x0)  * kx
    myix=myix.astype(int)
    myix=myix+border  # array of indexes values
    print "myix =", myix,myix.size
    maskix=(0<=myix)&(myix<imgw)
    print "maskix=",maskix
    print "myix[maskix]=",myix[maskix], myix[maskix].size,myix.size

    myiy= (y - y0) * ky
    myiy=myiy.astype(int)
    myiy=myiy+border
    print "myiy =", myiy,myiy.size
    maskiy=(0<=myiy)&(myiy<imgh)
    
    print "maskiy=",maskiy
    print "myiy[maskiy]=",myiy[maskiy], myiy[maskiy].size,myiy.size

    maskglob=maskix&maskiy

    print "maskglob=",maskglob,maskglob.size
    print "myix[maskglob]=",myix[maskglob],myix[maskglob].size
    print "myiy[maskglob]=",myiy[maskglob],myiy[maskglob].size
    #img[myix[maskglob]][myiy[maskglob]]=np.random.random_integers(1024)
    for a,b in zip(myiy[maskglob],myix[maskglob]): #myix[maskglob],myiy[maskglob]):
        img[a][b] += 1

    #img=math.log(img)
    #duplicated,unique_ind	= np.unique(myix[maskix],return_index=True)

    #print "duplicated,unique_ind",duplicated,duplicated.size,unique_ind, unique_ind.size
    print "timing grid_density_gaussian_filter = ", time.clock() - to
    
    return ndi.gaussian_filter(img, sigma=(r,r),order=0)  ## gaussian convolution
    #ndi.gaussian_filter(img, (r,r))  ## gaussian convolution

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# compute
# file      = input file (simulation name)
# component = selected components
# times     = selected times
def compute(file,out,component,times,xrange,sigma):
    t0=time.clock()
    uns=CunsIn(file,component,times)  # create UNS object
    cpt=0
    while(uns.nextFrame("")):
        cpt+=1
        comp=component
        # get data
        ok,pos=uns.getArrayF(comp,"pos") # positions
        print ok

        # get data
        ok,rho=uns.getArrayF(comp,"rho") # rho
        if ( ok!=True) :
            rho=np.ones(pos.size)

        # get data
        ok,mass=uns.getArrayF(comp  ,"mass") # gas density

        ok,timex=uns.getValueF("time")

        print ok
        print "Loading time : " , time.clock()-t0
        # reshape array in x,y,z arrays
        pos=np.reshape(pos,(-1,3))        # pos reshaped in a 2D array [nbody,3]
        x = pos[:,0]                      # x coordinates
        y = pos[:,1]                      # y coordinates
        z = pos[:,2]                      # z coordinates

        nbody = x.size
        ## if (rho.size > 0) :
        ##     rho = np.log(rho)
        ## else:
        ##     rho = z

        # center according COM
        xcom=np.average(x.astype(np.float64),weights=mass)
        ycom=np.average(y.astype(np.float64),weights=mass)
        zcom=np.average(z.astype(np.float64),weights=mass)
        print "COM = ",xcom,ycom,zcom

        x -= xcom
        y -= ycom
        z -= zcom

        if xrange>0:
            xmin = -xrange
            xmax =  xrange
            ymin = -xrange
            ymax =  xrange
        else:
            xmin = x.min()
            xmax = x.max()
            ymin = y.min()
            ymax = y.max()


        # data points xrange
        data_ymin = ymin
        data_ymax = ymax
        data_xmin = xmin
        data_xmax = xmax
        # view area xrange
        view_ymin = ymin
        view_ymax = ymax
        view_xmin = xmin
        view_xmax = xmax

        # get visible data points
        xlvis = x 
        ylvis = y 
        print "x.size= ", x.size, x[10], x[100]
        xx=x.size
        #for i in range(0,x.size):
        #    if view_xmin < x[i] < view_xmax and view_ymin < y[i] < view_ymax:
        #        xlvis.append(x[i])
        #        ylvis.append(y[i])

        fig = plt.figure(figsize=(20,20),dpi=150)#figsize=(1024,1024))
        
        DPI = fig.get_dpi()
        print "DPI:", DPI
        DefaultSize = fig.get_size_inches()
        print "Default size in Inches", DefaultSize
        print "Which should result in a %i x %i Image"%(DPI*DefaultSize[0], DPI*DefaultSize[1])

        # gaussian filter
        t0 = time.clock()
        zd = grid_density_gaussian_filter(view_xmin, view_ymin, view_xmax, view_ymax, 1024,1024, x, y, sigma,rho)
        plt.title(uns.getFileName())
        print "Filename = ",uns.getFileName()

        # plot component text
        ## xtext=(view_xmin+0.85*(view_xmax-view_xmin))
        ## ytext=(view_ymin+0.90*(view_ymax-view_ymin))
        ## plt.text(xtext,ytext,component,style='italic',color='black',bbox={'facecolor':'white', 'alpha':0.2, 'pad':10})

        # plot time text
        xtext=(view_xmin+0.05*(view_xmax-view_xmin))
        ytext=(view_ymin+0.90*(view_ymax-view_ymin))
        plt.text(xtext,ytext,"time:"+"%.3f"%timex+"\n"+component,style='italic',color='black',bbox={'facecolor':'white', 'alpha':0.2, 'pad':10})
        print ">> ",plt.get_cmap()
        plt.imshow(zd , origin='lower',cmap='jet',interpolation='nearest', norm = matplotlib.colors.LogNorm(), extent=[view_xmin, view_xmax, view_ymin, view_ymax])

        # to plot contour, comment out the following line
        #plt.contour(zd, origin='lower',levels=[0,1,2,3,4,5],cmap='winter',norm = matplotlib.colors.LogNorm(), extent=[view_xmin, view_xmax, view_ymin, view_ymax])

        if (out==''):
            plt.show()
        else:
            outfile=out+"."+"%05d"%cpt+".png"
            print ">> ",outfile
            #fig = plt.figure(figsize=(10,10),dpi=300)
            plt.savefig(outfile, bbox_inches=0)
        plt.close(fig)

if __name__=='__main__':
    main(sys.argv[0:])

    

