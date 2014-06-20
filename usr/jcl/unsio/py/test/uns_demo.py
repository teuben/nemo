# -----------------------------------------------------------------------
#!/usr/bin/python
# -----------------------------------------------------------------------
# The following program shows how to use UNSIO library from
# a Python program
#
# This program reads an unsio compatible snapshot from the command line
# and save it in gadget2 format
#
# Syntaxe : uns_demo.py -i <inputfile> -o <outfile> -c <components> -t <times>
#
# inputfile  : UNS input snapshot
# outfile    : gadget2 out filename
# components : specify just one or a coma separeted list of components
#              among => disk,gas,stars,halo,bulge,bndry
#              example : -c disk,stars
# times      : selected time
#              for all time steps, use "-t all"
#              for one specific time, use "-t 10"
#              for a range of times, use "-t 5:6"
#
# -----------------------------------------------------------------------
# For more information about how to use UNSIO, visit:
# http://projets.lam.fr/projects/unsio/
# -----------------------------------------------------------------------
#  Copyright Jean-Charles LAMBERT (CeSAM)- 2008-2014
#
#  e-mail:   Jean-Charles.Lambert@lam.fr                                      
#  address:  Centre de donneeS Astrophysique de Marseille (CeSAM)         
#            Laboratoire d'Astrophysique de Marseille                          
#            Pole de l'Etoile, site de Chateau-Gombert                         
#            38, rue Frederic Joliot-Curie                                     
#            13388 Marseille cedex 13 France                                   
#            CNRS U.M.R 6110
# -----------------------------------------------------------------------

# unsio module loading
# ==> do not forget to update PYTHONPATH environment variable with
#     py_unsio location path
from py_unsio import *
import numpy as np
# cmd line
import sys, getopt

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# printAndSaveProp
# print properties for the component comp given as argument
# and save them
def printAndSaveProp(uns,unso,comp):
    info="""\
    ----------------------------------------------
    Component : [%s]
    ----------------------------------------------
    """
    print info % (comp)
    # return a 1D numpy data array with mass
    ok,mass=uns.getArrayF(comp,"mass")
    if ok:
        print "mass =",mass
        unso.setArrayF(comp,"mass",mass) # save mass
        
    # return a 1D numpy data array with pos
    ok,pos=uns.getArrayF(comp,"pos")
    if ok:
        print "pos =",pos
        unso.setArrayF(comp,"pos",pos) # save pos
        
    # return a 1D numpy data array with vel
    ok,vel=uns.getArrayF(comp,"vel")
    if ok:
        print "vel =",vel
        unso.setArrayF(comp,"vel",vel) # save vel

    # return a 1D numpy data array with rho
    ok,rho=uns.getArrayF(comp,"rho")
    if ok:
        print "rho =",rho
        unso.setArrayF(comp,"rho",rho) # save rho

    # return a 1D numpy data array with temperature
    ok,temp=uns.getArrayF(comp,"temp")
    if ok:
        print "temp =",temp
        unso.setArrayF(comp,"temp",temp) # save temperature

    # return a 1D numpy data array with hsml
    ok,hsml=uns.getArrayF(comp,"hsml")
    if ok:
        print "hsml =",hsml
        unso.setArrayF(comp,"hsml",hsml) # save hsml

    # return a 1D numpy data array with particles age
    ok,age=uns.getArrayF(comp,"age")
    if ok:
        print "age =",age
        unso.setArrayF(comp,"age",age) # save age

    # return a 1D numpy data array with mettalicity
    ok,metal=uns.getArrayF(comp,"metal")
    if ok:
        print "metal =",metal
        unso.setArrayF(comp,"metal",metal) # save mettalicity


    # return a 1D numy data array with id
    ok,indexes=uns.getArrayI(comp,"id")
    if ok:
        print "indexes =", indexes
        unso.setArrayI(comp,"id",indexes) # save id

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# compute
# loop on all selected time steps and print out properties
# for every components from "comp" variable
# and save them in a gadegt2 file
def compute(file,out,comp,times):
    print "file=",file, " outfile=",out," comp=",comp, " times=",times

    # instantiate a CunsIn object, here we request to load "all" components
    uns=CunsIn(file,"all",times)

    # load frame
    cpt=0
    while (uns.nextFrame("")):  # load every snasphots

        nameout=out+".%05d"%cpt # create out filename
        print "Filename Out =",nameout
        # instantiate a CunsOut object in "gadget2" format 
        unso=CunsOut(nameout,"gadget2")

        ok,tsnap=uns.getValueF("time") # return snasphot time
        print "Snapshot time : ","%.03f"%tsnap
        unso.setValueF("time",tsnap) # save snapshot time
        
        # loop on all components stored in comp variable
        for onecomp in (comp.split(",")):
            printAndSaveProp(uns,unso,onecomp) # print properties for the component

        unso.save()  # save snasphot
        cpt+=1 # one more frame
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# main
def main(argv):
    prog = argv.pop(0) # program name
    infile=''
    out=''
    comp='gas,disk,stars,halo,bulge,bndry'
    times='all'

    try:
        opts,args=getopt.getopt(argv,"hi:o:c:t:",["in=","out=","comp=","times="])

    except getopt.GetoptError:
        print "\nUnknown parameters, please check ....\n\n"
        printHelp(prog)
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            printHelp(prog)
            sys.exit()
        elif opt in ("-i", "--in"):
            infile = arg
        elif opt in ("-o", "--out"):
            out = arg
        elif opt in ("-c", "--comp"):
            comp = arg
        elif opt in ("-t", "--times"):
            times = arg


    if (infile != '' and out != ''):
        compute(infile,out,comp,times)
    else:
        print "\nYou must provide input and output files.....\n\n"
        printHelp(prog)
        
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# printHelp
def printHelp(prog):
    help= """\
    --------------------------------------------------
    Print out components properties of an UNS file
    and save them in a gadget2 file
    --------------------------------------------------
    
    Syntaxe : %s  -i <inputfile> -o <outfile> -c <components> -t <times>
    Example : %s  -i gtr118_1912 -o myfile.g2 -c gas,stars
    
    Notes :
        inputfile  : UNS input snapshot

        outfile    : gadget2 out filename
        
        components : specify just one or a coma separeted list of components
                     among => disk,gas,stars,halo,bulge,bndry
                     exemple : -c disk,stars
                     
        times      : selected time
                     for all time steps, use "-t all"
                     for one specific time, use "-t 10"
                     for a range of times, use "-t 5:6"

    """
    print help % (prog,prog)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# main
if __name__ == "__main__":
    main(sys.argv[0:])
