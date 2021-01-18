/*
 *  rudimentary support to write a template script
 *  for tabplot and tabhist
 *
 *  Typical usage:
 *       str = pyplot_init("sample.py");
 *       pyplot_plot(str,"mytable.tab",xcol,ycol,dycol=None,xrange,yrange)
 *       pyplot_hist(str,"mytable.tab",xcol,xrange,bins)
 *       pyplot_close(str)
 */

#include <nemo.h>

stream pyplot_init(string fname)
{
  stream str = stropen(fname,"w");
  dprintf(0,"Writing template python plotting script template %s\n",fname);
  fprintf(str,"#! /usr/bin/env python\n");
  fprintf(str,"#  template script written by NEMO::pytable\n");
  fprintf(str,"#\n");
  fprintf(str,"import os\n");
  fprintf(str,"import sys\n");
  fprintf(str,"import numpy as np\n");
  fprintf(str,"import matplotlib.pyplot as plt\n");
  fprintf(str,"\n");
  return str;
}

void pyplot_close(stream str)
{
  strclose(str);
}

//  A4 is 8.27" x 11.69" = 210 x 297 mm 
//  USA   8.5   x 11     = 216 x 279 mm
void pyplot_plot(stream str, string tabname, int *xcol, int *ycol, int *dycol, real *xrange, real *yrange)
{
  fprintf(str,"tab1 = '%s'\n",tabname);
  fprintf(str,"data1 = np.loadtxt(tab1).T\n");
  fprintf(str,"print(data1.shape)\n");
  fprintf(str,"#plt.figure(dpi=300,figsize=(20/2.54,20/2.54))\n");
  fprintf(str,"plt.figure()\n"); 
  fprintf(str,"plt.plot(data1[%d],data1[%d],label=tab1)\n",xcol[0]-1,ycol[0]-1);
  fprintf(str,"plt.scatter(data1[%d],data1[%d],label=tab1)\n",xcol[0]-1,ycol[0]-1);
  fprintf(str,"#plt.errorbar(data1[%d],data1[%d],data1[%d],label=tab1)\n",xcol[0]-1,ycol[0]-1,dycol[0]-1);  
  fprintf(str,"plt.xlabel('X')\n");
  fprintf(str,"plt.ylabel('Y')\n");
  fprintf(str,"#plt.xlim([0,1])\n");
  fprintf(str,"#plt.ylim([0,1])\n");
  fprintf(str,"plt.title('tabplot')\n");
  fprintf(str,"plt.legend()\n");
  fprintf(str,"plt.savefig('pyplot.png')\n");  
  fprintf(str,"plt.show()\n");
}


void pyplot_hist(stream str, string tabname, int *xcol, real *xrange, int bins)
{
  fprintf(str,"tabname = '%s'\n",tabname);
  fprintf(str,"data = np.loadtxt(tabname).T\n");
  fprintf(str,"print(data.shape)\n");
  fprintf(str,"#plt.figure(dpi=300,figsize=(20/2.54,20/2.54))\n");   
  fprintf(str,"plt.figure()\n");
  fprintf(str,"plt.hist(data[%d],%d)\n",xcol[0]-1,bins);
  fprintf(str,"plt.xlabel('X')\n");
  fprintf(str,"plt.ylabel('N')\n");
  fprintf(str,"#plt.xlim([0,1])\n");
  fprintf(str,"#plt.ylim([0,1])\n");
  fprintf(str,"plt.title('tabhist')\n");
  fprintf(str,"plt.savefig('pyplot.png')\n");  
  fprintf(str,"plt.show()\n");
}


