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
  return str;
}

void pyplot_close(stream str)
{
  strclose(str);
}

void pyplot_plot(stream str, string tabname, int *xcol, int *ycol, int *dycol, real *xrange, real *yrange)
{
  fprintf(str,"tabname = '%s'\n",tabname);
  fprintf(str,"data = np.loadtxt(tabname).T\n");
  fprintf(str,"print(data.shape)\n");
  fprintf(str,"plt.figure()\n");
  fprintf(str,"# plt.plot(data[%d],data[%d])\n",   xcol[0]-1,ycol[0]-1);
  fprintf(str,"plt.scatter(data[%d],data[%d])\n",xcol[0]-1,ycol[0]-1);
  fprintf(str,"plt.show()\n");
}


void pyplot_hist(stream str, string tabname, int *xcol, real *xrange, int bins)
{
  fprintf(str,"tabname = '%s'\n",tabname);
  fprintf(str,"data = np.loadtxt(tabname).T\n");
  fprintf(str,"print(data.shape)\n");
  fprintf(str,"plt.figure()\n");
  fprintf(str,"plt.hist(data[%d],%d)\n",xcol[0]-1,bins);
  fprintf(str,"plt.show()\n");
}


