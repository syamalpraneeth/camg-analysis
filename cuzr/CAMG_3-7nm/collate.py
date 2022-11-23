#!/usr/bin/env python3

##################################################################
##								##                               
##								##
##Dependencies	:						##
##Influences	:						##
##################################################################
## ver.	: 2019--, Syamal Praneeth Chilakalapudi, KIT, INT	##                            
##Author Email    :syamalpraneeth@gmail.com			##
##################################################################

import os, sys
from ovito.io import *
from ovito.modifiers import *
import numpy
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import csv

def csdict(val):
  rtstr= []
  if val == "All": rtstr = "blue"
  elif val == "Core": rtstr = 'mediumvioletred'
  elif val == "Interface": rtstr = 'antiquewhite'
  return rtstr;

def plot_collated(source,ylab,c,name,tag,only):
  if only == 0: csval=['All','Core','Interface']
  else: csval=['All']

  b=5
  a=b*len(csval)
  fig, ax= plt.subplots(1,len(csval),sharey=True,figsize=(a,b))
  fig.subplots_adjust(bottom=0.2)

  markr = ['o','s']
  mks = [9,12]
  for m in c:
    for s,mkrsize,mkrshape in zip(siz,mks,markr):
      csind = 0
      # if isinstance(ax, list):
      for axs,cs in zip(ax,csval):
      # else: for cs in csval:
      # for axs, cs in zip([ax], csval):
        # if type(ax) != list: axs=ax
        with open(s+'/tmp.'+source+'_coll_'+s+'_'+name.replace(" ", "_")+'_'+m+'_'+tag+'_'+cs,'r') as csvfile:
          plots = csv.reader(csvfile, delimiter=',')
          x, y, l = [], [], []
          for row in plots:
            x.append(float(row[0]))
            l.append(str(row[-1]))
        p=[0,0]
        plt.sca(axs)
        fmtstr = StrMethodFormatter('{x:0>5.2f}')
        axs.yaxis.set_major_formatter(fmtstr)

        axs.plot(l, x, mkrshape+':', c='black', mfc=csdict(cs), mec='k', label=s, ms=mkrsize, zorder=5)

        if 'yl' in locals(): yl = max(max(x),yl)
        else: yl = max(x)
        if 'yl2' in locals(): yl2 = min(min(x),yl2)
        else: yl2 = min(x)
        delta = 0.1*(yl-yl2)
        ylv = yl+delta
        ylv2= yl2-delta

        if axs == ax[0]:   axs.set_ylabel(ylab, fontsize=18)
        axs.set_facecolor('white')
        axs.grid(color='k', alpha=0.1, zorder=-2)

        if source != 'vol': axs.set_ylim(ylv2,ylv)
        axs.tick_params(axis="x",direction="in",left='off',labelsize=14)
        axs.tick_params(axis="y",direction="in",labelsize=14)
        # axs.legend(fontsize=12, bbox_to_anchor=(0.6, -0.02, 0.4, 1), bbox_transform=axs.transAxes)
        axs.legend(fontsize=12,loc='lower left')

        if cs!='All': axs.text(0.7, 0.7, cs, transform=axs.transAxes, fontsize='15')
        else: axs.text(0.7, 0.65, r"Entire"+"\n"+"sample", transform=axs.transAxes, fontsize='15')

        plt.setp(axs.spines.values(), linewidth=1.5)
        plt.xticks(rotation=45)  # ,weight='bold')

        txy=.94*yl
        csind+=1

  fig.tight_layout(rect=[0, 0, 1, 0.95])

  plt.savefig('collated/'+source+'_coll_'+name.replace(" ", "_")+'.png',dpi=400)
  plt.clf()
  plt.close()

cases = ['1e10']
siz = ['3nm','7nm']
# name=siz[0]+'Cu50Zr50 CAMG vs MG'
name='Cu50Zr50 3nm vs 7nm'
tag = 'asp'

plot_collated('ico','FI Fraction (at %)',cases,name,tag,0)
plot_collated('ico-like','ILO Fraction (at %)',cases,name,tag,0)
plot_collated('nico','Non-ILO Fraction (at %)',cases,name,tag,0)
plot_collated('vol',r'Volume per Atom ($\AA ^{3}$/atom)',cases,name,tag,0)
plot_collated('nicovol',r'Free Volume ($\AA ^{3}$/atom)',cases,name,tag,0)
plot_collated('pote','P.E./ Atom (eV/atom)',cases,name,tag,0)
plot_collated('pote-min',r'$P.E._{Min}$/ Atom (eV/atom)',cases,name,tag,0)
