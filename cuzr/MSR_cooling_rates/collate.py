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
  b=5
  a=b*len(c)
  fig, ax= plt.subplots(1,len(c),sharey=True,figsize=(a,b))

  fig.subplots_adjust(bottom=0.2)
  if only == 0: csval=['All','Core','Interface']
  else: csval=['All']
  i = 0
  # yl = 0
  # yl2 = 0
  # ymx = 0
  # ymi = 0
  for m in c:
    for s in siz:
      csind = 0
      for cs in csval:
        with open(m+'/tmp.'+source+'_coll_'+name.replace(" ", "_")+'_'+m+'_'+tag+'_'+cs,'r') as csvfile:
          plots = csv.reader(csvfile, delimiter=',')
          x = []
          y = []
          l = []
          for row in plots:
            x.append(float(row[0]))
            # y.append(float(row[1]))
            l.append(str(row[-1]))
        # print(x)
        p=[0,0]
        plt.sca(ax[i])
        fmtstr = StrMethodFormatter('{x:0>5.2f}')
        ax[i].yaxis.set_major_formatter(fmtstr)
        # if csind == 0: ax[i].plot(l, x, 'o-', label=s + ' ' + cs,ms=7,zorder=5)
        # if csind == 1: ax[i].plot(l, x, '^-', label=s + ' ' + cs,ms=5,zorder=3)
        # if csind == 2: ax[i].plot(l, x, 's-', label=s + ' ' + cs,ms=5,zorder=4)
        if source == 'vol' and csind==0: ax[i].plot([0,len(x)-1], [19.2,19.2], ':',label='Liquid Melt',color='k',ms=7,zorder=5)
        # print(csdict(cs))
        if csind == 0: ax[i].plot(l, x, 'o:',c='black', mfc=csdict(cs),mec='k',  label=cs,ms=9,zorder=5)
        if csind == 1: ax[i].plot(l, x, '^:', c='black' , mfc=csdict(cs), mec='k', label=cs,ms=7,zorder=3)
        if csind == 2: ax[i].plot(l, x, 's:', c='black', mfc=csdict(cs), mec='k', label=cs,ms=7,zorder=4)
        #ax[i].plot(l,y,'o-', label=s+' Cu-centered')
        # print(yl2)
        if 'yl' in locals(): yl = max(max(x),yl)
        else: yl = max(x)
        if 'yl2' in locals(): yl2 = min(min(x),yl2)
        else: yl2 = min(x)
        delta = 0.1*(yl-yl2)
        ylv = yl+delta
        ylv2= yl2-delta
        # print('min',yl2,'max',yl)
        if i == 0:   ax[i].set_ylabel(ylab, fontsize=16)
        ax[i].set_facecolor('white')
        ax[i].grid(color='k', alpha=0.1, zorder=-2)

        # ax[i].set_ylabel('.', color=(0, 0, 0, 0),labelpad=1)
        # ax[i].set_ylabel('.', color=(0, 0, 0, 0),labelpad=1)
        if source != 'vol': ax[i].set_ylim(ylv2,ylv)
        ax[i].tick_params(axis="x",direction="in",left='off',labelsize=14)
        ax[i].tick_params(axis="y",direction="in",labelsize=14)
        # ax[i].grid(color='k', zorder=-2)
        plt.setp(ax[i].spines.values(), linewidth=1.5)
        plt.xticks(rotation=45) #,weight='bold')
        # txy=(yl+ylv2)/2
        txy=.94*yl
        # if source =='vol': ax[i].text(2, 19, 'Volume in Liquid State', fontsize='13')
        # if cs =='All': ax[i].text(4, txy, "{:.0e}".format(float(m)) + r'$\frac{K}{s}$', weight='bold', fontsize='13')
        if i ==0 :
          if source != 'vol' and source != 'nicovol':
            if source != 'ico-like' in source != 'ico-like': plt.legend(loc=0,fontsize=13)
            else: plt.legend(loc=3,fontsize=13)
          else: plt.legend(loc=2,fontsize=13)
        csind+=1
        # del yl,yl2
    i+=1

  for m,a in zip(c,ax):
    qr = m.split("e")
    # print(qr,txy)
    # if cs == 'All':
    a.text(0.7, 0.7, r'${10}^{' + str(qr[1]) + '}$' + r' $K/s$',transform=a.transAxes, fontsize='15')

  fig.tight_layout(rect=[0, 0, 1, 0.95])
  plt.savefig('collated/'+source+'_coll_'+name.replace(" ", "_")+'.png',dpi=400)
  plt.clf()
  plt.close()

cases = ['1e10','1e12','1e14']
siz = ['3nm']
name=siz[0]+' Cu50Zr50 CAMG vs MG'
tag = 'asp'

plot_collated('ico','FI Fraction (at %)',cases,name,tag,0)
plot_collated('ico-like','ILO Fraction (at %)',cases,name,tag,0)
plot_collated('nico','Non-ILO Fraction (at %)',cases,name,tag,0)
plot_collated('vol',r'Volume per Atom ($\AA ^{3}$/atom)',cases,name,tag,0)
plot_collated('nicovol',r'Free Volume ($\AA ^{3}$/atom)',cases,name,tag,0)
plot_collated('pote','P.E./ Atom (eV/atom)',cases,name,tag,0)
# plot_collated('pote-min',r'$P.E._{Min}$/ Atom (eV/atom)',cases,name,tag,0)
