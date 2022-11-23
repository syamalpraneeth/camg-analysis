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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import csv

def csdict(val):
  rtstr= []
  if val == "All": rtstr = "red"
  elif val == "Core": rtstr = 'mediumvioletred'
  elif val == "Interface": rtstr = 'antiquewhite'
  return rtstr;

casedict = {"2GPa":"darkorange", "5GPa":"royalblue"}

def plot_collated(source,ylab,c,name,tag,only):
  tagdict = {"asp":"As Prepared", "ann":"Annealed"}
  b=5
  a=b*len(c)
  fig, ax= plt.subplots(1,1,sharey=True,figsize=(2*b,b))

  fig.subplots_adjust(bottom=0.2)
  if only == 0: csval=['All','Core','Interface']
  else: csval=['All']
  i = 0
  for m in c:
    csind = 0
    for cs in csval:
      for t in tag:
        with open(m+'/tmp.'+source+'_coll_'+name.replace(" ", "_")+'_'+m+'_'+str(t)+'_All','r') as csvfile:
          plots = csv.reader(csvfile, delimiter=',')
          x = []
          l = []
          for row in plots:
            if source=='strain': x.append(float(row[2]))
            else: x.append(float(row[0]))
            l.append(str(row[-2]))
        p=[0,0]
        # plt.sca(ax[i])
        # l.pop(0)              #uncomment to allow only data of 2-7nm clusters be plotted
        # x.pop(0)
        if source == 'pote':
          #print(x,m,s,cs,t)
          ax.plot(l, [-4.9728 for j in range(0,len(l))], '--', color='gray',zorder=5) #x = [i/(-4.979) for i in x]
          ax.text(0.0, -4.9755, 'MG As Prepared', color='gray', fontsize='10')
        if t == tag[0]: ax.plot(l, x, 'o:', c='black', mfc=casedict[m],mec='k',  label=m+' '+tagdict[t],ms=9,zorder=6)
        else: ax.plot(l, x, 's:', c='black', mfc=casedict[m],mec='k',  label=m+' '+tagdict[t],ms=10,zorder=5)
        if 'yl' in locals(): yl = max(max(x),yl)
        else: yl = max(x)
        if 'yl2' in locals(): yl2 = min(min(x),yl2)
        else: yl2 = min(x)
        delta = 0.1*(yl-yl2)
        ylv = yl+delta
        ylv2= yl2-delta
        if i == 0:   ax.set_ylabel(ylab, fontsize=16)
        ax.set_facecolor('white')
        ax.grid(color='k', alpha=0.1, zorder=-2)

        if source != 'vol': ax.set_ylim(ylv2,ylv)
        ax.tick_params(axis="x",direction="in",left='off',labelsize=14)
        ax.tick_params(axis="y",direction="in",labelsize=14)
        plt.setp(ax.spines.values(), linewidth=1.5)
        plt.xticks(rotation=45)
        txy=.94*yl
        if source == 'pote': plt.legend(loc=4,fontsize=12)
        elif source == 'porosity': plt.legend(loc=2,fontsize=12)
        elif source == 'ico-like' or source == 'ico-like': plt.legend(loc=0,fontsize=12)
        else: plt.legend(loc=3,fontsize=12)

        if source == 'strain':
          if m == '5GPa':
            if t == 'asp':
              x0 = []
              x0 = [j for j in x]
            elif t == 'ann':
              xrel = [(i - j)*100 / j for i, j in zip(x, x0)]
              axins1 = inset_axes(ax, "45%",  # width = 30% of parent_bbox
                                  height=1.4,  # height : 1 inch
                                  loc=1,  # bbox_to_anchor=(-0.05,0,1,1),
                                  bbox_transform=ax.transAxes)
              axins1.plot(l, xrel, '*g', ms=10)
              axins1.grid(color='k', alpha=0.1, zorder=-2)
              axins1.set_ylabel('Increase with \n annealing (%)')
              axins1.set_ylim(0.8*min(xrel),1.1*max(xrel))
              for tick in axins1.get_xticklabels():
                tick.set_rotation(45)
    csind+=1
    i+=1

  # for m,a in zip(c,ax):
  #   a.text(0.7, 0.60, m,transform=a.transAxes, fontsize='15')

  fig.tight_layout(rect=[0, 0, 1, 0.95])
  # plt.savefig('collated/'+source+'_coll_'+name.replace(" ", "_")+'_'+tag+'.png',dpi=400)
  plt.savefig('collated/'+source+'_coll_'+name.replace(" ", "_")+'.png',dpi=400)
  plt.clf()
  plt.close()

cases = ['5GPa']
name='Cu50Zr50 random NG seg'
tag = ['asp'] #,'ann']

#plot_collated('strain', r'Atoms with $\eta \geq 2$ (at %)', cases, name, tag, 1)
# plot_collated('strain', r'Atomic readjustment (at % w. $\eta \geq 2$)', cases, name, tag, 1)
# plot_collated('density','Density (units)',cases,name,tag,1)
# plot_collated('porosity','Porosity (% unfilled Vol.)',cases,name,tag,1)
plot_collated('pote','P.E./ Atom (eV/atom)',cases,name,tag,1)
# plot_collated('ico-like','ILO Fraction (at %)',cases,name,tag,1)

# for t in tag:
#   plot_collated('strain','Strain Fraction (at %)',cases,name,t,1)
  # plot_collated('density','Density (units)',cases,name,t,1)

# plot_collated('ico','FI Fraction (at %)',cases,name,tag,1)
# plot_collated('ico-like','ILO Fraction (at %)',cases,name,tag,1)
# plot_collated('nico','Non-ILO Fraction (at %)',cases,name,tag,0)
# plot_collated('vol',r'Volume per Atom ($\AA ^{3}$/atom)',cases,name,tag,0)
# plot_collated('nicovol',r'Free Volume ($\AA ^{3}$/atom)',cases,name,tag,0)
# plot_collated('pote','P.E./ Atom (eV/atom)',cases,name,t,0)
##plot_collated('pote-min',r'$P.E._{Min}$/ Atom (eV/atom)',cases,name,tag,0)
