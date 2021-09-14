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

import os, csv
import numpy
import pandas as pd
import math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
#from plotting import *

def cscheme(c):
  if c==1: ret='black'
  elif c==2: ret='blue'
  elif c==3: ret='green'
  elif c==4: ret='darkorange'
  elif c==5: ret='red'
  elif c==6: ret='crimson'
  elif c==7: ret='white'
  # elif c==7: ret='lightcyan'
  elif c==8: ret='mediumvioletred'
  elif c==9: ret='antiquewhite'
  else: raise ValueError('Color value between 1-6 only')
  return ret;

class datclass:
    def __init__(self):
        self.a=[1,2]
        self.b=[2,3]
        self.c=[3,4]
        self.d=1
        self.h=[0,0,0,0]
        self.v=1

def sort_voro(ind_mat,vol):
  a = pd.read_csv('ico.txt',header=None)
  b = pd.read_csv('crys.txt',header=None)
  c = pd.read_csv('mixed.txt',header=None)

  cnt_ico=0
  cnt_crys=0
  cnt_mix=0
  cnt_etc=0    
  itr = 0

  ico=[]
  nico=[]

  for r in ind_mat:
    # print(r.tolist())
    indval = r.tolist()
    # print(sum(indval[6:]))
    if sum(indval[6:])==0:   #setting number of 1,2, and >6 edged faces to zero
      if r.tolist()[2:6] in a.values.tolist():
        cnt_ico+=1
        ico.append(vol[itr])
      elif r.tolist()[2:6] in b.values.tolist():
        cnt_crys+=1
        nico.append(vol[itr])
      elif r.tolist()[2:6] in c.values.tolist():
        cnt_mix+=1
        nico.append(vol[itr])
      else:
          cnt_etc+=1
          nico.append(vol[itr])
    else:
      cnt_etc+=1
      nico.append(vol[itr])
    itr+=1
  ht = [cnt_ico, cnt_crys, cnt_mix, cnt_etc]
  ht_p = [x*100/len(ind_mat) for x in ht]

  print(numpy.shape(nico),'non ICO in',numpy.shape(ind_mat),'not-unique indices')
  # print('heights',ht)
  # print('Sum hts',sum(ht_p))
  # print(ht_p)
  return ht_p,nico;

def voro_plot_sort(name,m1,m2,m3,l1,l2,l3,species,cs,typ):
    h1 = m1.h
    h2 = m2.h
    h3 = m3.h

    if species==0:
        case='All'
    elif species==1 or species ==3:
        case='Cu'
    elif species==2 or species ==4:
        case='Zr'
    else:
        raise ValueError('species value only 0,1 or 2!')
    if typ==1:
        comm='As Prepared'
        tag='asp'
    elif typ==2:
        comm='Annealed'
        tag='ann'
    else:
        raise ValueError('typ value only 1 or 2!')
    if cs==0:
        csstr='All'
    elif cs==1:
        csstr='Core'
    elif cs==2:
        csstr='Interface'
    else:
        raise ValueError('cs value only 0, 1 or 2!')

    fig, (ax,ax2,ax3) = plt.subplots(3,1, figsize=(5,5),sharex=True)
    ax.set_facecolor(cscheme(cs + 7))
    ax2.set_facecolor(cscheme(cs + 7))
    ax3.set_facecolor(cscheme(cs + 7))

    if cs == 1:
      ax.patch.set_alpha(0.2)
      ax2.patch.set_alpha(0.2)
      ax3.patch.set_alpha(0.2)
    # plt.grid(color='m', alpha=0.1, zorder=-2)
    b=['ICO-Like','Crys-Like','Mixed','Other']
    w=0.1
    y1 = numpy.arange(len(b))+1.2

    fig.text(0.03,0.32,'Atomic Fraction (at %)', fontsize='16',rotation='vertical')
    plt.ylabel(' ', fontsize='22')
    #plt.ylabel('Atomic Fraction (at %)', fontsize='15')
    plt.xlabel('Voronoi Index Type', fontsize='16')

    for axs in [ax,ax2,ax3]:
      plt.sca(axs)
      if cs == 0:
        axs.grid(color='k', alpha=0.1, zorder=-2)
      else:
        axs.grid(color='w', zorder=-2)
      # axs.grid(color='w', zorder=-2)
      plt.setp(axs.spines.values(), linewidth=1.5)
      p1 = axs.bar(y1, h1, width=0.2, label=l1, color=cscheme(1), edgecolor='black', zorder=3)
      p2 = axs.bar(y1 + w, h2, width=0.2, label=l2, color=cscheme(2), edgecolor='black', zorder=3)
      p3 = axs.bar(y1 + 2 * w, h3, width=0.2, label=l3, color=cscheme(3), edgecolor='black', hatch=".", zorder=3)

      axs.yaxis.set_label_position('left')
      axs.tick_params(axis="y", direction="in", labelsize=14)
    ax.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax.tick_params(axis="x", direction="in", left='on', bottom=False, labelsize=14)
    ax2.tick_params(axis="x", direction="in", left='on', bottom=False ,labelsize=14)
    ax3.tick_params(axis="x", direction="in", left='on', bottom=True, labelsize=14)

    plt.xticks(y1+2.5*w, b,rotation=0) #,fontsize=25)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])#

    qr = name.split(" ")[-1].split("e")
    if species == 0:
      ax3.legend(fontsize=12, loc=10,
                 bbox_to_anchor=(0.02,1.45,1,1))
    else:
      ax3.legend(fontsize=12, title=case,
                bbox_to_anchor=(0.02, 1.45, 1, 1))

    # yl = max(max(h1, h2, h3, h4, h5, h6))*1.10
    yl = (max(h1[3], h2[3], h3[3])) * 1.05
    ylmin = (min(h1[3], h2[3], h3[3])) * 0.97
    yl3 = (max(h1[2], h2[2], h3[2])) * 1.10
    yl2 = (max(h1[0], h2[0], h3[0])) * 1.05
    yl2min = (min(h1[0], h2[0], h3[0])) * 0.95

    d = .015
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal
    ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    kwargs.update(transform=ax3.transAxes)  # switch to the bottom axes
    ax3.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax3.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    # print('YMAX', yl)
    ax3.set_ylim(0,yl3) #10
    ax2.set_ylim(0.9*yl2min, yl2)
    ax.set_ylim(ylmin,yl)
    #ax.set_ylim([0,90])
    # ax.set_ylabel(' ',fontsize=18)
    # plt.subplots_adjust(bottom=0.4, top=0.8)
    plt.subplots_adjust(hspace=0.05)

    if cs==0: plt.savefig('voronoi-sort_'+ name.replace(" ", "_")+'_'+tag+'_'+case+'.png',dpi=400)
    elif cs>0: plt.savefig('voronoi-sort_'+ name.replace(" ", "_")+'_'+tag+'_'+case+'_'+csstr+'.png',dpi=400)

    plt.clf()
    plt.close()
    
def plot_atmhisto_nico(name,m1,m2,m3,l1,l2,l3,species,cs,typ):
    h1,n1 = m1.v,m1.d
    h2,n2 = m2.v,m2.d
    h3,n3 = m3.v,m3.d

    if species==0: case='All'
    elif species==1 or species ==3: case='Cu'
    elif species==2 or species ==4: case='Zr'
    else: raise ValueError('species value only 0,1 or 2!')
    if cs==0: csstr='All'
    elif cs==1: csstr='Core'
    elif cs==2: csstr='Interface'
    else: raise ValueError('cs value only 0, 1 or 2!')
    if typ==1:
        comm='As Prepared'
        tag='asp'
    elif typ==2:
        comm='Annealed'
        tag='ann'
    else: raise ValueError('typ value only 1 or 2!')

    bin1=math.floor(n1/50)
    bin2=math.floor(n2/50)
    bin3=math.floor(n3/50)

    fig,ax= plt.subplots(figsize=(10,5),sharey=True)
    if cs == 0:
      plt.grid(color='k', alpha=0.1, zorder=-2)
    else:
      plt.grid(color='w', zorder=-2)
    # plt.grid(color='w', alpha=0.1, zorder=-2)
    plt.setp(ax.spines.values(), linewidth=1.5)
    # plt.grid(color='m', alpha=0.1, zorder=-2)
    plt.setp(ax.spines.values(), linewidth=1.5)
    plt.sca(ax)

    bin1 = bin2 = bin3 = bin4 = bin5 = bin6 = 100

    ht1,b1 = numpy.histogram(h1,bins=bin1)
    ht2,b2 = numpy.histogram(h2,bins=bin1)
    ht3,b3 = numpy.histogram(h3,bins=bin1)

    a1,c1=ht1/n1,(b1[:-1] + b1[1:])/2
    a2,c2=ht2/n2,(b2[:-1] + b2[1:])/2
    a3,c3=ht3/n3,(b3[:-1] + b3[1:])/2

    # print('Dis tot', ht1.sum(axis=0)/n1,ht2.sum(axis=0)/n2,ht3.sum(axis=0)/n3,ht4.sum(axis=0)/n4,ht5.sum(axis=0)/n5,ht6.sum(axis=0)/n6)

    plt.plot(c1, a1, color=cscheme(1), label=l1, zorder=3, lw=2)
    plt.plot(c2, a2, color=cscheme(2), label=l2, zorder=3, lw=2)
    plt.plot(c3, a3, color=cscheme(3), label=l3, zorder=3, lw=2)

    qr = name.split(" ")[-1].split("e")
    # ax.legend(fontsize=12, title=r'${10}^{' + str(qr[1]) + '}$' + r'$K/s$' + r'$\frac{K}{s}$')
    ax.legend(fontsize=12)


    ax.set_xlabel('Volume per Atom (Units)', fontsize=18)
    ax.set_ylabel('Normalised Counts ', fontsize=18)
    ax.yaxis.set_label_position('left')
    ax.tick_params(axis="both", direction="in", bottom=True, top=True, left=True, right=True, labelsize=15)
    plt.subplots_adjust(bottom=0.4, top=0.8)

    y1 = max(max(a1), max(a2), max(a3))
    y2 = min(min(a1), min(a2), min(a3))
    x1 = max(max(c1), max(c2), max(c3))
    x2 = min(min(c1), min(c2), min(c3))
    # print('Max Y',yl)
    yl1 = y1 * 1.05
    yl2 = y2 * 0.90
    xl1 = x1 * 1.05
    xl2 = x2 * 0.90
    ax.set_ylim(yl2, yl1)
    ax.set_xlim(xl2, xl1)

    plt.text(15, 0.8 * yl1, r'Cu', weight='bold', fontsize='13')
    plt.text(23, 0.6 * yl1, r'Zr', weight='bold', fontsize='13')

    # axins1 = inset_axes(ax,
    #             width="40%", # width = 30% of parent_bbox
    #             height=1.5, # height : 1 inch
    #             loc=1)
    #
    # axins2 = inset_axes(ax,
    #                 width="30%", # width = 30% of parent_bbox
    #                 height=1.2, # height : 1 inch
    #                 loc=4)

    # if species==0:
    #     x1, x2, y1, y2 = 12,24,2.5e-2,3.8e-2
    #     axins1.set_xlim(x1, x2)
    #     axins1.set_ylim(y1, y2)
    # if species==1 or species ==3:
    #     x1, x2, y1, y2 = 12.5,15.5,1.8e-2,3.5e-2
    #     axins1.set_xlim(x1, x2)
    #     axins1.set_ylim(y1, y2)
    # elif species==2 or species == 4:
    #     x1, x2, y1, y2 = 20.5,22.5,2.5e-2,4.2e-2
    #     axins1.set_xlim(x1, x2)
    #     axins1.set_ylim(y1, y2)
    # axins2.xaxis.tick_top()
    # axins2.xaxis.set_label_position('top')
    # axins2.tick_params('x',labelrotation=45)

    rng=[l1,l2,l3]

    # axins1.plot(c1,a1, color='black', label=l1)
    # axins1.plot(c2,a2, color='pink', label=l2)
    # axins1.plot(c3,a3, color='lime', label=l3)
    # axins1.plot(c4,a4, color='cyan', label=l4)
    # axins1.plot(c5,a5, color='crimson', label=l5)
    # axins1.plot(c6,a6, color='blue', label=l6)

    #v1=ht1*c1/n1
    #v2=ht2*c2/n2
    #v3=ht3*c3/n3
    #v4=ht4*c4/n4
    #v5=ht5*c5/n5
    #v6=ht6*c6/n6
    v1 = numpy.array(h1)
    v2 = numpy.array(h2)
    v3 = numpy.array(h3)

    tot_v = [v1.sum(axis=0)/n1,v2.sum(axis=0)/n2,v3.sum(axis=0)/n3]
    # tot_v = [h1.sum(axis=0)/n1,h2.sum(axis=0)/n2,h3.sum(axis=0)/n3,h4.sum(axis=0)/n4,h5.sum(axis=0)/n5,h6.sum(axis=0)/n6]

    # axins2.plot(rng,tot_v,label=r'$\sum^{N_{nonICO}}(V_{atom})$')
    # axins2.legend()
    #
    # axins1.set_yscale("log")
    # axins1.set_xscale("log")
    # mark_inset(ax, axins1, loc1=2, loc2=4, fc="none", ec="0.5")

    fig.tight_layout(rect=[0, 0.05, 1, 0.96])
    if cs==0: plt.savefig('vol-atom-nico_'+ name.replace(" ", "_")+'_'+tag+'_'+case+'.png',dpi=400)
    elif cs>0: plt.savefig('vol-atom-nico_'+ name.replace(" ", "_")+'_'+tag+'_'+case+'_'+csstr+'.png',dpi=400)

    plt.clf()
    plt.close()
    
    if species == 0:
      rows=zip(tot_v,rng)
      # fil = 'tmp.nicovol_coll_' + name.replace(" ", "_") + '_' + tag + '_' + case + '_' + csstr
      fil = 'tmp.nicovol_coll_' + name.replace(" ", "_") + '_' + tag + '_' + csstr
      with open(fil, 'w') as f:
        writer = csv.writer(f)
        for row in rows:
          writer.writerow(row)
