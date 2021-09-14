#!/usr/bin/env python3

import numpy, csv
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from analysis import *
from sort_voro import *

def cscheme(c):
  if c==1: ret='black'
  elif c==2: ret='blue'
  elif c==3: ret='green'
  elif c==4: ret='darkorange'
  elif c==5: ret='red'
  elif c==6: ret='crimson'
  # elif c==7: ret='lightcyan'
  elif c==7: ret='white'
  elif c==8: ret='mediumvioletred'
  elif c==9: ret='antiquewhite'
  else: raise ValueError('Color value between 1-9 only')
  return ret;

def plot_rdf(name, case, delsub, cs, typ):
  if typ == 1:
    comm = 'As Prepared'
    tag = 'asp'
  elif typ == 2:
    comm = 'Annealed'
    tag = 'ann'
  else:
    raise ValueError('typ value only 1 or 2!')
  if cs == 0:
    csstr = 'All'
  elif cs == 1:
    csstr = 'Core'
  elif cs == 2:
    csstr = 'Interface'
  else:
    raise ValueError('cs value only 0, 1 or 2!')
  
  i = 0
  x, y1, y2, y3 = [], [], [], []
  fil = open('tmp.rdf_' + case+'_'+csstr, 'r')
  for line in fil:
    values = [float(s) for s in line.split()]
    x.append(values[0])
    if delsub == 1:
      y1.append(values[-3])
      y2.append(values[-2])
      y3.append(values[-1])
    else:
      y1.append(values[1])
      y2.append(values[2])
      y3.append(values[3])

  fil.close()

  fig, ax = plt.subplots(figsize=(10, 5))  # ,sharey=True)
  ax.set_facecolor(cscheme(cs + 7))
  if cs == 1: ax.patch.set_alpha(0.2)
  if cs == 0:
    ax.grid(color='k',alpha=0.1, zorder=-2)
  else:
    ax.grid(color='w',zorder=-2)
  # if delsub == 0:
  #   fig.text(0.4, 0.94, 'RDF ' + name[:-4] + ' (' + comm + ') :' + case, ha='center', fontsize=15)
  # elif delsub == 1:
  #   fig.text(0.4, 0.94, 'RDF ' + name[:-4] + ' (' + comm + ') :' + case, ha='center', fontsize=15)
  #   if cs == 0:
  #     fig.text(0.75, 0.94, 'All atoms', ha='left', fontsize=15)
  #   elif cs == 1:
  #     fig.text(0.76, 0.94, ' ' + csstr + ' atoms', ha='left', color="red", fontsize=15)
  #   elif cs == 2:
  #     fig.text(0.76, 0.94, ' ' + csstr + ' atoms', ha='left', color="blue", fontsize=15)
  # fig.text(0.1,0.85,"{:.0e}".format(float(name.split(" ")[-1]))+r'$\frac{K}{s}$',weight='bold',fontsize='15')

  # ax.grid(color='m', alpha=0.1, zorder=-2)
  plt.setp(ax.spines.values(), linewidth=1.5)

  plt.xlabel(r'Pair Separation Distance ($\AA$)', fontsize=15)
  plt.ylabel('g(r)', fontsize=15)

  p1 = ax.plot(x, y1, label='Cu-Cu', zorder=3, lw=2)
  p2 = ax.plot(x, y2, label='Cu-Zr', zorder=3, lw=2)
  p3 = ax.plot(x, y3, label='Zr-Zr', zorder=3, lw=2)

  ax.tick_params(axis="both", direction="in", bottom=True, top=True, left=True, right=True, labelsize=14)
  fig.tight_layout(rect=[0, 0.03, 1, 0.95])

  # qr = name.split(" ")[-1].split("e")
  # ax.legend(fontsize=12, title= r'${10}^{' + str(qr[1]) + '}$' + r'$K/s$'+','+case)
  ax.legend(fontsize=12)
  # ax.legend(fontsize=12, title="{:.0e}".format(float(name.split(" ")[-1])) + r'$\frac{K}{s}$'+','+case)

  if cs == 0:
    plt.savefig('rdf_' + name.replace(" ", "_") + '_' + case + '_' + tag + '.png',dpi=400)
  elif cs > 0:
    plt.savefig('rdf_' + name.replace(" ", "_") + '_' + case + '_' + tag + '_' + csstr + '.png',dpi=400)

  plt.close()

def plot_allrdf(name, case1, case2, case3, case4, cs, typ):
  if typ == 1:
    comm = 'As Prepared'
    tag = 'asp'
  elif typ == 2:
    comm = 'Annealed'
    tag = 'ann'
  else:
    raise ValueError('typ value only 1 or 2!')
  if cs == 0:
    csstr = 'All'
  elif cs == 1:
    csstr = 'Core'
  elif cs == 2:
    csstr = 'Interface'
  else:
    raise ValueError('cs value only 0, 1 or 2!')

  fig , axs = plt.subplots(1,3, figsize=(15, 5))
  # fig = plt.figure(figsize=plt.figaspect(0.33))
  # axs = [fig.add_subplot(1,3,1, projection='3d'), fig.add_subplot(1,3,2, projection='3d'), fig.add_subplot(1,3,3, projection='3d')]
  for ax in axs:
    ax.set_facecolor(cscheme(cs+7))
    if cs == 1: ax.patch.set_alpha(0.2)

  # plt.grid(color='w', zorder=-2)
  for ax in axs:
    # ax.grid(color='k', alpha=0.1, zorder=-2)
    if cs == 0:
      ax.grid(color='k', alpha=0.1, zorder=-2)
    else:
      ax.grid(color='w', zorder=-2)
    # ax.grid(color='w', zorder=-2)
    plt.setp(ax.spines.values(), linewidth=1.5)

  w=0.5
  i = 0
  for case in [case1, case2, case3, case4]:
    x, y1, y2, y3 = [], [], [], []
    fil = open('tmp.rdf_' + case + '_' + csstr, 'r')

    for line in fil:
      values = [float(s) for s in line.split()]
      x.append(values[0])

      if 'MG' in case or case == 'NG':
        y1.append(values[1])
        y2.append(values[2])
        y3.append(values[3])
      else:
        y1.append(values[-3])
        y2.append(values[-2])
        y3.append(values[-1])

    x = numpy.array(x)
    y1 = numpy.array(y1)
    y2 = numpy.array(y2)
    y3 = numpy.array(y3)

    i+=1
    p1 = axs[0].plot(x, y1+(4-i)*w, label='_Hidden', color=cscheme(i), zorder=3, lw=2)
    p2 = axs[1].plot(x, y2+(4-i)*w, label='_Hidden', color=cscheme(i), zorder=3, lw=2)
    p3 = axs[2].plot(x, y3+(4-i)*w, label='_Hidden', color=cscheme(i), zorder=3, lw=2)
    # p1 = axs[0].plot( x, y1, (6-i)*w, label='_Hidden', color=cscheme(i), fillstyle ='full', zorder=3, lw=1)
    # p2 = axs[1].plot(x, y2, (6-i)*w, label='_Hidden', color=cscheme(i), zorder=3, lw=1)
    # p3 = axs[2].plot(x, y3, (6-i)*w, label='_Hidden', color=cscheme(i), zorder=3, lw=1)

    for ax in axs:
      ax.text(0, 0.08+(4 - i)*w, case, fontsize='10')
      # ax.view_init(azim=275, elev=100)
      ax.tick_params(axis="x",which='minor',bottom=True, top=True, labelsize=14)

  axs[0].tick_params(axis="both", direction="in", bottom=True, top=True, left=True, right=False, labelsize=15)
  axs[1].tick_params(axis="both", direction="in", bottom=True, top=True, left=False, right=False, labelsize=15)
  axs[2].tick_params(axis="both", direction="in", bottom=True, top=True, left=False, right=True, labelsize=15)
  fil.close()

  axs[1].set_xlabel(r'Pair Separation Distance ($\AA$)', fontsize=17)
  axs[0].set_ylabel('g(r) (relative units)', fontsize=17)

  fig.tight_layout(rect=[0, 0.03, 1, 0.95])

  qr = name.split(" ")[-1].split("e")
  # ax.legend(fontsize=12, title= r'${10}^{' + str(qr[1]) + '}$' + r'$K/s$'+','+case)
  # axs[0].legend(fontsize=14, title= r'${10}^{' + str(qr[1]) + '}$' + r'$K/s$' + ', ' + 'Cu-Cu')
  # axs[1].legend(fontsize=14, title= r'${10}^{' + str(qr[1]) + '}$' + r'$K/s$' + ', ' + 'Cu-Zr')
  # axs[2].legend(fontsize=14, title= r'${10}^{' + str(qr[1]) + '}$' + r'$K/s$' + ', ' + 'Zr-Zr')
  axs[0].legend(fontsize=14, title='Cu-Cu')
  axs[1].legend(fontsize=14, title= 'Cu-Zr')
  axs[2].legend(fontsize=14, title= 'Zr-Zr')

  plt.savefig('rdf_' + name.replace(" ", "_") +  '_' + tag + '_' + csstr + '.png', dpi=400)
  plt.close()

def pe_plot_species(name, c1, c2, c3, c4, cs, typ):
  if cs == 0:
    csstr = 'All'
  elif cs == 1:
    csstr = 'Core'
  elif cs == 2:
    csstr = 'Interface'
  else:
    raise ValueError('cs value only 0, 1 or 2!')

  #print("works pe plot 1")

  i = 0
  rng = [c1, c2, c3, c4]
  for c in rng:
    x1, y1 = [], []
    if cs == 0: fil = open('tmp.histo_' + c, 'r')
    if cs > 0: fil = open('tmp.histo_' + c + '_' + csstr, 'r')
    next(fil)
    next(fil)
    for line in fil:
      values = [float(s) for s in line.split()]
      x1.append(values[0])
      y1.append(values[1])
    fil.close()
    if i == 0:
      x = numpy.zeros((len(x1), len(rng)))
      y = numpy.zeros((len(y1), len(rng)))
    x[:, i] = numpy.array(x1) #.reshape(len(x1))
    y[:, i] = numpy.array(y1) #.reshape(len(y1))
    i += 1

  #print("works pe plot 2")

  tot_atms = y.sum(axis=0)
  #print('totalatoms',tot_atms)
  y_norm = y / tot_atms
  pot_e = (numpy.multiply(x, y_norm).sum(axis=0))

  if typ == 1:
    comm = 'As Prepared'
    tag = 'asp'
  elif typ == 2:
    comm = 'Annealed'
    tag = 'ann'
  else:
    raise ValueError('typ value only 1 or 2!')

  #print("works pe plot 2.5")

  fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
  for axis in [ax,ax2]:
    axis.set_facecolor(cscheme(cs + 7))
    if cs == 1: axis.patch.set_alpha(0.2)
    if cs == 0:
      axis.grid(color='k', alpha=0.1, zorder=-2)
    else:
      axis.grid(color='w', zorder=-2)
    # axis.grid(color='w',zorder=-2)
    plt.setp(axis.spines.values(), linewidth=1.5)

    plt.sca(axis)
    p1 = axis.plot(x[:, 0], y_norm[:, 0], color=cscheme(1), label=c1,zorder=3, lw=2)
    p2 = axis.plot(x[:, 1], y_norm[:, 1], color=cscheme(2), label=c2,zorder=3, lw=2)
    p3 = axis.plot(x[:, 2], y_norm[:, 2], color=cscheme(3), label=c3,zorder=3, lw=2)
    p4 = axis.plot(x[:, 3], y_norm[:, 3], color=cscheme(4), label=c4,zorder=3, lw=2)

  # qr = name.split(" ")[-1].split("e")
  # ax.legend(fontsize=12, bbox_to_anchor=(1, 0.6), loc='upper right', title= r'${10}^{' + str(qr[1]) + '}$' + r'$K/s$')
  ax.legend(fontsize=12, bbox_to_anchor=(1, 0.6))

  ax.set_xlabel('P.E. per Zr Atom (eV)', fontsize=17)
  ax2.set_xlabel('P.E. per Cu Atom (eV)', fontsize=17)
  ax.set_ylabel('Normalised Counts', fontsize=17)
  ax.yaxis.set_label_position('left')
  ax.tick_params(axis="both", direction="in", bottom=True, top=True, left=True, right=False, labelsize=16)
  ax2.tick_params(axis="both", direction="in", bottom=True, top=True, left=False, right=True, labelsize=16)

  ax.set_ylim(0.0, 0.1)
  # ax2.set_ylim(0.0,1.1)
  ax.set_xlim(-7, -5)  # outliers only
  ax2.set_xlim(-4, -2)  # most of the data
  ax.spines['right'].set_visible(False)
  ax2.spines['left'].set_visible(False)

  axins1 = inset_axes(ax,
                      width="50%",  # width = 30% of parent_bbox
                      height=1.,  # height : 1 inch
                      loc=1)
  axins2 = inset_axes(ax2,
                      width="50%",  # width = 30% of parent_bbox
                      height=1.,  # height : 1 inch
                      loc=1)
  axins3 = inset_axes(ax2,
                      width="40%",  # width = 30% of parent_bbox
                      height=1.5,  # height : 1 inch
                      bbox_to_anchor=(0.0771,0.11,0.88,0.3),
                      bbox_transform=ax2.transAxes,
                      loc=4)

  for axis in [axins1,axins2]:
    axis.plot(x[:, 0], y_norm[:, 0], color=cscheme(1), label=c1)
    axis.plot(x[:, 1], y_norm[:, 1], color=cscheme(2), label=c2)
    axis.plot(x[:, 2], y_norm[:, 2], color=cscheme(3), label=c3)
    axis.plot(x[:, 3], y_norm[:, 3], color=cscheme(4), label=c4)

  # axins3.plot(rng, pot_e, label=r'$\sum^{N_{norm}}(E^{pot})$')
  fig.tight_layout()
  plt.subplots_adjust(wspace=0.15)

  d = .015
  kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
  ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
  ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
  kwargs.update(transform=ax2.transAxes)
  ax2.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal
  ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal

  x1, x2, y1, y2 = -6.5, -6.3, 0.058, 0.084
  axins1.set_xlim(x1, x2)
  axins1.set_ylim(y1, y2)
  x1, x2, y1, y2 = -3.6, -3.4, 0.058, 0.084
  axins2.set_xlim(x1, x2)
  axins2.set_ylim(y1, y2)
  mark_inset(ax, axins1, loc1=2, loc2=4, fc="none", ec="0.5",zorder=3)
  mark_inset(ax2, axins2, loc1=2, loc2=4, fc="none", ec="0.5",zorder=3)

  x1, y1 = [], []
  fil2 = open('tmp.pote_coll_Cu50Zr50_MG_asp_' + csstr, 'r')
  for line in fil2:
    values = [s for s in line.split(',')]
    y1.append(float(values[0]))
    x1.append(values[1].split('\n')[0])

  print(x1)
  print(y1)
  axins3.plot([1,2,3,4], y1, 'o-', label=r'$\sum^{N_{norm}}(E^{pot})$')
  axins3.set_ylabel(r'Avg. P.E. (eV/atom)',fontsize=12)
  plt.xticks([1,2,3,4], x1, rotation=45)

  if cs == 0:
    plt.savefig('pe-atom_' + name.replace(" ", "_") + '_' + tag + '.png',dpi=400)
  elif cs > 0:
    plt.savefig('pe-atom_' + name.replace(" ", "_") + '_' + tag + '_' + csstr + '.png',dpi=400)

  plt.clf()
  plt.close()

  # axins3.plot(rng, pot_e, label=r'$\sum^{N_{norm}}(E^{pot})$')
  # rows=zip(pot_e,rng)
  # fil = 'tmp.pote_coll_' + name.replace(" ", "_") + '_' + tag + '_' + csstr
  # with open(fil, 'w') as f:
  #   writer = csv.writer(f)
  #   for row in rows:
  #     writer.writerow(row)

def voro_plot_species(name, m1, m2, m3, m4, l1, l2, l3, l4, species, cs, typ):
  voro_plot_sort(name, m1, m2, m3, m4, l1, l2, l3, l4, species, cs, typ)
  h1, b1 = m1.a, m1.b
  h2, b2 = m2.a, m2.b
  h3, b3 = m3.a, m3.b
  h4, b4 = m4.a, m4.b
  h1, b1, h2, b2, h3, b3, h4, b4 = reassemble(h1, b1, h2, b2, h3, b3, h4, b4)
  if species == 0:
    case = 'All'
  elif species == 1 or species == 3:
    case = 'Cu'
  elif species == 2 or species == 4:
    case = 'Zr'
  else:
    raise ValueError('species value only 0,1 or 2!')
  if typ == 1:
    comm = 'As Prepared'
    tag = 'asp'
  elif typ == 2:
    comm = 'Annealed'
    tag = 'ann'
  else:
    raise ValueError('typ value only 1 or 2!')
  if cs == 0:
    csstr = 'All'
  elif cs == 1:
    csstr = 'Core'
  elif cs == 2:
    csstr = 'Interface'
  else:
    raise ValueError('cs value only 0, 1 or 2!')

  fig, ax = plt.subplots(figsize=(5, 6))
  ax.set_facecolor(cscheme(cs + 7))
  if cs == 1: ax.patch.set_alpha(0.2)
  if cs == 0:
    plt.grid(color='k',alpha=0.1, zorder=-2)
  else:
    plt.grid(color='w',zorder=-2)
  # plt.grid(color='w', zorder=-2)
  # plt.grid(color='m', alpha=0.1, zorder=-2)
  plt.setp(ax.spines.values(), linewidth=1.5)

  w = 0.15
  y1 = numpy.arange(len(b1)) + 2

  # ax.set_ylabel('.', color=(0, 0, 0, 0), labelpad=10)
  plt.ylabel('Atomic Fraction (at %)',fontsize='16')
  plt.xlabel('Polyhedron Index',fontsize='16')


  # fig.text(0.4, 0.96, r'VP for ' + name[:-4] + ' (' + comm + ') : ' + case, ha='center', fontsize=20)
  # if cs == 0:
  #   fig.text(0.75, 0.96, ' atoms', ha='left', fontsize=20)
  # elif cs == 1:
  #   fig.text(0.75, 0.96, ' ' + csstr + ' atoms', ha='left', color="red", fontsize=20)
  # elif cs == 2:
  #   fig.text(0.75, 0.96, ' ' + csstr + ' atoms', ha='left', color="blue", fontsize=20)
  # fig.text(0.5, 0.01, 'Polyhedron Index', ha='center', fontsize='15')
  # fig.text(0.01, 0.5, 'Atomic Fraction (at %)', va='center', rotation='vertical', fontsize='15')
  # fig.text(0.85,0.68,"{:.0e}".format(float(name.split(" ")[-1]))+r'$\frac{K}{s}$',weight='bold',fontsize='13')

  plt.sca(ax)

  p1 = ax.bar(y1, h1, width=0.15, label=l1, color=cscheme(1), edgecolor='black',zorder=3)
  p2 = ax.bar(y1 + w, h2, width=0.15, label=l2, color=cscheme(2), edgecolor='black',zorder=3)
  p3 = ax.bar(y1 + 2 * w, h3, width=0.15, label=l3, color=cscheme(3), edgecolor='black', hatch=".",zorder=3)
  p4 = ax.bar(y1 + 3 * w, h4, width=0.15, label=l4, color=cscheme(4), edgecolor='black', hatch="++",zorder=3)

  plt.xticks(y1, b1, rotation=90)

  qr = name.split(" ")[-1].split("e") #bbox_to_anchor=(1, 0.6),
  # if species == 0: ax.legend(fontsize=13, loc=1, title=r'${10}^{' + str(qr[1]) + '}$' + r'$K/s$')
  # else: ax.legend(fontsize=13, loc=1, title=r'${10}^{' + str(qr[1]) + '}$' + r'$K/s$' + case)
  if species == 0: ax.legend(fontsize=13, loc=1)
  else: ax.legend(fontsize=13, loc=1, title=case)
  # if species==0: ax.legend(fontsize=12,title="{:.0e}".format(float(name.split(" ")[-1]))+r'$\frac{K}{s}$')
  # else: ax.legend(fontsize=12,title="{:.0e}".format(float(name.split(" ")[-1]))+r'$\frac{K}{s}$, '+case+'')


  #    ax.invert_yaxis()
  yl = max(max(h1, h2, h3, h4)) * 1.05
  ax.set_ylim(0, yl)
  # ax.set_ylabel(' ', fontsize=12)
  ax.yaxis.set_label_position('left')
  ax.tick_params(axis="x", direction="in", left='off', labelsize=14)
  ax.tick_params(axis="y", direction="in", bottom="on", top='on', labelsize=14)
  plt.subplots_adjust(bottom=0.4, top=0.8)

  # fig.tight_layout(rect=[0.1, 0.03, 1, 0.95])  #
  fig.tight_layout()
  if cs == 0:
    plt.savefig('voronoi_' + name.replace(" ", "_") + '_' + tag + '_' + case + '.png',dpi=400)
  elif cs > 0:
    plt.savefig('voronoi_' + name.replace(" ", "_") + '_' + tag + '_' + case + '_' + csstr + '.png',dpi=400)

  plt.clf()
  plt.close()

def voro_plot_case(name, m1, m2, m3, case, cs, typ):
  h1, b1 = m1.a, m1.b
  h2, b2 = m2.a, m2.b
  h3, b3 = m3.a, m3.b
  #    h1,b1,h2,b2,h3,b3,h4,b4=reassemble(h1,b1,h2,b2,h3,b3,h3,b3)

  if typ == 1:
    comm = 'As Prepared'
    tag = 'asp'
  elif typ == 2:
    comm = 'Annealed'
    tag = 'ann'
  else:
    raise ValueError('typ value only 1 or 2!')
  if cs == 0:
    csstr = 'All'
  elif cs == 1:
    csstr = 'Core'
  elif cs == 2:
    csstr = 'Interface'
  else:
    raise ValueError('cs value only 0, 1 or 2!')

  fig, ax = plt.subplots(1, 3, sharex=True, figsize=(12, 12))
  y1 = numpy.arange(len(b1)) + 1
  y2 = numpy.arange(len(b2)) + 1
  y3 = numpy.arange(len(b3)) + 1

  ax[0].set_ylabel('.', color=(0, 0, 0, 0), labelpad=10)
  #    fig.suptitle('VP for '+str(siz)+'nm $Cu_{50}Zr_{50}$: '+case+' '+comm, fontsize=20) #, pad=20)
  fig.text(0.4, 0.96, r'VP for: ' + name[:-4] + ': ' +case + ' (' + comm+ ')', ha='center', fontsize=20)
  #    if cs==0: fig.text(0.6,0.9,' atoms', ha='left', fontsize=20)
  if cs == 1:
    fig.text(0.65, 0.96, ' ' + csstr + ' atoms', ha='left', color="red", fontsize=20)
  elif cs == 2:
    fig.text(0.65, 0.96, ' ' + csstr + ' atoms', ha='left', color="blue", fontsize=20)
  fig.text(0.01, 0.5, 'Polyhedron Index', va='center', rotation='vertical', fontsize='22')
  fig.text(0.5, 0.01, 'Atomic Fraction (at %)', ha='center', fontsize='20')
  fig.text(0.9,0.95,"{:.0e}".format(float(name.split(" ")[-1]))+r'$\frac{K}{s}$',weight='bold',fontsize='15')

  plt.sca(ax[0])
  p1 = ax[0].barh(y1, h1, height=0.25, label='All', color='black', edgecolor='black')
  plt.yticks(y1, b1)  # ,rotation=90)
  plt.legend(handles=[p1], loc='lower right', frameon=False, fontsize=20)
  ax[0].invert_yaxis()
  ax[0].set_xlim([0, 13])
  ax[0].set_xlabel('around all atoms', fontsize=18)
  ax[0].xaxis.set_label_position('top')
  ax[0].tick_params(axis="y", direction="in", left='off', labelsize=20)
  ax[0].tick_params(axis="x", direction="in", bottom="on", top='on', labelsize=20)
  for i, v in enumerate(h1):
    ax[0].text(v + 0.25, i + 1, str(v), color='blue', fontweight='bold')
  plt.subplots_adjust(bottom=0.4, top=0.8)

  plt.sca(ax[1])
  p2 = ax[1].barh(y2, h2, height=0.25, label='Cu', color='black', edgecolor='black')
  plt.yticks(y2, b2)  # ,rotation=90)
  plt.legend(handles=[p2], loc='lower right', frameon=False, fontsize=20)
  ax[1].invert_yaxis()
  ax[1].set_xlim([0, 13])
  ax[1].set_xlabel('around Cu atoms', fontsize=18)
  ax[1].xaxis.set_label_position('top')
  ax[1].tick_params(axis="y", direction="in", left='off', labelsize=20)
  ax[1].tick_params(axis="x", direction="in", bottom="on", top='on', labelsize=20)
  for i, v in enumerate(h2):
    ax[1].text(v + 0.25, i + 1, str(v), color='blue', fontweight='bold')
  plt.subplots_adjust(bottom=0.4)  # , top=0.8)

  plt.sca(ax[2])
  p3 = ax[2].barh(y3, h3, height=0.25, label='Zr', color='black', edgecolor='black')
  plt.yticks(y3, b3)  # ,rotation=90)
  plt.legend(handles=[p3], loc='lower right', frameon=False, fontsize=20)
  ax[2].invert_yaxis()
  ax[2].set_xlim([0, 13])
  ax[2].set_xlabel('around Zr atoms', fontsize=18)
  ax[2].xaxis.set_label_position('top')
  ax2 = ax[2].twiny()
  ax2.set_xlabel("hi", fontsize=40, color='white', labelpad=4)
  ax2.xaxis.set_label_position('bottom')
  ax2.tick_params(axis="x", top='off', labelsize=0)
  ax[2].tick_params(axis="y", direction="in", left='off', labelsize=20)
  ax[2].tick_params(axis="x", direction="in", bottom="on", top='on', labelsize=20)
  for i, v in enumerate(h3):
    ax[2].text(v + 0.25, i + 1, str(v), color='blue', fontweight='bold')
  plt.subplots_adjust(bottom=0.4)  # , top=0.8)

  fig.tight_layout(rect=[0, 0.03, 1, 0.95])  #
  if cs == 0:
    plt.savefig('voronoi_' + name.replace(" ", "_") + tag + '_' + case + '.png')
  elif cs > 0:
    plt.savefig('voronoi_' + name.replace(" ", "_") + tag + '_' + case + '_' + csstr + '.png')

  plt.clf()
  plt.close()

  rows = zip(b1,h1,b2,h2,b3,h3)
#write out tmp files with voronoi data for each case: All, Core, Interface in six columns in total
  fil = 'tmp.voronoi_' + name.replace(" ", "_") + '_' + tag + '_' + case
  with open(fil, 'w') as f:
    writer = csv.writer(f)
    writer.writerow('Index,All,Index,Cu,Index,Zr')
    for row in rows:
      writer.writerow(row)

def plot_atmhisto(name,m1, m2, m3, m4, l1, l2, l3, l4, species, cs, typ):
  plot_atmhisto_nico(name,m1, m2, m3, m4, l1, l2, l3, l4, species, cs, typ)
  h1, n1 = m1.c, m1.d
  h2, n2 = m2.c, m2.d
  h3, n3 = m3.c, m3.d
  h4, n4 = m4.c, m4.d

  if species == 0:
    case = 'All'
  elif species == 1 or species == 3:
    case = 'Cu'
  elif species == 2 or species == 4:
    case = 'Zr'
  else:
    raise ValueError('species value only 0,1 or 2!')
  if cs == 0:
    csstr = 'All'
  elif cs == 1:
    csstr = 'Core'
  elif cs == 2:
    csstr = 'Interface'
  else:
    raise ValueError('cs value only 0, 1 or 2!')
  if typ == 1:
    comm = 'As Prepared'
    tag = 'asp'
  elif typ == 2:
    comm = 'Annealed'
    tag = 'ann'
  else:
    raise ValueError('typ value only 1 or 2!')

  bin1 = math.floor(n1 / 50)
  bin2 = math.floor(n2 / 50)
  bin3 = math.floor(n3 / 50)
  bin4 = math.floor(n4 / 50)

  #    y1,x1=numpy.histogram(amg.flatten(),bins=bin1,density="True")
  #    c1 = (x1[:-1] + x1[1:]) / 2
  fig, ax = plt.subplots(figsize=(10, 5), sharey=True)

  ax.set_facecolor(cscheme(cs + 7))
  if cs == 1: ax.patch.set_alpha(0.2)
  if cs == 0:
    plt.grid(color='k',alpha=0.1, zorder=-2)
  else:
    plt.grid(color='w',zorder=-2)
  # plt.grid(color='w', zorder=-2)

  # plt.grid(color='m', alpha=0.1, zorder=-2)
  plt.setp(ax.spines.values(), linewidth=1.5)
  plt.sca(ax)

  bin1 = bin2 = bin3 = bin4 = bin5 = bin6 = 100

  ht1, b1 = numpy.histogram(h1, bins=bin1)
  ht2, b2 = numpy.histogram(h2, bins=bin1)
  ht3, b3 = numpy.histogram(h3, bins=bin1)
  ht4, b4 = numpy.histogram(h4, bins=bin1)
  a1, c1 = ht1 / n1, (b1[:-1] + b1[1:]) / 2
  a2, c2 = ht2 / n2, (b2[:-1] + b2[1:]) / 2
  a3, c3 = ht3 / n3, (b3[:-1] + b3[1:]) / 2
  a4, c4 = ht4 / n4, (b4[:-1] + b4[1:]) / 2

  # print('His tot', ht1.sum(axis=0) / n1, ht2.sum(axis=0) / n2, ht3.sum(axis=0) / n3, ht4.sum(axis=0) / n4,
  #       ht5.sum(axis=0) / n5, ht6.sum(axis=0) / n6)

  plt.plot(c1, a1, color=cscheme(1), label=l1,zorder=3,lw=2)
  plt.plot(c2, a2, color=cscheme(2), label=l2,zorder=3,lw=2)
  plt.plot(c3, a3, color=cscheme(3), label=l3,zorder=3,lw=2)
  plt.plot(c4, a4, color=cscheme(4), label=l4,zorder=3,lw=2)

  # qr = name.split(" ")[-1].split("e")
  # ax.legend(fontsize=12, bbox_to_anchor=(1, 0.6), loc='upper right', title= r'${10}^{' + str(qr[1]) + '}$' + r'$K/s$')
  ax.legend(fontsize=12, bbox_to_anchor=(1, 0.6), loc='upper right')

  ax.set_xlabel('Volume per Atom (Units)', fontsize=18)
  ax.set_ylabel('Normalised Counts ', fontsize=18)
  ax.yaxis.set_label_position('left')
  ax.tick_params(axis="both", direction="in", bottom=True, top=True, left=True, right=True, labelsize=15)
  plt.subplots_adjust(bottom=0.4, top=0.8)

  y1=max(max(a1),max(a2),max(a3),max(a4))
  y2=min(min(a1),min(a2),min(a3),min(a4))
  x1=max(max(c1),max(c2),max(c3),max(c4))
  x2=min(min(c1),min(c2),min(c3),min(c4))
  #print('Max Y',yl)
  yl1=y1*1.05
  yl2=y2*0.90
  xl1 = x1 * 1.05
  xl2 = x2 * 0.90
  ax.set_ylim(yl2, yl1)
  ax.set_xlim(xl2, xl1)
  # ax.set_yscale("log")

  plt.text(15,0.8*yl1,r'Cu',weight='bold',fontsize='13')
  plt.text(23,0.6*yl1,r'Zr',weight='bold',fontsize='13')

  rng = [l1, l2, l3, l4]

  # v1 = ht1 * c1 / n1
  # v2 = ht2 * c2 / n2
  # v3 = ht3 * c3 / n3
  # v4 = ht4 * c4 / n4
  # v5 = ht5 * c5 / n5
  # v6 = ht6 * c6 / n6
  v1 = h1 / n1
  v2 = h2 / n2
  v3 = h3 / n3
  v4 = h4 / n4

  tot_v = [v1.sum(axis=0), v2.sum(axis=0), v3.sum(axis=0), v4.sum(axis=0)]

  fig.tight_layout(rect=[0, 0.03, 1, 0.95])
  if cs == 0 and species==0:
    plt.savefig('vol-atom_' +  name.replace(" ", "_") + '_' + tag + '_' + case + '.png',dpi=400)
  elif cs > 0 and species==0:
    plt.savefig('vol-atom_' +  name.replace(" ", "_") + '_' + tag + '_' + case + '_' + csstr + '.png',dpi=400)
  plt.clf()
  plt.close()

  if species == 0:
    rows=zip(tot_v,rng)
    # fil = 'tmp.vol_coll_' + name.replace(" ", "_") + '_' + tag + '_' + case + '_' + csstr
    fil = 'tmp.vol_coll_' + name.replace(" ", "_") + '_' + tag + '_' + csstr
    with open(fil, 'w') as f:
      writer = csv.writer(f)
      for row in rows:
        writer.writerow(row)

    rows2=zip(c1,a1,c2,a2,c3,a3,c4,a4)
    fil = 'tmp.atvolume_' + name.replace(" ", "_") + '_' + tag + '_' + csstr
    with open(fil, 'w') as f:
      writer = csv.writer(f)
      writer.writerow([l1,'Counts',l2,'Counts',l3,'Counts',l4,'Counts'])
      for row in rows2:
        writer.writerow(row)

def plot_chains(name, m1, m2, m3, m4, l1, l2, l3, l4, species, cs, typ):
  if species == 0:
    case = 'All'
  elif species == 1 or species == 3:
    case = 'Cu'
  elif species == 2 or species == 4:
    case = 'Zr'
  else:
    raise ValueError('species value only 0,1 or 2!')
  if typ == 1:
    comm = 'As Prepared'
    tag = 'asp'
  elif typ == 2:
    comm = 'Annealed'
    tag = 'ann'
  else:
    raise ValueError('typ value only 1 or 2!')
  if cs == 0:
    csstr = 'All'
  elif cs == 1:
    csstr = 'Core'
  elif cs == 2:
    csstr = 'Interface'
  else:
    raise ValueError('cs value only 0, 1 or 2!')

  # print(max(m1.csa),max(m2.csa),max(m3.csa),max(m4.csa),max(m5.csa),max(m6.csa))
  ht1, b1 = numpy.histogram(m1.csa, max(m1.csa)) #change c1 to m1's property which holds cluster size array
  ht2, b2 = numpy.histogram(m2.csa, max(m2.csa))
  ht3, b3 = numpy.histogram(m3.csa, max(m3.csa))
  ht4, b4 = numpy.histogram(m4.csa, max(m4.csa))

  # y1, c1 = ht1 / m1.d, (b1[:-1] + b1[1:]) / 2
  # y2, c2 = ht2 / m2.d, (b2[:-1] + b2[1:]) / 2
  # y3, c3 = ht3 / m3.d, (b3[:-1] + b3[1:]) / 2
  # y4, c4 = ht4 / m4.d, (b4[:-1] + b4[1:]) / 2
  # y5, c5 = ht5 / m5.d, (b5[:-1] + b5[1:]) / 2
  # y6, c6 = ht6 / m6.d, (b6[:-1] + b6[1:]) / 2
  # y1, c1 = ht1 / m1.d, [x for x in range(1,len(b1))] #works only in this situation where bins are integers
  # y2, c2 = ht2 / m2.d, [x for x in range(1,len(b2))]
  # y3, c3 = ht3 / m3.d, [x for x in range(1,len(b3))]
  # y4, c4 = ht4 / m4.d, [x for x in range(1,len(b4))]
  # y5, c5 = ht5 / m5.d, [x for x in range(1,len(b5))]
  # y6, c6 = ht6 / m6.d, [x for x in range(1,len(b6))]
  y1, c1 = ht1 / sum(m1.csa), [x for x in range(1,len(b1))] #works only in this situation where bins are integers
  y2, c2 = ht2 / sum(m2.csa), [x for x in range(1,len(b2))]
  y3, c3 = ht3 / sum(m3.csa), [x for x in range(1,len(b3))]
  y4, c4 = ht4 / sum(m4.csa), [x for x in range(1,len(b4))]

  # print('SUM M1 CSA',sum(m1.csa),'TOT M1 atoms',m1.d)
  # print('Max csa m1',max(m1.csa),'Num bins', len(c1))
  fig, ax = plt.subplots(figsize=(5, 5), sharey=True)
  ax.set_facecolor(cscheme(cs + 7))
  if cs == 1: ax.patch.set_alpha(0.2)
  if cs == 0:
    plt.grid(color='k',alpha=0.1, zorder=-2)
  else:
    plt.grid(color='w',zorder=-2)
  # plt.grid(color='w', zorder=-2)
  # plt.grid(color='m', alpha=0.1, zorder=-2)
  plt.setp(ax.spines.values(), linewidth=1.5)

  # for i in [3,4,5,6,7,8,9,10]:
  # maxch = max(max(m1.csa),max(m2.csa),max(m3.csa),max(m4.csa),max(m5.csa),max(m6.csa))
  maxch = 100
  iter = math.floor(maxch/10)
  lst = [x for x in range(3,10,1)]+[x for x in range(20,40,10)]
  markers = ['o', 'v', '^', '<', '>', 's', 'p', 'P',  'h', 'H', '+', 'x', 'X', '*', '1', '2', '3', '4', '8',
             'D', 'd', '|', '_']
  it = 0
  for i in lst:
    it+=1
  # for i in range(0,maxch,iter):
    # chain_tot = [sum([y for y in x if y>i])/sum(x) for x in [m1.csa,m2.csa,m3.csa,m4.csa,m5.csa,m6.csa]]
    # chain_tot = [len([y for y in x.csa if y >= i])*100/len(x.csa) for x in [m1, m2, m3, m4, m5, m6]] #before it was x.d
    chain_tot = [len([y for y in x.csa if y <= i])*100/len(x.csa) for x in [m1, m2, m3, m4]] #before it was x.d

    # chain_tot = [numpy.sum(x[i:]) for x in [y1,y2,y3,y4,y5,y6]]
    plt.plot([1,2,3,4],chain_tot,'k-',marker=markers[it],label=str(i),lw=1, zorder=3)
    # if cs == 0: plt.plot([1,2,3,4,5,6],chain_tot,'k-',marker=markers[it],label=str(i),lw=1, zorder=3)
    # else: plt.plot([1,2,3,4],chain_tot[2:],'k-',marker=markers[it],label=str(i),lw=1, zorder=3)
  # plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
  plt.ylabel(r'$N_{Strings}$ (%)',fontsize=16)
  # ax.legend(title="{:.0e}".format(float(name.split(" ")[-1]))+r'$\frac{K}{s}$'
                  # +'\n Minimum \n Chain Size',loc=0,fontsize=12) #bbox_to_anchor=(.95, 1))
  # ax.legend(title='Minimum \nChain Size',loc=0,fontsize=10,bbox_to_anchor=(1.025, 1.02))
  # qr = name.split(" ")[-1].split("e")  # bbox_to_anchor=(1, 0.6),
  # qrvl=r'${10}^{' + str(qr[1]) + '}$' + r'$K/s$'
  ax.legend(title='\nMaximum \nString Size',loc=0,fontsize=10,bbox_to_anchor=(1.035, 1.02))

  # if cs == 0: plt.xticks([1,2,3,4,5,6],[l1,l2,l3,l4,l5,l6], rotation=45) #labelsize=15)
  # else: plt.xticks([1,2,3,4],[l3,l4,l5,l6], rotation=45) #labelsize=15)
  plt.xticks([1, 2, 3, 4], [l1, l2, l3, l4], rotation=45)  # labelsize=15)

  ax.yaxis.set_label_position('left')
  ax.tick_params(axis="x", direction="in", left='on', labelsize=14)
  ax.tick_params(axis="y", direction="in", bottom="on", top='on', labelsize=14)
  ax.yaxis.set_label_position('left')
  ax.tick_params(axis="x", direction="in", left='on', labelsize=12)
  ax.tick_params(axis="y", direction="in", bottom="on", top='on', labelsize=12)
  fig.tight_layout(rect=[0, 0, 0.98, 0.95])

  if cs == 0:
    plt.savefig('chains_ico_' + name.replace(" ", "_") + '_' + tag + '_' + case + '.png',dpi=400)
  elif cs > 0:
    plt.savefig('chains_ico_' + name.replace(" ", "_") +'_' + tag + '_' + case + '_' + csstr + '.png',dpi=400)
  plt.clf()
  plt.close()

  fig, ax = plt.subplots(1, 1, figsize=(5,5), sharey=True)
  ax.set_facecolor(cscheme(cs + 7))
  if cs == 1: ax.patch.set_alpha(0.2)
  if cs == 0:
    plt.grid(color='k',alpha=0.1, zorder=-2)
  else:
    plt.grid(color='w',zorder=-2)
  # plt.grid(color='w', zorder=-2)
  plt.setp(ax.spines.values(), linewidth=1.5)
  # plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

  # s1 = sum([x * y for x, y in zip(c1, ht1)]) / sum(m1.csa) #these will all add up to one
  # s2 = sum([x * y for x, y in zip(c2, ht2)]) / sum(m2.csa)
  # s3 = sum([x * y for x, y in zip(c3, ht3)]) / sum(m3.csa)
  # s4 = sum([x * y for x, y in zip(c4, ht4)]) /sum(m4.csa)
  # s5 = sum([x * y for x, y in zip(c5, ht5)]) /sum(m5.csa)
  # s6 = sum([x * y for x, y in zip(c6, ht6)]) /sum(m6.csa)

  # print('heights',ht3)
  # # print('bins',c3)
  # print('bins',b1,'len',len(b1))

  # for i in [3,4,5,6,7,8,9,10]:
  # # for i in [10]:
  #   s3 = sum([x * y for x, y in zip(c3, y3) if x <= i])
  #   s4 = sum([x * y for x, y in zip(c4, y4) if x <= i])
  #   s5 = sum([x * y for x, y in zip(c5, y5) if x <= i])
  #   s6 = sum([x * y for x, y in zip(c6, y6) if x <= i])
  #   plt.plot([1,2,3,4],[s3,s4,s5,s6],'o-',label=str(i),lw=2, zorder=3)
  # # plt.plot([1,2,3,4],[s3,s4,s5,s6],lw=2,zorder=3)
  # #  ax.legend(title="{:.0e}".format(float(name.split(" ")[-1])) + r'$\frac{K}{s}$', loc=0, fontsize=12)  # bbox_to_anchor=(.95, 1))
  # ax.legend(title='Minimum \nChain Size', loc=0, fontsize=10, bbox_to_anchor=(1.025, 1.02))
  # # print(s3,s4,s5,s6)
  # plt.xticks([1,2,3,4], [l3, l4, l5, l6], rotation=45,fontsize=12)
  # plt.yticks(fontsize=12)
  # plt.ylabel('Atomic Fraction',fontsize=15)
  # ax.yaxis.set_label_position('left')
  # ax.tick_params(axis="x", direction="in", left='on', labelsize=12)
  # ax.tick_params(axis="y", direction="in", bottom="on", top='on', labelsize=12)
  # fig.tight_layout(rect=[0, 0.03, 0.98, 0.95])


  # axins1 = inset_axes(ax,"40%",  # width = 30% of parent_bbox
  #                     height=1.2,  # height : 1 inch
  #                     loc=9, bbox_to_anchor=(-0.05,0,1,1),
  #                     bbox_transform=ax.transAxes)
  #                     bbox_to_anchor = (0.9,0.9,1,1))

  p1 = [x * y *100 for x, y in zip(c1, y1)]
  p2 = [x * y *100 for x, y in zip(c2, y2)]
  p3 = [x * y *100 for x, y in zip(c3, y3)]
  p4 = [x * y *100 for x, y in zip(c4, y4)]
  # plt.plot([1,2,3,4],[s3,s4,s5,s6],'o-',label=str(i),lw=2, zorder=3)

  # print('sum')
  # print(s3)
  itr = 1
  for i,j,k in zip([c1,c2,c3,c4],[p1,p2,p3,p4],[l1,l2,l3,l4]):
  # for i,j,k in zip([c3,c4,c5,c6],[p3,p4,p5,p6],[l3,l4,l5,l6]):
    itr += 1
    ax.plot(i,j, 'o-', label=str(k),marker=markers[itr-2], color=cscheme(itr-1), lw=2, zorder=itr+1)
    # ax.legend(title="{:.0e}".format(float(name.split(" ")[-1])) + r'$\frac{K}{s}$', loc=0, fontsize=10) #, bbox_to_anchor=(0.95, 1.02))
    qr = name.split(" ")[-1].split("e")  # bbox_to_anchor=(1, 0.6),
    # ax.legend(fontsize=10, loc=0, title=r'${10}^{' + str(qr[1]) + '}$' + r'$K/s$')
    ax.legend(loc=0, fontsize=10)  # , bbox_to_anchor=(0.95, 1.02))
  # print(s3,s4,s5,s6)
  # plt.xticks([1,2,3,4], [l3, l4, l5, l6], rotation=45,fontsize=12)
  # ax.set_yticks(fontsize=12)
  ax.set_ylabel(r'$N_{atoms}$ (%)',fontsize=16)
  ax.set_xlabel('Strings Size',fontsize=16)
  # ax.set_yscale("log")
  fig.tight_layout(rect=[0, 0.03, 0.95, 0.95])

  # if cs>0: ax.set_xlim(0,150)
  # else: ax.set_xlim(0,1500)
  # ax.set_ylim(0,0.03)
  ax.yaxis.set_label_position('left')
  ax.tick_params(axis="x", direction="in", left='on', labelsize=12)
  ax.tick_params(axis="y", direction="in", bottom="on", top='on', labelsize=12)

  if cs == 0:
    plt.savefig('wtdsum_ico_' + name.replace(" ", "_") + '_' + tag + '_' + case + '.png',dpi=400)
  elif cs > 0:
    plt.savefig('wtdsum_ico_' + name.replace(" ", "_") +'_' + tag + '_' + case + '_' + csstr + '.png',dpi=400)

  plt.clf()
  plt.close()

  fig, ax = plt.subplots(1, 1, figsize=(5, 5), sharey=True)
  ax.set_facecolor(cscheme(cs + 7))
  if cs == 1: ax.patch.set_alpha(0.2)
  if cs == 0:
    plt.grid(color='k', alpha=0.1, zorder=-2)
  else:
    plt.grid(color='w', zorder=-2)
  # plt.grid(color='w', zorder=-2)
  plt.setp(ax.spines.values(), linewidth=1.5)

  avgclus = [sum(x) / len(x) for x in [m1.csa, m2.csa, m3.csa, m4.csa]]
  ax.plot([1, 2, 3, 4], avgclus, 'o:', c='black', mfc='blue',mec='k', markersize=9)
  ax.set_ylabel('Avg String Size', fontsize=16)
  ax.set_xticks([1, 2, 3, 4])
  ax.set_xticklabels([l1, l2, l3, l4], rotation=45, fontsize=14)
  ax.tick_params(axis="y", direction="in", bottom="on", top='on', labelsize=14)

  fig.tight_layout(rect=[0, 0.03, 0.95, 0.95])

  if cs == 0:
    plt.savefig('avg_chains_' + name.replace(" ", "_") + '_' + tag + '_' + case + '.png',dpi=400)
  elif cs > 0:
    plt.savefig('avg_chains_' + name.replace(" ", "_") +'_' + tag + '_' + case + '_' + csstr + '.png',dpi=400)

  plt.clf()
  plt.close()