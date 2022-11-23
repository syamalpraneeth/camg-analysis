#!/usr/bin/env python3

import os, sys
from ovito.io import *
from ovito.modifiers import *
import numpy
import matplotlib.pyplot as plt
import csv
sys.path.append('..')
from analysis import *
from plotting import *
from energy_minim2 import *

src='/home/mj0054/Documents/work/simulations/projects/'
name='Cu50Zr50 random NG seg 5GPa' #syntax: 3nm Cu50Zr50 CAMG vs MG

# f1,l1=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_1nm','1nm'
# f2,l2=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_1.5nm','1.5nm'
# f3,l3=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_2nm','2nm'
# f4,l4=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_2.25nm','2.25nm'
# f5,l5=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_2.5nm','2.5nm'
# f6,l6=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_3nm','3nm'
# f7,l7=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_5nm','5nm'
# f8,l8=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_7nm','7nm'

f1,l1=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_1nm','1nm'
f2,l2=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_1.5nm','1.5nm'
f3,l3=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_1.8nm','1.8nm'
f4,l4=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_2nm','2nm'
f5,l5=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_2.25nm','2.25nm'
f6,l6=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_2.49nm','2.49nm'
f7,l7=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_2.5nm','2.5nm'
f8,l8=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_2.75nm','2.75nm'
f9,l9=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_2.8nm','2.8nm'
f10,l10=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_3nm','3nm'
f11,l11=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_3.1nm','3.1nm'
f12,l12=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_5nm','5nm'
f13,l13=src+'nanoglass/50-50/1e10/segregated/random/5GPa/data/dump.ng_unl_7nm','7nm'

def cleardata(source,typ):
  if typ == 1: tag = 'asp'
  elif typ == 2: tag = 'ann'
  for f in glob.glob('tmp.'+source+'_coll_' + name.replace(" ", "_") + '_' + tag + '*'):
    with open("deleted_"+source+"_files.txt", "a") as text_file:
      text_file.write(f+'\n')
    if os.path.exists(f): os.remove(f)

"""
#pe/atom function
# def pe_atom(f1,f2,f3,f4,f5,f6,f7,f8,l1,l2,l3,l4,l5,l6,l7,l8,typ):
def pe_atom(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, typ):
  cleardata('pote',typ)
  for cs in range(0,1):
    # if cs == 0: minimize(name,src,f1,l1,f2,l2,f3,l3,f4,l4,f5,l5,f6,l6)
    ovito_pote(name,f1,l1,cs,0,typ)
    ovito_pote(name,f2,l2,cs,0,typ)
    ovito_pote(name,f3,l3,cs,1,typ)
    ovito_pote(name,f4,l4,cs,1,typ)
    ovito_pote(name,f5,l5,cs,1,typ)
    ovito_pote(name,f6,l6,cs,1,typ)
    ovito_pote(name,f7,l7,cs,1,typ)
    ovito_pote(name,f8,l8,cs,1,typ)
    ovito_pote(name,f9,l9,cs,1,typ)
    ovito_pote(name,f10,l10,cs,1,typ)
    ovito_pote(name,f11,l11,cs,1,typ)
    ovito_pote(name,f12,l12,cs,1,typ)
    ovito_pote(name,f13,l13,cs,1,typ)

    # pe_plot_species(name,l1,l2,l3,l4,l5,l6,l7,l8,cs,typ)

#rdf function
# def rdf(f1,f2,f3,f4,f5,f6,f7,f8,l1,l2,l3,l4,l5,l6,l7,l8,typ):
def rdf(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, typ):

  cleardata('porosity',typ)
  cleardata('density',typ)
  for cs in range(0,1):
    ovito_rdf(name,f1,l1,cs,0,typ)
    ovito_rdf(name,f2,l2,cs,0,typ)
    ovito_rdf(name,f3,l3,cs,1,typ)
    ovito_rdf(name,f4,l4,cs,1,typ)
    ovito_rdf(name,f5,l5,cs,1,typ)
    ovito_rdf(name,f6,l6,cs,1,typ)
    ovito_rdf(name,f7,l7,cs,1,typ)
    ovito_rdf(name,f8,l8,cs,1,typ)
    ovito_rdf(name,f9,l9,cs,1,typ)
    ovito_rdf(name,f10,l10,cs,1,typ)
    ovito_rdf(name,f11,l11,cs,1,typ)
    ovito_rdf(name,f12,l12,cs,1,typ)
    ovito_rdf(name,f13,l13,cs,1,typ)

    plot_rdf(name,l1,0,cs,typ)
    plot_rdf(name,l2,0,cs,typ)
    plot_rdf(name,l3,1,cs,typ)
    plot_rdf(name,l4,1,cs,typ)
    plot_rdf(name,l5,1,cs,typ)
    plot_rdf(name,l6,1,cs,typ)
    plot_rdf(name,l7,1,cs,typ)
    plot_rdf(name,l8,1,cs,typ)
    plot_rdf(name,l9,1,cs,typ)
    plot_rdf(name,l10,1,cs,typ)
    plot_rdf(name,l11,1,cs,typ)
    plot_rdf(name,l12,1,cs,typ)
    plot_rdf(name,l13,1,cs,typ)

    # plot_allrdf(name,l1,l2,l3,l4,l5,l6,l7,l8,cs,typ)
  # os._exit(1)
"""

#voronoi
# def voronoi(f1,f2,f3,f4,f5,f6,f7,f8,l1,l2,l3,l4,l5,l6,l7,l8,typ):
def voronoi(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, typ):
  cleardata('strain',typ)
  for cs in range(0,1):
    d1a, d1b, d1c = voro_data_film(name, f1, l1, cs, typ)
    #voro_plot_case(name,d1a,d1b,d1c,l1,0,1)
    # d2a,d2b,d2c = d1a,d1b,d1c
    d2a, d2b, d2c = voro_data_film(name, f2, l2, cs, typ)
    # voro_plot_case(name, d2a, d2b, d2c, l2, 0, 1)
    # d3a,d3b,d3c = d1a,d1b,d1c
    # d4a,d4b,d4c = d1a,d1b,d1c
    # d5a,d5b,d5c = d1a,d1b,d1c
    d3a,d3b,d3c = voro_data_film(name, f3,l3,cs, typ)
    d4a,d4b,d4c = voro_data_film(name, f4,l4,cs, typ)
    d5a,d5b,d5c = voro_data_film(name, f5,l5,cs, typ)
    d6a,d6b,d6c = voro_data_film(name, f6,l6,cs, typ)
    d7a,d7b,d7c = voro_data_film(name, f7,l7,cs, typ)
    d8a,d8b,d8c = voro_data_film(name, f8,l8,cs, typ)
    d9a,d9b,d79c = voro_data_film(name, f9,l9,cs, typ)
    d10a,d10b,d10c = voro_data_film(name, f10,l10,cs, typ)
    d11a,d11b,d11c = voro_data_film(name, f11,l11,cs, typ)
    d12a,d12b,d12c = voro_data_film(name, f12,l12,cs, typ)
    d13a,d13b,d13c = voro_data_film(name, f13,l13,cs, typ)

    #voro_plot_species(name, m1, m2, m3, m4, m5, m6, l1, l2, l3, l4, l5, l6, species, cs, typ):
    #voro_plot_species(name,d1a,d2a,d3a,d4a,d5a,l1,l2,l3,l4,l5,0,cs,typ) #All
    # voro_plot_species(name,d1b,d2b,d3b,d4b,d5b,d6b,l1,l2,l3,l4,l5,l6,1,cs,1) #Cu
    # voro_plot_species(name,d1c,d2c,d3c,d4c,d5c,d6c,l1,l2,l3,l4,l5,l6,2,cs,1) #Zr
    # voro_plot_case(name,d3a,d3b,d3c,l3,cs,1)
    # voro_plot_case(name,d4a,d4b,d4c,l4,cs,1)
    # voro_plot_case(name,d5a,d5b,d5c,l5,cs,1)
    # voro_plot_case(name,d6a,d6b,d6c,l6,cs,1)
    # collate_ico(name,d1a,d2a,d3a,d4a,d5a,d6a,d7a,d8a,d1b,d2b,d3b,d4b,d5b,d6b,d7b,d8b,l1,l2,l3,l4,l5,l6,l7,l8,cs,typ) #All and Cu only
    collate_ico(name,d2a,d3a,d4a,d5a,d6a,d7a,d8a,d9a,d2b,d3b,d4b,d5b,d6b,d7b,d8b,d9b,l2,l3,l4,l5,l6,l7,l8,l9,cs,typ) #All and Cu only
    # plot_chains(name,d1a,d2a,d3a,d4a,d5a,d6a,l1,l2,l3,l4,l5,l6,0,cs,1) #All
    # plot_atmhisto(name,d1a,d2a,d3a,d4a,d5a,d6a,l1,l2,l3,l4,l5,l6,0,cs,1) #All
    # plot_atmhisto(name,d1b,d2b,d3b,d4b,d5b,d6b,l1,l2,l3,l4,l5,l6,3,cs,1) #Cu
    # plot_atmhisto(name,d1c,d2c,d3c,d4c,d5c,d6c,l1,l2,l3,l4,l5,l6,4,cs,1) #Zr

for typ in range(1,2):
  if typ==2:
    flist=[]
    for f in [f1,f2,f3,f4,f5,f6,f7,f8]:
      ind = f.find('/data')
      f = f[:ind] + '/anneal' + f[ind:]
      f= f.replace(f.split("/")[-1],'') + 'dump.ng_annealed_' + f.split("/")[-1].split("_")[-1]
      flist.append(f)
    [f1,f2,f3,f4,f5,f6,f7,f8] = flist

  # pe_atom(f1,f2,f3,f4,f5,f6,f7,f8,l1,l2,l3,l4,l5,l6,l7,l8,typ)
  # rdf(f1,f2,f3,f4,f5,f6,f7,f8,l1,l2,l3,l4,l5,l6,l7,l8,typ)
  # voronoi(f1,f2,f3,f4,f5,f6,f7,f8,l1,l2,l3,l4,l5,l6,l7,l8,typ)
  # pe_atom(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, typ)
  # rdf(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, typ)
  voronoi(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, typ)

os._exit(1)
f1,l1=src+'nanoglass/50-50/1e10/segregated/random/5GPa/dsc/continue/','2nm'
f2,l2=src+'nanoglass/50-50/1e10/segregated/random/5GPa/dsc/continue/','3nm'
f3,l3=src+'nanoglass/50-50/1e10/segregated/random/5GPa/dsc/continue/','5nm'
f4,l4=src+'nanoglass/50-50/1e10/segregated/random/5GPa/dsc/continue/','7nm'
cp_fit(f1,l1)
cp_fit(f2,l2)
cp_fit(f3,l3)
cp_fit(f4,l4)
