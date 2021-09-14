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
name='3nm Cu50Zr50 CAMG vs MG 1e14' #syntax: 3nm Cu50Zr50 CAMG vs MG
f1,l1=src+'bulk/msr/50-50/1e14/data/dump.msr_19','MG'
f2,l2=src+'bulk/msr/50-50/1e14/anneal/data/dump.msr_19','MG_ht'
#f3,l3=src+'nanoglass/50-50/1e14/random/5GPa/unload/data/dump.ng_unl_3nm','NG'
f3,l3=src+'/nanoglass/50-50/1e14/3nm/hcp/5GPa/restart/data/dump.ng_unl_3nm','NG'
f4,l4=src+'cibd/cibd_multi/50-50/50K/1e14/data/dump.dep.cibd_3nm_60meV-lay3','60meV'
f5,l5=src+'cibd/cibd_multi/50-50/50K/1e14/data/dump.dep.cibd_3nm_300meV-lay3','300meV'
f6,l6=src+'cibd/cibd_multi/50-50/50K/1e14/data/dump.dep.cibd_3nm_600meV-lay3','600meV'
# f5,l5=src+'cibd/cibd_multi/50-50/50K/penghui/2020_10_30/restart_300-600/data/dump.dep.cibd_3nm_600meV-lay4','600meV'
# f6,l6=src+'cibd/cibd_multi/50-50/50K/penghui/2020_12_04/data/dump.dep.cibd_3nm_6000meV-lay4','6000meV'

"""
#pe/atom AS_PREPARED
for f in glob.glob('tmp.pote_coll_' + name.replace(" ", "_") + '_asp_*'):
  if os.path.exists(f): os.remove(f)
for cs in range(0,3):
  # if cs == 0: minimize(name,src,f1,l1,f2,l2,f3,l3,f4,l4,f5,l5,f6,l6)
  ovito_pote(name,f1,l1,cs,0,1)
  ovito_pote(name,f2,l2,cs,0,1)
  ovito_pote(name,f3,l3,cs,1,1)
  ovito_pote(name,f4,l4,cs,1,1)
  ovito_pote(name,f5,l5,cs,1,1)
  ovito_pote(name,f6,l6,cs,1,1)
  pe_plot_species(name,l1,l2,l3,l4,l5,l6,cs,1)

#os._exit(1)

#rdf AS_PREPARED
for cs in range(0,3):
  ovito_rdf(f1,l1,cs,0)
  ovito_rdf(f2,l2,cs,0)
  ovito_rdf(f3,l3,cs,1)
  ovito_rdf(f4,l4,cs,1)
  ovito_rdf(f5,l5,cs,1)
  ovito_rdf(f6,l6,cs,1)

  plot_rdf(name,l1,0,cs,1)
  plot_rdf(name,l2,0,cs,1)
  plot_rdf(name,l3,1,cs,1)
  plot_rdf(name,l4,1,cs,1)
  plot_rdf(name,l5,1,cs,1)
  plot_rdf(name,l6,1,cs,1)

  plot_allrdf(name,l1,l2,l3,l4,l5,l6,cs,1)

# os._exit(1)
"""

#voronoi AS_PREPARED
d1a,d1b,d1c = voro_data_bulk(f1,l1,0)
#voro_plot_case(name,d1a,d1b,d1c,l1,0,1)
# d2a,d2b,d2c = d1a,d1b,d1c
d2a, d2b, d2c = voro_data_bulk(f2, l2, 0)
# voro_plot_case(name, d2a, d2b, d2c, l2, 0, 1)

for cs in range(0,3):
  # d3a,d3b,d3c = d1a,d1b,d1c
  # d4a,d4b,d4c = d1a,d1b,d1c
  # d5a,d5b,d5c = d1a,d1b,d1c
  # d6a,d6b,d6c = d1a, d1b, d1c
  d3a,d3b,d3c = voro_data_film(f3,l3,cs)
  d4a,d4b,d4c = voro_data_film(f4,l4,cs)
  d5a,d5b,d5c = voro_data_film(f5,l5,cs)
  d6a,d6b,d6c = voro_data_film(f6,l6,cs)

#voro_plot_species(name, m1, m2, m3, m4, m5, m6, l1, l2, l3, l4, l5, l6, species, cs, typ):
  voro_plot_species(name,d1a,d2a,d3a,d4a,d5a,d6a,l1,l2,l3,l4,l5,l6,0,cs,1) #All
  voro_plot_species(name,d1b,d2b,d3b,d4b,d5b,d6b,l1,l2,l3,l4,l5,l6,1,cs,1) #Cu
  voro_plot_species(name,d1c,d2c,d3c,d4c,d5c,d6c,l1,l2,l3,l4,l5,l6,2,cs,1) #Zr
  # voro_plot_case(name,d3a,d3b,d3c,l3,cs,1)
  # voro_plot_case(name,d4a,d4b,d4c,l4,cs,1)
  # voro_plot_case(name,d5a,d5b,d5c,l5,cs,1)
  # voro_plot_case(name,d6a,d6b,d6c,l6,cs,2)
  collate_ico(name,d1a,d2a,d3a,d4a,d5a,d6a,d1b,d2b,d3b,d4b,d5b,d6b,l1,l2,l3,l4,l5,l6,cs,1) #All and Cu only
  plot_chains(name,d1a,d2a,d3a,d4a,d5a,d6a,l1,l2,l3,l4,l5,l6,0,cs,1) #All
  plot_atmhisto(name,d1a,d2a,d3a,d4a,d5a,d6a,l1,l2,l3,l4,l5,l6,0,cs,1) #All
  # plot_atmhisto(name,d1b,d2b,d3b,d4b,d5b,d6b,l1,l2,l3,l4,l5,l6,3,cs,1) #Cu
  # plot_atmhisto(name,d1c,d2c,d3c,d4c,d5c,d6c,l1,l2,l3,l4,l5,l6,4,cs,1) #Zr