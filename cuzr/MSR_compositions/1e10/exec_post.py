#!/usr/bin/env python3

import os, sys
from ovito.io import *
from ovito.modifiers import *
import numpy
import matplotlib.pyplot as plt
import csv

from analysis import *
from plotting import *
from energy_minim2 import *

src='/home/mj0054/Documents/work/simulations/projects/bulk/msr/'
name='CuZr MG'
f1,l1=src+'46-54/1e10/bigger/data/dump.msr_20','46-54'
f2,l2=src+'50-50/1e10/bigger/data/dump.msr_20','50-50'
f3,l3=src+'55-45/1e10/bigger/data/dump.msr_20','55-45'

cs=0

#"""
#pe/atom AS_PREPARED
for f in glob.glob('tmp.pote_coll_' + name.replace(" ", "_") + '_asp_*'):
  if os.path.exists(f): os.remove(f)
ovito_pote(name,f1,l1,cs,0,1)
ovito_pote(name,f2,l2,cs,0,1)
ovito_pote(name,f3,l3,cs,0,1)
pe_plot_species(name,l1,l2,l3,cs,1)
#"""
# os._exit(1)

#"""
#rdf AS_PREPARED
ovito_rdf(f1,l1,cs,0)
ovito_rdf(f2,l2,cs,0)
ovito_rdf(f3,l3,cs,0)

plot_rdf(name,l1,0,cs,1)
plot_rdf(name,l2,0,cs,1)
plot_rdf(name,l3,0,cs,1)
plot_allrdf(name,l1,l2,l3,cs,1)
#"""

# os._exit(1)

#voronoi AS_PREPARED
d1a,d1b,d1c = voro_data_bulk(f1,l1,0)
#  d2a,d2b,d2c = d1a,d1b,d1c
#  d3a,d3b,d3c = d1a,d1b,d1c
#  d4a,d4b,d4c = d1a,d1b,d1c
d2a,d2b,d2c = voro_data_bulk(f2,l2,cs)
d3a,d3b,d3c = voro_data_bulk(f3,l3,cs)

voro_plot_species(name, d1a, d2a, d3a, l1, l2, l3, 0, cs, 1)  # All
voro_plot_species(name, d1b, d2b, d3b, l1, l2, l3, 1, cs, 1)  # Cu
voro_plot_species(name, d1c, d2c, d3c, l1, l2, l3, 2, cs, 1)  # Zr

# voro_plot_case(name,d1a,d1b,d1c,l1,0,1)
# voro_plot_case(name,d2a,d2b,d2c,l2,0,1)
# voro_plot_case(name,d3a,d3b,d3c,l3,cs,1)
# voro_plot_case(name,d4a,d4b,d4c,l4,cs,1)

# collate_ico(name, d1a, d2a, d3a, d1b, d2b, d3b, l1, l2, l3, cs, 1)  # All and Cu only
# plot_chains(name, d1a, d2a, d3a, l1, l2, l3, 0, cs, 1)  # All
plot_atmhisto(name, d1a, d2a, d3a, l1, l2, l3, 0, cs, 1)  # All
# plot_atmhisto(name,d1b,d2b,d3b,d4b,d5b,d6b,l1,l2,l3,l4,3,cs,1) #Cu
# plot_atmhisto(name,d1c,d2c,d3c,d4c,d5c,d6c,l1,l2,l3,l4,4,cs,1) #Zr
