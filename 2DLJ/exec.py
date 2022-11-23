#!/usr/bin/env python3

import os, sys
import ovito
from ovito.io import *
from ovito import scene
from ovito.modifiers import *
from ovito.vis import *
import scipy
import matplotlib.pyplot as plt
import matplotlib

import numpy as np
import freud

import math, csv
# from ovito.vis import Viewport, CoordinateTripodOverlay
from PySide2.QtCore import *
from PySide2 import QtCore
from PySide2.QtGui import *
import lammps_logfile as lmplg
import pandas as pd
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.spatial import Voronoi, voronoi_plot_2d

pipeline = import_file('dump.atom')
# rdf1 = CoordinationAnalysisModifier(cutoff=3.0, number_of_bins=200, partial=True)
# pipeline.modifiers.append(rdf1)
data = pipeline.compute()
# XY = data.tables['coordination-rdf'].xy()
# X = np.array(XY[:,0]).reshape(len(XY[:,0]),1)
# Y = np.delete(XY,0,axis=1) #*region_volume/cell_volume
# rdf_data = np.append(X,Y,axis=1)

#with open('lj.rdf') as f:
#    lines = f.readlines()
#print(lines[4:])

a = np.loadtxt('lj.rdf')
# print(a[:,1])
# b = np.array(a)
# print(b)

fig, ax = plt.subplots(figsize=(5, 5))  # ,sharey=True)
ax.patch.set_alpha(0.2)
ax.grid(color='k', alpha=0.1, zorder=-2)
plt.setp(ax.spines.values(), linewidth=1.5)
plt.xlabel(r'Pair Separation Distance (LJ units)', fontsize=15)
plt.ylabel('g(r)', fontsize=15)
p1 = ax.plot(a[:,1], a[:,2], 'r', zorder=3, lw=2)

ax.tick_params(axis="both", direction="in", bottom=True, top=True, left=True, right=True, labelsize=14)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
# ax.legend(fontsize=12,)
plt.savefig('rdf.png', dpi=400)
# plt.show()
plt.close()

positions = data.particles['Position']
# print(data.particles['Position'])
pts = []
for xyz in positions:
    pts.append([xyz[0],xyz[1],0])
# print(points)

points=np.array(pts)
# points = np.hstack((ptsarr, np.zeros((ptsarr.shape[0], 1))))    #https://freud.readthedocs.io/en/v1.1.0/examples/module_intros/Voronoi-Voronoi.html

L = 20
box = freud.box.Box.square(L)
voro = freud.locality.Voronoi()
cells = voro.compute((box, points)).polytopes

plt.figure()
ax = plt.gca()
voro.plot(ax=ax)
ax.scatter(points[:, 0], points[:, 1], s=10, c="k")
# plt.show()
plt.savefig('voronoi.png', dpi=400)
plt.close()

# vor = Voronoi(points)
# fig, ax = plt.subplots(figsize=(5, 5))  # ,sharey=True)
# voronoi_plot_2d(vor, ax, show_vertices=False, line_colors='black', line_width=1, line_alpha=0.4, point_size=4)
# plt.ylim(-10,10)
# plt.ylim(-10,10)
# plt.tick_params(axis="both", bottom=False, top=False, left=False, right=False, labelsize=14)
# ax = plt.gca()
# ax.axes.xaxis.set_visible(False)
# ax.axes.yaxis.set_visible(False)
# plt.savefig('voronoi.png', dpi=400)
# # plt.show()
# plt.close()