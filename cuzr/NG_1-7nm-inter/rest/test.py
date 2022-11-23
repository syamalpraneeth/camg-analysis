#!/usr/bin/env python3

import os, sys
import ovito
from ovito.io import *
from ovito import scene
from ovito.modifiers import *
from ovito.vis import *
import numpy, scipy
import matplotlib.pyplot as plt
import pandas as pd
import math, csv
# from ovito.vis import Viewport, CoordinateTripodOverlay
from PySide2.QtCore import *
from PySide2 import QtCore
from PySide2.QtGui import *
import lammps_logfile as lmplg
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline

fn = 'log.file'
log = lmplg.File(fn)

S = log.get("Step")
T = log.get("Temp")
U = log.get("TotEng")
P = log.get("Press") * 100000  # bars to Pascal conversion
V = log.get("Volume") * 1e-30  # A^3 to m^3   #get volume
# enth = U + numpy.multiply(P,dV)
enth = U #2 # + numpy.multiply(P,dV)/(1.60218e-19) #Joules to eV

enthdf = pd.DataFrame({'Temp':T,'Enthalpy':enth,'Step':S})
# enthdf = enthdf.set_index('step')

delind = []
stn = 75
for hts in range(1,21):
  for x in range(stn + 1, stn + 41):
    delind.append(x)
  jump=76+41
  stn+=jump
print(delind)

new = enthdf.drop(enthdf.index[delind])
print(enthdf)
print(new)
# os._exit(1)

X=new['Temp']
Y=new['Enthalpy']

# enthalpy output to PNG
fig, ax = plt.subplots(figsize=(5, 5))  # ,sharey=True)
ax.set_facecolor('white')
plt.setp(ax.spines.values(), linewidth=1.5)
ax.grid(color='k', alpha=0.1, zorder=-2)
ax.plot(X, Y, marker=".", label='rest', mec='k')
ax.set_ylabel('Total Enthalpy (eV)', fontsize=14)
ax.set_xlabel('Temperature (K)', fontsize=14)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.legend()
plt.savefig('enth_.png', dpi=400)
plt.close()
new.to_csv('enth.csv')

# # numerical derivative of enthalpy to get C_p
# d_T    = numpy.gradient(T)
# d_enth = numpy.gradient(enth)
# C_p = d_enth/d_T
# C_p = numpy.diff(enth) / numpy.diff(T2)
# T3 = numpy.delete(T2, 0)
#
# # print(len(C_p),len(T3),len(enth))
# # f = InterpolatedUnivariateSpline(T, enth, k=1)
# # dfdx = f.derivative()
# # dydx = dfdx(x)
#
# # initial estimates
# cs_est = C_p[0]
# dc_est = abs(C_p[0] - C_p[-1])
# ts_est = 800  # I am putting this by hand
# ss_est = 1000
# dH_est = 100
# sp_est = 1000
# tp_est = 800  # I am putting this by hand
#
# df = pd.DataFrame({'T': T3, 'Cp': C_p, 'enth': numpy.delete(enth, 0)})
#
# with open('enthalpy_data-' + val + '.txt', 'w') as f:
#   dfAsString = df.to_string(header=False, index=False)
#   f.write(dfAsString)
#
# # df["model"] = cp_analytic(df["T"],cs_est,dc_est,ts_est,ss_est,dH_est,sp_est,tp_est)
# #
# # popt, pcov = curve_fit(f=cp_analytic,
# #                        xdata=df["T"],
# #                        ydata=df["Cp"],
# #                        p0=(cs_est,dc_est,ts_est,ss_est,dH_est,sp_est,tp_est)
# #                        )
# #
# # df["fit"] = cp_analytic(df["T"],popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6])
#
# # R2 = numpy.sum((df["fit"] - df["Cp"].mean())**2) / numpy.sum((df["Cp"] - df["Cp"].mean())**2)
#
# fig, ax = plt.subplots(figsize=(5, 5))  # ,sharey=True)
# plt.setp(ax.spines.values(), linewidth=1.5)
# ax.set_facecolor('white')
# ax.grid(color='k', alpha=0.1, zorder=-2)
# ax.plot(df["T"], df["Cp"], marker=".", mec='k', label=val)
# # ax.plot(df["T"],df["fit"],"r-",label='fit')
# ax.set_ylabel(r'Heat Capacity C$_p$', fontsize=14)
# ax.set_xlabel('Temperature (K)', fontsize=14)
# fig.tight_layout(rect=[0, 0.03, 1, 0.95])
# plt.legend()
# plt.savefig('cp_' + val + '.png', dpi=400)
# plt.close()