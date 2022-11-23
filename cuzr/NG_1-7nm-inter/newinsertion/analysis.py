##!/usr/bin/env python3

import os, sys
import ovito
from ovito.io import *
from ovito import scene
from ovito.modifiers import *
from ovito.vis import *
import numpy, scipy

import math, csv
# from ovito.vis import Viewport, CoordinateTripodOverlay
from PySide2.QtCore import *
from PySide2 import QtCore
from PySide2.QtGui import *
from sort_voro import *
import lammps_logfile as lmplg
import pandas as pd
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline

def rem_surf(pipeline,case):
  # if case!= 'NG':
  #   string = '(Position.Z>50 || Position.Z<25) || ((abs(Position.X)>70) || (abs(Position.Y)>70))'
  #   pipeline.modifiers.append(ExpressionSelectionModifier(expression=string))
  #   pipeline.modifiers.append(DeleteSelectedModifier())
  #   data = pipeline.compute()

  iter = 0
  while iter < 2:
    surf = ConstructSurfaceModifier(radius=3, smoothing_level=10, select_surface_particles=True)
    pipeline.modifiers.append(surf)
    # pipeline.modifiers.append(ExpressionSelectionModifier(expression='SurfaceDistance<3'))
    data = pipeline.compute()
    surf.vis.enabled = False
    pipeline.modifiers.append(DeleteSelectedModifier())
    iter += 1
  #cluster analysis
  pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=3.5, sort_by_size=True))
  #choose largest cluster
  pipeline.modifiers.append(ExpressionSelectionModifier(expression='Cluster!=1'))
  pipeline.modifiers.append(DeleteSelectedModifier())
  return pipeline;

# This function is called by OVITO on every viewport update.
def scale_render(args):
  print('Hello from inside scale render')
  # Parameters:
  bar_length = 100  # Simulation units (e.g. Angstroms)
  bar_color = QColor(0, 0, 0)
  label_text = "{} nm".format(bar_length / 10)
  label_color = QColor(255, 255, 255)

  if args.is_perspective:
      raise Exception("This overlay only works with non-perspective viewports.")

  # Compute length of bar in screen space
  screen_length = args.project_size((0,0,0), bar_length)

  # Define geometry of bar in screen space
  height = 0.07 * args.painter.window().height()
  margin = 0.45 * args.painter.window().height()
  rect = QRectF(0.45*margin, margin, screen_length, height)

  # Render bar rectangle
  args.painter.fillRect(rect, bar_color)

  # Render text label
  font = args.painter.font()
  font.setPixelSize(height)
  args.painter.setFont(font)
  args.painter.setPen(QPen(label_color))
  args.painter.drawText(rect, Qt.AlignCenter, label_text)

def ovito_voro(name, fil, val, cs, typ, species, delsub):
  if species == 0:
    case = 'All'
  elif species == 1 or species == 3:
    case = 'Cu'
  elif species == 2 or species == 4:
    case = 'Zr'
  else:
    raise ValueError('species value only 0,1 or 2!')
  if cs == 0:
    csstr = ' '
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
  fn = fil.split("/")[-1].split(".")[1]

  pipeline = import_file(fil)
  # Set atomic radii (required for polydisperse Voronoi tessellation).
  atom_types = pipeline.source.data.particles_.particle_types_
  for x in atom_types:
    if x % 2 == 0:
      atom_types.type_by_id(x).radius = 1.55
    else:
      atom_types.type_by_id(x).radius = 1.35

  data = pipeline.compute()
  at_typ_ind = data.particles['Particle Type']
  print(fn, 'Voronoi:', csstr)
  print('Initial num atoms:', len(at_typ_ind))

  # Enable core and Interface viz
  if delsub == 1:
    color_mod = ColorCodingModifier(property='i_int', start_value=1, end_value=3,
                                    gradient=ColorCodingModifier.Magma())
    pipeline.modifiers.append(color_mod)
    slice_mod = SliceModifier(normal=(0,0,1))
    pipeline.modifiers.append(slice_mod)

  tripod = CoordinateTripodOverlay()
  tripod.size = 0.07
  tripod.style = CoordinateTripodOverlay.Style.Solid

  cell_vis = pipeline.source.data.cell.vis
  cell_vis.enabled = False

  # Viz system
  if cs == 0 and species == 0:
    pipeline.add_to_scene()
    vp = Viewport()
    vp.type = Viewport.Type.Ortho
    vp.zoom_all()
    # if delsub == 1:
    #   vp.camera_pos = (-0.770151, 0.773092, 83.9268)
    #   vp.camera_dir = (-0.999997, 0.00251444, 4.57373e-15)
    #   vp.fov = 120.403
    # else:
    #   vp.camera_pos = (36.8437, 36.8438, 36.845)
    #   vp.camera_dir = (0, 1, 0)
    #   vp.fov = 83.3431

    if val == '60meV':
      vp.overlays.append(PythonViewportOverlay(function = scale_render))

    vp.overlays.append(tripod)
    vp.render_image(size=(800, 600), filename=fn + '_' + tag + "_ortho.png", background=(1, 1, 1))
    pipeline.remove_from_scene()

  if delsub == 1: slice_mod.enabled = False  # disable slice, which was used for the visual
  if delsub == 1 and cs==0 and species ==0:  ovito_strain(name, pipeline,fil,typ,val) #and val != 'NG' and '1e10' in fil.split("/"):


  # Set up the Voronoi analysis modifier.
  voro = VoronoiAnalysisModifier(
    compute_indices=True,
    use_radii=True,
    edge_threshold=0.1
  )
  pipeline.modifiers.append(voro)

  if delsub == 1:
    # Delete all substrate atoms
    pipeline.modifiers.append(ExpressionSelectionModifier(expression='ParticleType < 3'))
    data = pipeline.compute()
    pipeline.modifiers.append(DeleteSelectedModifier())
    data = pipeline.compute()

    atom_types = data.particles['Particle Type']
    print('Num non substrate atoms:', len(atom_types))

    #Retain only non-surface atoms
    pipeline = rem_surf(pipeline,val)
    data = pipeline.compute()

    # Retain only core or Interface atoms
    if cs == 1:
      pipeline.modifiers.append(ExpressionSelectionModifier(expression='i_int != 2'))
      pipeline.modifiers.append(DeleteSelectedModifier())
      data = pipeline.compute()
    elif cs == 2:
      pipeline.modifiers.append(ExpressionSelectionModifier(expression='i_int != 3'))
      pipeline.modifiers.append(DeleteSelectedModifier())
      data = pipeline.compute()
    atom_types = data.particles['Particle Type']
    print('Number case: ', csstr, 'atoms', len(atom_types))

  if species == 0:
    data = pipeline.compute()
    atom_types = data.particles['Particle Type']
    print('Num slab atoms: ', len(atom_types))

  else:
    pipeline.modifiers.append(ExpressionSelectionModifier(expression='ParticleType!=' + str(species)))
    pipeline.modifiers.append(DeleteSelectedModifier())
    data = pipeline.compute()
    atom_types = data.particles['Particle Type']
    print('Num slab ', case, 'atoms: ', len(atom_types))

  # Access computed Voronoi indices.
  # This is an (N) x (M) array, where M is the maximum face order.
  voro_indices = data.particles['Voronoi Index']
  atom_volume = data.particles['Atomic Volume']
  atom_types = data.particles['Particle Type']
  e_del = []
  data = pipeline.compute()
  atom_types = data.particles['Particle Type']
  print('Final num atoms:', len(atom_types))
  tot_selec_atms = len(voro_indices)

  if delsub ==1:
    slice_mod2 = SliceModifier(normal=(0,0,1))
    pipeline.modifiers.append(slice_mod2)
    # slice_mod2.enabled = True  # disable slice, which was used for the visual

  # Viz Slab
  pipeline.add_to_scene()
  vp = Viewport()
  vp.type = Viewport.Type.Ortho
  # vp.camera_pos = (-0.770151, 0.773092, 83.9268)
  # vp.camera_dir = (-0.999997, 0.00251444, 4.57373e-15)
  # vp.fov = 120.403
  vp.zoom_all()
  vp.overlays.append(tripod)
  if cs == 0 and species==0:
    vp.render_image(size=(800, 600), filename=fn + '_' + tag + "_ortho_del.png", background=(1, 1, 1))
  elif cs > 0 and species==0:
    vp.render_image(size=(800, 600), filename=fn + '_' + csstr + '_' + tag + "_ortho_del.png", background=(1, 1, 1))
  else: pass
  print('Voro_indices', numpy.shape(voro_indices), str(csstr))
  if delsub ==1: slice_mod2.enabled = False  # disable slice, which was used for the visual
  pipeline.remove_from_scene()

  choose_ico = ExpressionSelectionModifier(expression='(Coordination==12 && VoronoiIndex.5 == 12)')
  pipeline.modifiers.append(choose_ico)
  # choose_icolike = ExpressionSelectionModifier(expression='(Coordination==12 && VoronoiIndex.5 == 12) ||'
  #                                             '(VoronoiIndex.3 == 0 && VoronoiIndex.4 == 1 && VoronoiIndex.5 == 10 && VoronoiIndex.6 < 5 &&'
  #                                             'MaxFaceOrder==6) ||'
  #                                             '(VoronoiIndex.3 == 0 && VoronoiIndex.4 == 2 && VoronoiIndex.5 == 8 && VoronoiIndex.6 < 5 &&'
  #                                             'MaxFaceOrder==6) ')
  # pipeline.modifiers.append(choose_icolike)

  pipeline.modifiers.append(InvertSelectionModifier())
  pipeline.modifiers.append(DeleteSelectedModifier())
  data = pipeline.compute()
  bonds = CreateBondsModifier(cutoff=3.5)
  bonds.vis.shading = BondsVis.Shading.Flat
  data.particles.vis.enabled = False
  pipeline.modifiers.append(bonds)
  data = pipeline.compute()
  pipeline.modifiers.append(ClusterAnalysisModifier(neighbor_mode=ClusterAnalysisModifier.NeighborMode.Bonding, cutoff=3.5, sort_by_size=True))
  data = pipeline.compute()
  # Cluster length of icosahedron connected polyhedra
  clus = data.tables['clusters']
  # meanSize = numpy.mean(clus['Cluster Size'][:])
  # maxClusterSize = numpy.max(clus['Cluster Size'][:])
  # n, bins, patches = plt.hist(clus['Cluster Size'][:], 25, facecolor='g', histtype='step', alpha=0.75)
  csiz_arr = clus['Cluster Size'][:]
  print('size of cluster size array is:',len(csiz_arr))

  pipeline.add_to_scene()
  vp = Viewport()
  vp.type = Viewport.Type.Top
  vp.zoom_all()
  # vp.camera_pos = (0, 0, 0)
  # vp.camera_dir = (0, 0, -1)
  vp.fov = 90
  if cs == 0 and species==0:
    vp.render_image(size=(800, 800), filename=fn +'_chains-il_'+ case + '_' + tag + ".png", background=(0, 0, 0))
  elif cs > 0 and species==0:
    vp.render_image(size=(800, 800), filename=fn +'_chains-il_' + case + '_' + csstr + '_' + tag + ".png", background=(0, 0, 0),
                    frame=8)
  else: pass
  pipeline.remove_from_scene()

  return [voro_indices, numpy.array(atom_volume), csiz_arr, tot_selec_atms,pipeline];

def voro(name, fil, val, cs, typ, species, delsub):
  fn = fil.split("/")[-1].split(".")[1]
  siz = 8
  sz = siz - 1
  [a, vols, csa, N,pipeline] = ovito_voro(name, fil, val, cs, typ, species, delsub)

  #voro to MRO viz
  if cs == 0:
    csstr = 'All'
  elif cs == 1:
    csstr = 'Core'
  elif cs == 2:
    csstr = 'Interface'
  else:
    raise ValueError('cs value only 0, 1 or 2!')
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

  pipeline.add_to_scene()
  ovito.scene.save("ico-chains_"+fn+'_'+tag+'_'+ csstr +".ovito")
  pipeline.remove_from_scene()

  #Voro Histogram snippet
  ca = numpy.ascontiguousarray(a).view([('', a.dtype)] * a.shape[1])
  unique, indices, inverse = numpy.unique(ca, return_index=True, return_inverse=True)

  counts = numpy.bincount(inverse)
  sort_indices = numpy.argsort(counts)[::-1]
  sort2 = sorted(sort_indices[0:sz])

  print('Unique Voronoi Indices', numpy.shape(indices))

  height = []
  bars = []
  ht = []
  brs = []

  ind = a[indices[sort2]]
  cnts = counts[sort2]
  for i in range(0, len(cnts)):
    h = 100.0 * (cnts[i]) / (N)
    ht.append(float('%.3f' % h))
    brs.append('<' + str(ind[i, 2]) + ' ' + str(ind[i, 3]) + ' ' + str(ind[i, 4]) + ' ' + str(ind[i, 5]) + '>')

  hsort = datclass()
  fullmat = a[indices[sort_indices]]
  fullcounts = counts[sort_indices]
  #    print('full mat',numpy.shape(ca),numpy.shape(vols))
  #    print((a[:,2:6]))
  # hsort, vol_nico = sort_voro(a[:, 2:6], vols)
  hsort, vol_nico = sort_voro(a, vols)
  #    print('From inside voro: ',numpy.shape(vol_nico),'for',case)
  #    return ht,hsort,brs,b,vol_nico,N;
  return ht, brs, vols, N, hsort, vol_nico, csa;

class datclass:
  def __init__(self):
    self.a = [1, 2]
    self.b = [2, 3]
    self.c = [3, 4]
    self.d = 1
    self.h = [0, 0, 0, 0]
    self.v = 1
    self.csa = [3, 4]

def voro_data_bulk(name, fil, val, cs, typ):
  d1, d2, d3 = datclass(), datclass(), datclass()
  d1.a, d1.b, d1.c, d1.d, d1.h, d1.v, d1.csa = voro(name, fil, val, cs, typ, 0, 0)
  # d2.a, d2.b, d2.c, d2.d, d2.h, d2.v, d2.csa = voro(name, fil, val, cs, typ, 1, 0)
  # d3.a, d3.b, d3.c, d3.d, d3.h, d3.v, d3.csa = voro(name, fil, val, cs, typ, 2, 0)
  #    print('From inside voro_data_bulk: ',numpy.shape(d1.v),numpy.shape(d2.v),numpy.shape(d3.v))
  return [d1, d2, d3];

def voro_data_film(name, fil, val, cs, typ):
  d1, d2, d3 = datclass(), datclass(), datclass()
  d1.a, d1.b, d1.c, d1.d, d1.h, d1.v, d1.csa = voro(name, fil, val, cs, typ, 0, 1)
  # d2.a, d2.b, d2.c, d2.d, d2.h, d2.v, d2.csa = voro(name, fil, val, cs, typ, 3, 1)
  # d3.a, d3.b, d3.c, d3.d, d3.h, d3.v, d3.csa = voro(name, fil, val, cs, typ, 4, 1)
  return [d1, d2, d3];

def ovito_pote(name,fil, case, cs, delsub,typ):
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
  fn = fil.split("/")[-1].split(".")[1]

  pipeline = import_file(fil)
  # Set atomic radii (required for polydisperse Voronoi tessellation).
  atom_types = pipeline.source.data.particles_.particle_types_
  for x in atom_types:
    if x % 2 == 0:
      atom_types.type_by_id(x).radius = 1.55
    else:
      atom_types.type_by_id(x).radius = 1.35

  data = pipeline.compute()
  at_typ_ind = data.particles['Particle Type']
  print(fn, 'PE/Atom :', csstr)
  print('Initial num atoms:', len(at_typ_ind))

  tripod = CoordinateTripodOverlay()
  tripod.size = 0.07
  tripod.style = CoordinateTripodOverlay.Style.Solid

  if delsub == 1:
    # Delete all substrate atoms
    pipeline.modifiers.append(ExpressionSelectionModifier(expression='ParticleType < 3'))
    data = pipeline.compute()
    pipeline.modifiers.append(DeleteSelectedModifier())
    data = pipeline.compute()
    color_mod = ColorCodingModifier(property='i_int', start_value=1, end_value=3,
                                    gradient=ColorCodingModifier.Magma())
    pipeline.modifiers.append(color_mod)
    atom_types = data.particles['Particle Type']
    print('Num non-substrate atoms:', len(atom_types))

    # Viz deleted substrate atoms
    if cs == 0:
      pipeline.add_to_scene()
      vp = Viewport()
      vp.type = Viewport.Type.Perspective
      vp.zoom_all()
      # vp.camera_pos = (308.835, -411.782, 438.376)
      # vp.camera_dir = (-0.4819, 0.6424, -0.5956)
      # vp.fov = math.radians(35)
      vp.overlays.append(tripod)
      vp.render_image(size=(800, 600), filename=fn +"_"+tag+ "_pers.png", background=(1, 1, 1))
      pipeline.remove_from_scene()

    # Retain only non surface atoms in CAMGs
    pipeline = rem_surf(pipeline,case)
    data = pipeline.compute()

    # Retain only core or Interface atoms
    if cs == 1:
      pipeline.modifiers.append(ExpressionSelectionModifier(expression='i_int != 2'))  # Select non-core atoms
      data = pipeline.compute()
    elif cs == 2:
      pipeline.modifiers.append(ExpressionSelectionModifier(expression='i_int != 3'))  # Select non-Interface atoms
      data = pipeline.compute()
    pipeline.modifiers.append(DeleteSelectedModifier())
    data = pipeline.compute()
    atom_types = data.particles['Particle Type']
    print('Num case:' + csstr + ' atoms:', len(atom_types))

    pipeline.add_to_scene()
    vp = Viewport()
    vp.type = Viewport.Type.Perspective
    # vp.camera_pos = (308.835, -411.782, 438.376) #(247.653, -179.174, 375.894) #
    # vp.camera_dir = (-0.4819, 0.6424, -0.5956) # (-0.62888, 0.4386, -0.6412) #
    # vp.fov = math.radians(35)
    vp.zoom_all()
    vp.overlays.append(tripod)
    if cs == 0:
      vp.render_image(size=(800, 600), filename=fn + "_"+tag+ "_pers_del.png", background=(1, 1, 1))
    elif cs > 0:
      vp.render_image(size=(800, 600), filename=fn +"_"+tag+ '_' + csstr + "_pers_del.png", background=(1, 1, 1))
    pipeline.remove_from_scene()

    atom_types = data.particles['Particle Type']
    print('Num slab atoms:', len(atom_types))

  atom_types = data.particles['Particle Type']
  print('Final num atoms:', len(atom_types))

  # Output PE/atom histogram data to tmp file
  modifier = HistogramModifier(bin_count=100, property='f_fpe')
  pipeline.modifiers.append(modifier)
  if cs == 0: histfil = 'tmp.histo_' + case + '_' +tag
  if cs > 0: histfil = 'tmp.histo_' + case + '_' + csstr + '_' +tag
  export_file(pipeline, histfil, "txt/table", key="histogram[f_fpe]")

  dtfil = fil.replace(fil.split("/")[-1],'')[0:-5]+'data.clus_anneal_'+case
  f = open(dtfil, 'r')
  dat = f.readlines()
  clus=dat[2].split(" ")[0]
  f.close()
  potarr = numpy.array(data.particles['f_fpe'])
  print(potarr.sum(axis=0))
  pote= potarr.sum(axis=0)/len(atom_types)
  # rows=zip(str(pote),case)
  rows = [str(pote),',',case,',',clus]
  # print(rows)
  fil2 = 'tmp.pote_coll_' + name.replace(" ", "_") + '_' + tag + '_' + csstr
  with open(fil2, 'a') as f:
    # writer = csv.writer(f)
    for row in rows:
    #   print(row)
    #   writer.writerow(row)
      f.write(row)
    f.write('\n')

  # Viz slice of copper atoms
  pipeline.modifiers.append(ExpressionSelectionModifier(expression='ParticleType==2 || ParticleType ==4'))
  pipeline.modifiers.append(DeleteSelectedModifier())

  # pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff=5.0, number_of_bins=200))
  # data = pipeline.compute()
  # rdffil = 'tmp.rdf_' + case
  # numpy.savetxt(rdffil, data.tables['coordination-rdf'].xy())

  z, t = 20, 20
  pipeline.modifiers.append(SliceModifier(distance=z, normal=(0, 0, 1), slab_width=t))
  ol = TextLabelOverlay(
    text='Slice @ Z=' + str(z) + ', Slab width = ' + str(t) + r'$\AA$ for '+ case,
    alignment=Qt.AlignHCenter ^ Qt.AlignTop,
    offset_y=0.020,
    font_size=0.03,
    text_color=(0, 0, 0))

  pipeline.add_to_scene()
  vp = Viewport()
  vp.zoom_all()
  vp.overlays.append(ol)
  vp.type = Viewport.Type.Top
  # if delsub == 1:
  #   vp.camera_pos = (0,0,57) #(0, 0, 0)
  # else:
  #   vp.camera_pos = (36.84, 36.84, 36.845)
  # vp.camera_dir = (0, 0, -1)
  # if delsub == 1:
  #   vp.fov = 123
  # else:
  #   vp.fov = 83.34

  if cs == 0:
    vp.render_image(size=(800, 800), filename=fn + '_' + tag + "_slice.png", background=(1, 1, 1))
  elif cs > 0:
    vp.render_image(size=(800, 800), filename=fn + '_' + tag + "_slice_" + csstr + ".png", background=(1, 1, 1))
  pipeline.remove_from_scene()

  color_mod = ColorCodingModifier(property='f_fpe', start_value=-3.79, end_value=-2.73)
  pipeline.modifiers.append(color_mod)

  ol2 = ColorLegendOverlay(
    modifier=color_mod,
    title='P.E./ Cu atom:',
    alignment=Qt.AlignRight ^ Qt.AlignBottom,
    orientation=Qt.Vertical,
    offset_y=0.04,
    font_size=0.08,
    format_string='%.2f eV')

  pipeline.add_to_scene()
  vp = Viewport()
  vp.overlays.append(ol)
  vp.overlays.append(ol2)
  vp.type = Viewport.Type.Top
  if delsub == 1:
    vp.camera_pos = (0,0,57) #(0, 0, 0)
  else:
    vp.camera_pos = (36.84, 36.84, 36.845)
  vp.camera_dir = (0, 0, -1)
  if delsub == 1:
    vp.fov = 100
  else:
    vp.fov = 83.34
  if cs == 0:
    vp.render_image(size=(800, 800), filename=fn + '_' + tag + "_slice-color.png", background=(1, 1, 1))
  elif cs > 0:
    vp.render_image(size=(800, 800), filename=fn + '_' + tag + "_slice-color_" + csstr + ".png", background=(1, 1, 1))

#  print("works")
  pipeline.remove_from_scene()

  return []

def ovito_rdf(name, fil, case, cs, delsub, typ):
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
  fn = fil.split("/")[-1].split(".")[1]

  pipeline = import_file(fil)
  # Set atomic radii (required for polydisperse Voronoi tessellation).
  atom_types = pipeline.source.data.particles_.particle_types_
  for x in atom_types:
    if x % 2 == 0:
      atom_types.type_by_id(x).radius = 1.55
    else:
      atom_types.type_by_id(x).radius = 1.35

  data = pipeline.compute()
  at_typ_ind = data.particles['Particle Type']
  print(fn, 'RDF:', csstr)
  print('Initial num atoms:', len(at_typ_ind))

  if delsub == 1:
    pipeline.modifiers.append(ExpressionSelectionModifier(expression='ParticleType < 3'))
    data = pipeline.compute()
    pipeline.modifiers.append(DeleteSelectedModifier())  # Delete all substrate atoms
    data = pipeline.compute()

    # Retain only non surface atoms in CAMGs
    pipeline = rem_surf(pipeline,case)
    data = pipeline.compute()
    at_typ_ind = data.particles['Particle Type']
    print('Num slab atoms', len(at_typ_ind))

    if cs == 1:
      pipeline.modifiers.append(ExpressionSelectionModifier(expression='i_int != 2'))  # Select non-core atoms
      data = pipeline.compute()
    elif cs == 2:
      pipeline.modifiers.append(ExpressionSelectionModifier(expression='i_int != 3'))  # Select non-Interface atoms
      data = pipeline.compute()
    pipeline.modifiers.append(DeleteSelectedModifier())
    data = pipeline.compute()

    at_typ_ind = data.particles['Particle Type']
    print('Num slab' ,case, csstr, 'atoms:', len(at_typ_ind))

  cell_volume = data.cell.volume

  dtfil = fil.replace(fil.split("/")[-1],'')[0:-5]+'data.clus_anneal_'+case
  print(dtfil)
  f = open(dtfil, 'r')
  dat = f.readlines()
  clus=dat[2].split(" ")[0]
  f.close()

  if cs==0:
    atoms = data.particles['Particle Type']
    density = (63.546*len(atoms[atoms==3]) + 91.224*len(atoms[atoms==4]))/cell_volume #not region volume
    rdensity = density/((63.546*67584+91.224*68112)/(134.18)**3)
    # print(atoms[1])
    print('Tot atoms',len(atoms),'Cu atoms',len(atoms[atoms==3]),'Zr atoms',len(atoms[atoms==4]))
    rows = [str(rdensity),',',case,',',clus]
    fil2 = 'tmp.density_coll_'+name.replace(" ", "_")+'_'+tag+'_All'
    with open(fil2, 'a') as f:
      for row in rows:
        f.write(row)
      f.write('\n')

  #Volume calculation using surface mesh
  remesh = ConstructSurfaceModifier(radius=3, smoothing_level=10, select_surface_particles=True,identify_regions=True)
  pipeline.modifiers.append(remesh)
  # pipeline.modifiers.append(ConstructSurfaceModifier(radius=3, smoothing_level=10, select_surface_particles=True,identify_regions=True))

  data = pipeline.compute()

  region_volume = data.attributes["ConstructSurfaceMesh.filled_volume"]
  frac = data.attributes["ConstructSurfaceMesh.filled_fraction"]
  poro = (1 - frac) * 100
  print(cell_volume,region_volume,frac,poro, case)

  if cs==0:
    poro = (1-frac)*100
    rows = [str(poro),',',case,',',clus]
    fil2 = 'tmp.porosity_coll_'+name.replace(" ", "_")+'_'+tag+'_All'
    with open(fil2, 'a') as f:
      for row in rows:
        f.write(row)
      f.write('\n')

  rdf1 = CoordinationAnalysisModifier(cutoff=10.0, number_of_bins=200, partial=True)
  pipeline.modifiers.append(rdf1)
  data = pipeline.compute()

  # print(data.tables['coordination-rdf'].xy()) #*region_volume/cell_volume)
  XY = data.tables['coordination-rdf'].xy()
  X = numpy.array(XY[:,0]).reshape(len(XY[:,0]),1)
  Y = numpy.delete(XY,0,axis=1)*region_volume/cell_volume
  rdf_data = numpy.append(X,Y,axis=1)

  # print(X)
  # print(rdf_data)
  rdffil = 'tmp.rdf_' + case +'_'+csstr +'_'+ tag
  # print(data.tables['coordination-rdf'].y.component_names)
  # numpy.savetxt(rdffil, data.tables['coordination-rdf'].xy())
  numpy.savetxt(rdffil, rdf_data)

  rows2 = [str(frac), ',', case, ',', clus]
  fil2 = 'tmp.filledfrac_coll_' + name.replace(" ", "_") + '_' + tag + '_All'
  with open(fil2, 'a') as f:
    for row in rows2:
      f.write(row)
    f.write('\n')

  rdf1.enabled = False
  remesh.enabled = False
  return []

def union(p, q):  # function to give union of two lists
  ap = [x for x in q if x not in p]
  un = p + ap
  return un;

def inter(p, q):  # function to give intersection of two lists
  intr = [x for x in q if x in p]
  return intr;

def reassemble(h1, b1, h2, b2, h3, b3, h4, b4, h5, b5, h6, b6, h7, b7, h8, b8):  # rearrange all b entries
  un1 = union(b1, b2)
  un2 = union(b3, un1)
  un3 = union(b4, un2)
  un4 = union(b5, un3)
  un5 = union(b6, un4)
  un6 = union(b7, un5)
  un7 = union(b8, un6)

  in1 = inter(b1, b2)
  in2 = inter(b3, in1)
  in3 = inter(b4, in2)
  in4 = inter(b5, in3)
  in5 = inter(b6, in4)
  in6 = inter(b7, in5)
  in7 = inter(b8, in6)

  ncdf = [x for x in un7 if x not in in7]
  b = in7 + ncdf

  H1, H2, H3, H4, H5, H6, H7, H8 = [], [], [], [], [], [], [], []
  for x in b:
    if x in b1:
      ind = b1.index(x)
      H1.append(h1[ind])
    else:
      H1.append(0)
    if x in b2:
      ind = b2.index(x)
      H2.append(h2[ind])
    else:
      H2.append(0)
    if x in b3:
      ind = b3.index(x)
      H3.append(h3[ind])
    else:
      H3.append(0)
    if x in b4:
      ind = b4.index(x)
      H4.append(h4[ind])
    else:
      H4.append(0)
    if x in b5:
      ind = b5.index(x)
      H5.append(h5[ind])
    if x in b6:
      ind = b6.index(x)
      H6.append(h6[ind])
    if x in b7:
      ind = b7.index(x)
      H7.append(h7[ind])
    if x in b8:
      ind = b8.index(x)
      H8.append(h8[ind])
    else:
      H5.append(0)

  return [H1, b, H2, b, H3, b, H4, b, H5, b, H6, b, H7, b, H8, b];

def collate_ico(name, o1, o2, o3, o4, o5, o6, o7, o8, p1, p2, p3, p4, p5, p6, p7, p8, l1, l2, l3, l4, l5, l6, l7 ,l8, cs, typ):
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

  i_ico = "<0 0 12 0>"
  a1 = [o1.a[o1.b.index(i_ico)], o2.a[o2.b.index(i_ico)], o3.a[o3.b.index(i_ico)], o4.a[o4.b.index(i_ico)], o5.a[o5.b.index(i_ico)], o6.a[o6.b.index(i_ico)], o7.a[o7.b.index(i_ico)], o8.a[o8.b.index(i_ico)]]
  #a2 = [p1.a[p1.b.index(i_ico)], p2.a[p2.b.index(i_ico)], p3.a[p3.b.index(i_ico)], p4.a[p4.b.index(i_ico)], p5.a[p5.b.index(i_ico), p6.a[p6.b.index(i_ico), p7.a[p7.b.index(i_ico), p8.a[p8.b.index(i_ico)]]
  #swtiched off voro analysis for cu atoms
  a3 = [l1, l2, l3, l4, l5, l6, l7, l8]

#  rows = zip(a1, a2, a3)
  rows = zip(a1, a3)  #swtiched off voro analysis for cu atoms
  fil2 = 'tmp.ico_coll_' + name.replace(" ", "_") + '_' + tag + '_' + csstr
  with open(fil2, 'w') as f:
    writer = csv.writer(f)
    for row in rows:
      writer.writerow(row)

  b1 = [o1.h[0], o2.h[0], o3.h[0], o4.h[0], o5.h[0], o6.h[0], o7.h[0], o8.h[0]]
  #b2 = [p1.h[0], p2.h[0], p3.h[0], p4.h[0], p5.h[0], p6.h[0], p7.h[0], p8.h[0]]
  #rows2 = zip(b1, b2, a3)
  rows2 = zip(b1, a3)
  fil2 = 'tmp.ico-like_coll_' + name.replace(" ", "_") + '_' + tag + '_' + csstr
  with open(fil, 'w') as f:
    writer = csv.writer(f)
    for row in rows2:
      writer.writerow(row)

  c1 = [sum(o1.h[1:]), sum(o2.h[1:]), sum(o3.h[1:]), sum(o4.h[1:]), sum(o5.h[1:]), sum(o6.h[1:]), sum(o7.h[1:]), sum(o8.h[1:])]
  #c2 = [sum(p1.h[1:]), sum(p2.h[1:]), sum(p3.h[1:]), sum(p4.h[1:]), sum(p5.h[1:]), sum(p6.h[1:]), sum(p7.h[1:]), sum(p8.h[1:])]
  #rows3 = zip(c1, c2, a3)
  rows3 = zip(c1, a3)
  fil2 = 'tmp.nico_coll_' + name.replace(" ", "_") + '_' + tag + '_' + csstr
  with open(fil, 'w') as f:
    writer = csv.writer(f)
    for row in rows3:
      writer.writerow(row)

def ovito_strain(name, pipeline, fil, typ, case):
  if typ == 1:
    comm = 'As Prepared'
    tag = 'asp'
  elif typ == 2:
    comm = 'Annealed'
    tag = 'ann'
  else:
    raise ValueError('typ value only 1 or 2!')
  fn = fil.split("/")[-1].split(".")[1]

  file='/home/mj0054/Documents/work/simulations/projects/nanoglass/50-50/1e10/segregated/random/'+name.split(' ')[-1]+'/newinsertion/refconfig_'+case+'.txt'
  print('strain ref file: ',file)
  pip2 = import_file(file)
  strmod = AtomicStrainModifier(cutoff=3.8)

  strmod.reference = pip2.source
  strmod.reference.load(file)
  pipeline.modifiers.append(strmod)

  slice_mod2 = SliceModifier(normal=(0,0,1))
  pipeline.modifiers.append(slice_mod2)
  color_mod2 = ColorCodingModifier(property='Shear Strain', start_value=0, end_value=3,
                                    gradient=ColorCodingModifier.Jet())
  pipeline.modifiers.append(color_mod2)

  tripod = CoordinateTripodOverlay()
  tripod.size = 0.07
  tripod.style = CoordinateTripodOverlay.Style.Solid

  # Viz system
  pipeline.add_to_scene()
  vp = Viewport()
  vp.type = Viewport.Type.Ortho
  vp.zoom_all()
  # vp.camera_pos = (-0.770151, 0.773092, 83.9268)
  # vp.camera_dir = (-0.999997, 0.00251444, 4.57373e-15)
  # vp.fov = 120.403
  vp.overlays.append(tripod)
  overlay = ColorLegendOverlay(
    modifier = color_mod2,
    title = 'Strain',
    alignment = Qt.AlignRight ^ Qt.AlignBottom,
    orientation = Qt.Vertical,
    offset_y = 0.06,
    font_size = 0.12,
    label1 ="≥3",
    format_string = '%1.0f')
  vp.overlays.append(overlay)
  vp.render_image(size=(800, 600), filename= fn + '_' + tag + "_strain_ortho.png", background=(1, 1, 1)) #alpha=True)
  pipeline.remove_from_scene()
  color_mod2.enabled = False
  slice_mod2.enabled = False

  # XY = data.tables['coordination-rdf'].xy()
  data = pipeline.compute()
  # print(data.particles)
  strain = data.particles['Shear Strain']
  print('length of data',len(strain))

  dtfil = fil.replace(fil.split("/")[-1],'')[0:-5]+'data.clus_anneal_'+case
  f = open(dtfil, 'r')
  dat = f.readlines()
  clus=dat[2].split(" ")[0]
  f.close()

  rows=[]
  for st in [1,1.5,2,2.5]:
    highstrainpop = len(strain[strain>st])
    print(len(strain))
    print(strain[strain>=st])
    highstrainpct = highstrainpop*100/len(strain)
    print('% atoms with shear strain above '+str(st)+' : '+str(highstrainpct))
    # rows = [case], str(highstrainpct)]
    rows.append(str(highstrainpct)+',')

  rows.append(case+',')
  rows.append(clus)

  fil2 = 'tmp.strain_coll_'+name.replace(" ", "_")+'_'+tag+'_All'
  with open(fil2, 'a') as f:
    for row in rows:
      f.write(row)
    f.write('\n')
  return []

def cp_analytic(x,cs,dc,ts,ss,dH,sp,tp):
  # print(type(x))
  # print(len(x))
  # erfarg = x.apply(lambda i: (i-ts)/math.sqrt(2*(ss**2)))
  # exparg = x.apply(lambda i: -1*(i-tp)**2/(2*sp**2))
  erfarg = (x-ts)/math.sqrt(2*(ss**2))
  exparg = -1*(x-tp)**2/(2*sp**2)

  def sigm(erfarg):
    return erfarg.apply(lambda j: cs + (dc / 2) * (1 + scipy.special.erf(j)))
    # return cs + (dc / 2) * (1 + scipy.special.erf(erfarg))
  def gauss(exparg):
    return exparg.apply(lambda j: (dH / (sp * math.sqrt(2 * (math.pi)))) * math.exp(j))
    # return (dH / (sp * math.sqrt(2 * (math.pi)))) * math.exp(exparg)

  a = sigm(erfarg)
  b = gauss(exparg)

  # print(exparg)

  return a.add(b) #pd.Series.sum(sigm(erfarg),gauss(exparg)) #sigmoid function + gaussian
  # return pd.cs + (dc / 2) * (1 + scipy.special.erf(erfarg)) + (dH / (sp * math.sqrt(2 * (math.pi)))) * math.exp(exparg)
  #refer to supplementary of D. Danilov 2016 paper

def cp_fit(fil, val):
  fn = fil+'log.ng-dsc_'+ val +  '_' +fil.split("/")[-4]
  log = lmplg.File(fn)

  T = log.get("Temp")
  # T = numpy.linspace(-2*math.pi, 2*math.pi)
  U = log.get("TotEng")
  # U = numpy.sin(T)
  P = log.get("Press")*100000 #bars to Pascal conversion
  V = log.get("Volume")*1e-30 #A^3 to m^3   #get volume

  dV = numpy.diff(V,n=1)
  dV = numpy.delete(dV,0)
  T2 = numpy.delete(T,[0,1])
  U2 = numpy.delete(U,[0,1])
  P = numpy.delete(P,[0,1])
  V = numpy.delete(V,[0,1])

  # enth = U + numpy.multiply(P,dV)
  enth = U2 #+ numpy.multiply(P,dV)/(1.60218e-19) #Joules to eV

  # T2 = [i for i in range(1,10,1)]
  # # enth = [i*i for i in T2]
  # enth = [6*x+3 for x in T2]

  #enthalpy output to PNG
  fig, ax = plt.subplots(figsize=(5, 5))  # ,sharey=True)
  ax.set_facecolor('white')
  plt.setp(ax.spines.values(), linewidth=1.5)
  ax.grid(color='k', alpha=0.1, zorder=-2)
  ax.plot(T2,enth,marker=".",mec='k',label=val)
  ax.set_ylabel('Total Enthalpy (eV)',fontsize=14)
  ax.set_xlabel('Temperature (K)',fontsize=14)
  fig.tight_layout(rect=[0, 0.03, 1, 0.95])
  plt.legend()
  plt.savefig('enth_' + val + '.png', dpi=400)
  plt.close()

  #numerical derivative of enthalpy to get C_p
  # d_T    = numpy.gradient(T)
  # d_enth = numpy.gradient(enth)
  # C_p = d_enth/d_T
  C_p = numpy.diff(enth)/numpy.diff(T2)
  T3 = numpy.delete(T2,0)

  # print(len(C_p),len(T3),len(enth))
  # f = InterpolatedUnivariateSpline(T, enth, k=1)
  # dfdx = f.derivative()
  # dydx = dfdx(x)

  #initial estimates
  cs_est = C_p[0]
  dc_est = abs(C_p[0]-C_p[-1])
  ts_est = 800 #I am putting this by hand
  ss_est = 1000
  dH_est = 100
  sp_est = 1000
  tp_est = 800 #I am putting this by hand

  df = pd.DataFrame({'T':T3,'Cp':C_p,'enth':numpy.delete(enth,0)})

  with open('enthalpy_data-'+val+'.txt', 'w') as f:
    dfAsString = df.to_string(header=False, index=False)
    f.write(dfAsString)

  # df["model"] = cp_analytic(df["T"],cs_est,dc_est,ts_est,ss_est,dH_est,sp_est,tp_est)
  #
  # popt, pcov = curve_fit(f=cp_analytic,
  #                        xdata=df["T"],
  #                        ydata=df["Cp"],
  #                        p0=(cs_est,dc_est,ts_est,ss_est,dH_est,sp_est,tp_est)
  #                        )
  #
  # df["fit"] = cp_analytic(df["T"],popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6])

  # R2 = numpy.sum((df["fit"] - df["Cp"].mean())**2) / numpy.sum((df["Cp"] - df["Cp"].mean())**2)

  fig, ax = plt.subplots(figsize=(5, 5))  # ,sharey=True)
  plt.setp(ax.spines.values(), linewidth=1.5)
  ax.set_facecolor('white')
  ax.grid(color='k', alpha=0.1, zorder=-2)
  ax.plot(df["T"],df["Cp"],marker=".",mec='k',label=val)
  # ax.plot(df["T"],df["fit"],"r-",label='fit')
  ax.set_ylabel(r'Heat Capacity C$_p$',fontsize=14)
  ax.set_xlabel('Temperature (K)',fontsize=14)
  fig.tight_layout(rect=[0, 0.03, 1, 0.95])
  plt.legend()
  plt.savefig('cp_' + val + '.png', dpi=400)
  plt.close()