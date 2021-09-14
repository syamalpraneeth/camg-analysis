##!/usr/bin/env python3

import os, sys
import ovito
from ovito.io import *
from ovito import scene
from ovito.modifiers import *
from ovito.vis import *
import numpy
import math, csv
# from ovito.vis import Viewport, CoordinateTripodOverlay
from PySide2.QtCore import Qt
from PySide2 import QtCore
from sort_voro import *

def rem_surf(pipeline,case):
  if case!= 'NG':
    string = '(Position.Z>50 || Position.Z<25) || ((abs(Position.X)>70) || (abs(Position.Y)>70))'
    pipeline.modifiers.append(ExpressionSelectionModifier(expression=string))
    pipeline.modifiers.append(DeleteSelectedModifier())
    data = pipeline.compute()

  # iter = 0
  # while iter < 2:
  #   surf = ConstructSurfaceModifier(radius=3, smoothing_level=10, select_surface_particles=True)
  #   pipeline.modifiers.append(surf)
  #   # pipeline.modifiers.append(ExpressionSelectionModifier(expression='SurfaceDistance<3'))
  #   data = pipeline.compute()
  #   surf.vis.enabled = False
  #   pipeline.modifiers.append(DeleteSelectedModifier())
  #   iter += 1
  # #cluster analysis
  # pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=3.5, sort_by_size=True))
  # #choose largest cluster
  # pipeline.modifiers.append(ExpressionSelectionModifier(expression='Cluster!=1'))
  # pipeline.modifiers.append(DeleteSelectedModifier())
  return pipeline;

def ovito_voro(fil, val, cs, species, delsub):
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
  fn = fil.split("/")[-1].split(".")[1] + '_' + val

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
    slice_mod = SliceModifier()
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
    if delsub == 1:
      vp.camera_pos = (-0.770151, 0.773092, 83.9268)
      vp.camera_dir = (-0.999997, 0.00251444, 4.57373e-15)
      vp.fov = 120.403
    else:
      vp.camera_pos = (36.8437, 36.8438, 36.845)
      vp.camera_dir = (0, 1, 0)
      vp.fov = 83.3431
    vp.overlays.append(tripod)
    vp.render_image(size=(800, 600), filename=fn + "_ortho.png", background=(1, 1, 1))
    pipeline.remove_from_scene()

  if delsub == 1: slice_mod.enabled = False  # disable slice, which was used for the visual
  if delsub == 1 and cs==0 and val != 'NG' and '1e10' in fil.split("/"): ovito_strain(pipeline,fil,val)

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
    slice_mod2 = SliceModifier()
    pipeline.modifiers.append(slice_mod2)
    # slice_mod2.enabled = True  # disable slice, which was used for the visual

  # Viz Slab
  pipeline.add_to_scene()
  vp = Viewport()
  vp.type = Viewport.Type.Ortho
  vp.camera_pos = (-0.770151, 0.773092, 83.9268)
  vp.camera_dir = (-0.999997, 0.00251444, 4.57373e-15)
  vp.fov = 120.403
  vp.overlays.append(tripod)
  if cs == 0 and species==0:
    vp.render_image(size=(800, 600), filename=fn + '_' + case + "_ortho_del.png", background=(1, 1, 1))
  elif cs > 0 and species==0:
    vp.render_image(size=(800, 600), filename=fn + '_' + case + '_' + csstr + "_ortho_del.png", background=(1, 1, 1))
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
  vp.camera_pos = (0, 0, 0)
  vp.camera_dir = (0, 0, -1)
  vp.fov = 90
  if cs == 0 and species==0:
    vp.render_image(size=(800, 800), filename=fn +'_chains-il_'+ case + ".png", background=(0, 0, 0))
  elif cs > 0 and species==0:
    vp.render_image(size=(800, 800), filename=fn +'_chains-il_' + case + '_' + csstr + ".png", background=(0, 0, 0),
                    frame=8)
  else: pass
  pipeline.remove_from_scene()

  return [voro_indices, numpy.array(atom_volume), csiz_arr, tot_selec_atms,pipeline];


def voro(fil, val, cs, species, delsub):
  fn = fil.split("/")[-1].split(".")[1] + '_' + val
  siz = 8
  sz = siz - 1
  [a, vols, csa, N,pipeline] = ovito_voro(fil, val, cs, species, delsub)

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
  pipeline.add_to_scene()
  ovito.scene.save("ico-chains_"+fn+'_' + case + '_'+ csstr +".ovito")
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


def voro_data_bulk(fil, val, cs):
  d1, d2, d3 = datclass(), datclass(), datclass()
  d1.a, d1.b, d1.c, d1.d, d1.h, d1.v, d1.csa = voro(fil, val, cs, 0, 0)
  d2.a, d2.b, d2.c, d2.d, d2.h, d2.v, d2.csa = voro(fil, val, cs, 1, 0)
  d3.a, d3.b, d3.c, d3.d, d3.h, d3.v, d3.csa = voro(fil, val, cs, 2, 0)
  #    print('From inside voro_data_bulk: ',numpy.shape(d1.v),numpy.shape(d2.v),numpy.shape(d3.v))
  return [d1, d2, d3];


def voro_data_film(fil, val, cs):
  d1, d2, d3 = datclass(), datclass(), datclass()
  d1.a, d1.b, d1.c, d1.d, d1.h, d1.v, d1.csa = voro(fil, val, cs, 0, 1)
  d2.a, d2.b, d2.c, d2.d, d2.h, d2.v, d2.csa = voro(fil, val, cs, 3, 1)
  d3.a, d3.b, d3.c, d3.d, d3.h, d3.v, d3.csa = voro(fil, val, cs, 4, 1)
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
  fn = fil.split("/")[-1].split(".")[1] + '_' + case

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
      vp.camera_pos = (308.835, -411.782, 438.376)
      vp.camera_dir = (-0.4819, 0.6424, -0.5956)
      vp.fov = math.radians(35)
      vp.overlays.append(tripod)
      vp.render_image(size=(800, 600), filename=fn + "_pers.png", background=(1, 1, 1))
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
    vp.camera_pos = (308.835, -411.782, 438.376) #(247.653, -179.174, 375.894) #
    vp.camera_dir = (-0.4819, 0.6424, -0.5956) # (-0.62888, 0.4386, -0.6412) #
    vp.fov = math.radians(35)
    vp.overlays.append(tripod)
    if cs == 0:
      vp.render_image(size=(800, 600), filename=fn + "_pers_del.png", background=(1, 1, 1))
    elif cs > 0:
      vp.render_image(size=(800, 600), filename=fn + '_' + csstr + "_pers_del.png", background=(1, 1, 1))
    pipeline.remove_from_scene()

    atom_types = data.particles['Particle Type']
    print('Num slab atoms:', len(atom_types))

  atom_types = data.particles['Particle Type']
  print('Final num atoms:', len(atom_types))

  # Output PE/atom histogram data to tmp file
  modifier = HistogramModifier(bin_count=100, property='f_fpe')
  pipeline.modifiers.append(modifier)
  if cs == 0: histfil = 'tmp.histo_' + case
  if cs > 0: histfil = 'tmp.histo_' + case + '_' + csstr
  export_file(pipeline, histfil, "txt/table", key="histogram[f_fpe]")

  potarr = numpy.array(data.particles['f_fpe'])
  # print(potarr)
  pote= potarr.sum(axis=0)/len(atom_types)
  # rows=zip(str(pote),case)
  rows = [str(pote),',',case]
  # print(rows)
  fil = 'tmp.pote_coll_' + name.replace(" ", "_") + '_' + tag + '_' + csstr
  with open(fil, 'a') as f:
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
  vp.overlays.append(ol)
  vp.type = Viewport.Type.Top
  if delsub == 1:
    vp.camera_pos = (0,0,57) #(0, 0, 0)
  else:
    vp.camera_pos = (36.84, 36.84, 36.845)
  vp.camera_dir = (0, 0, -1)
  if delsub == 1:
    vp.fov = 123
  else:
    vp.fov = 83.34

  if cs == 0:
    vp.render_image(size=(800, 800), filename=fn + "_slice.png", background=(1, 1, 1))
  elif cs > 0:
    vp.render_image(size=(800, 800), filename=fn + "_slice_" + csstr + ".png", background=(1, 1, 1))
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
    vp.render_image(size=(800, 800), filename=fn + "_slice-color.png", background=(1, 1, 1))
  elif cs > 0:
    vp.render_image(size=(800, 800), filename=fn + "_slice-color_" + csstr + ".png", background=(1, 1, 1))

#  print("works")
  pipeline.remove_from_scene()

  return []


def ovito_rdf(fil, case, cs, delsub):
  if cs == 0:
    csstr = 'All'
  elif cs == 1:
    csstr = 'Core'
  elif cs == 2:
    csstr = 'Interface'
  else:
    raise ValueError('cs value only 0, 1 or 2!')
  fn = fil.split("/")[-1].split(".")[1] + '_' + case

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

  #Volume calculation using surface mesh
  remesh = ConstructSurfaceModifier(radius=3, smoothing_level=10, select_surface_particles=True,identify_regions=True)
  pipeline.modifiers.append(remesh)
  # pipeline.modifiers.append(ConstructSurfaceModifier(radius=3, smoothing_level=10, select_surface_particles=True,identify_regions=True))

  data = pipeline.compute()

  cell_volume = data.cell.volume
  region_volume = data.attributes["ConstructSurfaceMesh.filled_volume"]
  frac = data.attributes["ConstructSurfaceMesh.filled_fraction"]
  print(cell_volume,region_volume,frac)

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
  rdffil = 'tmp.rdf_' + case +'_'+csstr #'_'+str(typ)
  # print(data.tables['coordination-rdf'].y.component_names)
  # numpy.savetxt(rdffil, data.tables['coordination-rdf'].xy())
  numpy.savetxt(rdffil, rdf_data)

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


def reassemble(h1, b1, h2, b2, h3, b3, h4, b4):  # rearrange all b entries
  un1 = union(b1, b2)
  un2 = union(b3, un1)
  un3 = union(b4, un2)
  in1 = inter(b1, b2)
  in2 = inter(b3, in1)
  in3 = inter(b4, in2)

  ncdf = [x for x in un3 if x not in in3]

  b = in3 + ncdf

  H1, H2, H3, H4, H5, H6 = [], [], [], [], [], []
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

  return [H1, b, H2, b, H3, b, H4, b];


def collate_ico(name, o1, o2, o3, o4, p1, p2, p3, p4, l1, l2, l3, l4, cs, typ):
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
  a1 = [o1.a[o1.b.index(i_ico)], o2.a[o2.b.index(i_ico)], o3.a[o3.b.index(i_ico)], o4.a[o4.b.index(i_ico)]]
  a2 = [p1.a[o1.b.index(i_ico)], p2.a[o2.b.index(i_ico)], p3.a[o3.b.index(i_ico)], p4.a[p4.b.index(i_ico)]]
  a3 = [l1, l2, l3, l4]

  rows = zip(a1, a2, a3)
  fil = 'tmp.ico_coll_' + name.replace(" ", "_") + '_' + tag + '_' + csstr
  with open(fil, 'w') as f:
    writer = csv.writer(f)
    for row in rows:
      writer.writerow(row)

  b1 = [o1.h[0], o2.h[0], o3.h[0], o4.h[0]]
  b2 = [p1.h[0], o2.h[0], p3.h[0], p4.h[0]]
  rows2 = zip(b1, b2, a3)
  fil = 'tmp.ico-like_coll_' + name.replace(" ", "_") + '_' + tag + '_' + csstr
  with open(fil, 'w') as f:
    writer = csv.writer(f)
    for row in rows2:
      writer.writerow(row)

  c1 = [sum(o1.h[1:]), sum(o2.h[1:]), sum(o3.h[1:]), sum(o4.h[1:])]
  c2 = [sum(p1.h[1:]), sum(p2.h[1:]), sum(p3.h[1:]), sum(p4.h[1:])]
  rows3 = zip(c1, c2, a3)
  fil = 'tmp.nico_coll_' + name.replace(" ", "_") + '_' + tag + '_' + csstr
  with open(fil, 'w') as f:
    writer = csv.writer(f)
    for row in rows3:
      writer.writerow(row)


def ovito_strain(pipeline, fil, case):
  fn = fil.split("/")[-1].split(".")[1] + '_' + case

  # pipeline = import_file(fil)
  # # Set atomic radii (required for polydisperse Voronoi tessellation).
  # atom_types = pipeline.source.data.particles_.particle_types_
  # for x in atom_types:
  #   if x % 2 == 0:
  #     atom_types.type_by_id(x).radius = 1.55
  #   else:
  #     atom_types.type_by_id(x).radius = 1.35
  #
  # data = pipeline.compute()
  # at_typ_ind = data.particles['Particle Type']
  # print(fn, 'Strain:', csstr)
  #
  # if delsub == 1:
  #   pipeline.modifiers.append(ExpressionSelectionModifier(expression='ParticleType < 3'))
  #   data = pipeline.compute()
  #   pipeline.modifiers.append(DeleteSelectedModifier())  # Delete all substrate atoms
  #   data = pipeline.compute()

  file='/home/mj0054/Documents/work/simulations/projects/cibd/cibd_multi/50-50/50K/1e10/3nm/undeposited/data/dump.dep.cibd_3nm_6000meV-lay3-un'
  pip2 = import_file(file)
  strmod = AtomicStrainModifier(cutoff=3.8)

  strmod.reference = pip2.source
  strmod.reference.load(file)
  pipeline.modifiers.append(strmod)

  slice_mod2 = SliceModifier()
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
  vp.camera_pos = (-0.770151, 0.773092, 83.9268)
  vp.camera_dir = (-0.999997, 0.00251444, 4.57373e-15)
  vp.fov = 120.403
  vp.overlays.append(tripod)
  overlay = ColorLegendOverlay(
    modifier = color_mod2,
    title = 'Strain',
    alignment = Qt.AlignRight ^ Qt.AlignBottom,
    orientation = Qt.Vertical,
    offset_y = 0.06,
    font_size = 0.12,
    format_string = '%1.0f')
  vp.overlays.append(overlay)
  vp.render_image(size=(800, 600), filename=fn + "_strain_ortho.png", background=(1, 1, 1)) #alpha=True)
  pipeline.remove_from_scene()
  color_mod2.enabled = False
  slice_mod2.enabled = False
  return []
