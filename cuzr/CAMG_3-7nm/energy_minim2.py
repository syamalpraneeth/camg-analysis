#!/usr/bin/env python3

import os, sys, csv, glob
from lammps import lammps
import ovito
from ovito.io import *
sys.path.append('..')
from analysis import *
import numpy as numpy

def minimize(name,src,f1,l1,f2,l2,f3,l3,f4,l4,f5,l5,f6,l6):
  file_list=[f1,f2,f3,f4,f5,f6]
  lab_list=[l1,l2,l3,l4,l5,l6]

  args = "-screen out.lammps"
  words = args.split()

  for f in glob.glob('tmp.pote-min_coll_' + name.replace(" ", "_") + '_asp_*'):
    if os.path.exists(f): os.remove(f)

  l = lammps(cmdargs=words)
  for file,lab in zip(file_list,lab_list):
  #  with open (src+file) as fh:
  #     next(fh)
    fn = file.split("/")[-1].split(".")[-1]
    fl1 = file[len(src):]
    fl2 = file.split("/")[-1].split(".")[-1]
    subloc = fl1.replace(fl2,'').split(".")[0][:-4]
    print(fl1)
    print(fl2)
    print(subloc)
    print(src+subloc+"restart."+fn)

    l.command("read_restart "+src+subloc+"restart."+fn)
    if "MG" not in lab:
      l.commands_string("fix corint all property/atom i_int \n"
                        "set group cores i_int 2 \n"
                        "set group interface i_int 3")
      if "NG" not in lab:
        l.commands_string("fix layrs all property/atom i_lyr \n"
                          "set group buffr i_lyr 1 \n "
                          "set group statd i_lyr 2 \n "
                          "set group fixed i_lyr 3 \n"
                          "group substrate union buffr statd fixed \n"
                          "set group substrate i_int 1")
    else: pass
    l.command("compute c2 all pe/atom")
    l.command("neigh_modify every 1 delay 0 check yes")
    l.command("pair_style eam/fs")
    l.command("pair_coeff * * Cu-Zr_4.eam.fs Cu Zr Cu Zr")
    l.command("thermo_style custom step pe ke press temp")
    l.command("thermo_modify format float '% .6e'")
    l.command("timestep 0.01")  # 10 times more than dynamics tstep according to LAMMPS manual
    l.command("min_style fire")
    l.command("minimize 1e-10  1e-10  100  100")
    l.command("fix fpe all ave/atom 1 1 1 c_c2")
    l.command("run 0 pre yes post no")
    l.command("variable pe_norm equal $(pe)/$(atoms)")
    l.command("print 'normalized total minimized PE ${pe_norm}'")
    pe_norm = l.extract_variable('pe_norm')
    if "MG" in lab: l.command("write_dump all custom dump.minim."+name.split(" ")[0]+"_"+lab+" id type x y z vx vy vz f_fpe")
    elif "NG" in lab: l.command("write_dump all custom dump.minim."+name.split(" ")[0]+"_"+lab+" id type x y z vx vy vz i_int f_fpe")
    else: l.command("write_dump all custom dump.minim."+name.split(" ")[0]+"_"+lab+" id type x y z vx vy vz i_int i_lyr f_fpe")
    l.command("clear")
    # l.close()
    pipeline = import_file("dump.minim."+name.split(" ")[0]+"_"+lab)
    atom_types = pipeline.source.data.particles_.particle_types_
    for x in atom_types:
      if x % 2 == 0: atom_types.type_by_id(x).radius = 1.55
      else: atom_types.type_by_id(x).radius = 1.35
    # print('ATOM TYPES',numpy.max(atom_types))
    data = pipeline.compute()
    pe_atom = data.particles['f_fpe']
    print(numpy.sum(pe_atom), len(pe_atom),'before deletion')
    if numpy.max(atom_types)>2:
      pipeline.modifiers.append(ExpressionSelectionModifier(expression='ParticleType < 3'))
      pipeline.modifiers.append(DeleteSelectedModifier())
      data = pipeline.compute()
    pipeline = rem_surf("dump.minim."+name.split(" ")[0]+"_"+lab,pipeline,lab)
    data = pipeline.compute()

    if "MG" in lab: cs_case=['All']
    else: cs_case= 'All Core Interface'.split()

    for cs in cs_case:
      if cs == 'Core' or cs == 'Interface':
        cs_mod = ExpressionSelectionModifier(expression='i_int!=' + str(cs_case.index(cs) + 1))
        pipeline.modifiers.append(cs_mod)
        delcs_mod = DeleteSelectedModifier()
        pipeline.modifiers.append(delcs_mod)
        data = pipeline.compute()

        pe_atom = data.particles['f_fpe']
        print(cs,numpy.sum(pe_atom),len(pe_atom))
        pe_norm = numpy.sum(pe_atom)/len(pe_atom)

        delcs_mod.enabled = False
        cs_mod.enabled = False

      else:
        pe_atom = data.particles['f_fpe']
        print(cs,numpy.sum(pe_atom),len(pe_atom))
        pe_norm = numpy.sum(pe_atom)/len(pe_atom)

      rows = zip([pe_norm], [lab])
      with open('tmp.pote-min_coll_' + name.replace(" ", "_") + '_asp_' + cs, "a+") as f:
        writer = csv.writer(f)
        for row in rows:
          writer.writerow(row)
