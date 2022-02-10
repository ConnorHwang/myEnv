import glob
import numpy as np

filename=glob.glob("/p/lustre2/hwang8/helmholts/dfd_case1_highWe/data/*.vtu")
filename.sort()

min = 0
max = len(filename)-1

savename = "case8"
savedir = "/p/lustre1/hwang8/ensight_stl/"

index_init = 102500
index_fin  = 170000
index_inc  = 500
range_start = 166000
range_final = 170000
range_start_idx = int((range_start-index_init)/index_inc)
range_final_idx = int((range_final-index_init)/index_inc)

for ii in range(range_start_idx,range_final_idx+1):

  if ii != range_start_idx:
    # open a new file
    ensight.case.link_modelparts_byname("OFF")
    ensight.solution_time.monitor_for_new_steps("off")
    ensight.case.apply_context("OFF")
    ensight.case.replace("Case 1","Case 1")
    ensight.case.select("Case 1")

  #ensight.command.delay_refresh("ON")
  ensight.part.select_default()
  ensight.part.modify_begin()
  ensight.part.elt_representation("3D_border_2D_full")
  ensight.part.modify_end()
  ensight.data.binary_files_are("native")
  ensight.data.format("VTK")
  ensight.data.reader_option("'Debug Mode' OFF")
  ensight.data.shift_time(1.000000,0.000000,0.000000)
  ensight.solution_time.monitor_for_new_steps("off")
  #ensight.data.replace("/p/lustre2/hwang8/helmholts/dfd_case1_highWe/data/isovof.00169000.vtu")
  ensight.data.replace(filename[ii])
  
  ensight.variables.activate("vof")
  ensight.variables.activate("vof")
  ensight.isos.select_default()
  ensight.part.modify_begin()
  ensight.isos.variable("vof")
  ensight.part.modify_end()
  ensight.isos.select_default()
  ensight.part.modify_begin()
  ensight.isos.component(0,0,0)
  ensight.part.modify_end()
  ensight.part.select_begin(1)
  ensight.isos.begin()
  ensight.isos.variable("vof")
  ensight.isos.end()
  ensight.isos.create()
  ensight.part.select_begin(1)
  ensight.part.visible("OFF")
  ensight.part.select_begin(2)
  ensight.savegeom.format("STL")
  ensight.savegeom.parameters("vof")
  ensight.savegeom.binary("ON")
  #ensight.savegeom.save_geometric_entities("/p/lustre2/hwang8/helmholts/dfd_case1_highWe/ensight_stl3")
  
  if ii < 10:
    filepath = savedir+savename+"."+"000"+str(ii)
  elif ii < 100:
    filepath = savedir+savename+"."+"00"+str(ii)
  elif ii < 1000:
    filepath = savedir+savename+"."+"0"+str(ii)
  else:
    filepath = savedir+savename+"."+str(ii)
  
  ensight.savegeom.save_geometric_entities(filepath)

  print('current index: ' + str(ii) + ', range: ' + str(range_start_idx) + '~' + str(range_final_idx))
