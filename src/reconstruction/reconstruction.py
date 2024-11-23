#!/usr/bin/env python3

import numpy as np
import vtk as vtk
from vtk.util import numpy_support
from vtk.numpy_interface import dataset_adapter
from modal.modal import modal


class reconstruction(modal):

  def __init__(self):
    print("Calling class: reconstruction")

  def get_variable_name(self, config):
    value_name = config['kind_flowvar'] 
    if config['kind_flowvar_replace'][0]:
      value_name = config['kind_flowvar_replace'][1]
    return value_name

  def get_kind_reconst(self, config, num_rank, value_name):
    kind_reconst=[]
    if config['flag_allmodes_reconstruct']:
      for n in range(0, num_rank):
        kind_reconst.append( value_name + '_R' + str(n+1))
    else:
      mode_dmd_reconstruct = config['mode_dmd_reconstruct']
      for n in range(0, len(mode_dmd_reconstruct)):
        kind_reconst.append( value_name + '_' + str(mode_dmd_reconstruct[n][0]))
    return kind_reconst

  def get_mode_list(self, config, num_rank):
    mode_dmd_reconstruct = config['mode_dmd_reconstruct']
    mode_s_list=[]
    mode_e_list=[]
    if config['flag_allmodes_reconstruct']:
      #--All modes reconsrtucted
      print('All modes are reconstructed')
      for n in range(0, num_rank):
        mode_s_list.append(n)
        mode_e_list.append(n+1)
    else:
      #--Selective modes reconsrtucted
      print('Selected mode is reconstructed')
      for n in range(0,len(mode_dmd_reconstruct)):
        mode_start = mode_dmd_reconstruct[n][1]
        mode_end   = mode_dmd_reconstruct[n][2]
        # Modified settings
        if mode_start > mode_end:
          mode_end_tmp = mode_end
          mode_end     = mode_start
          mode_start   = mode_end_tmp
        #if mode_start > num_rank:
        #  mode_start = num_rank
        #if mode_end > num_rank:
        #  mode_end = num_rank
        #if mode_start == mode_end:
        #  mode_end = mode_start + 1
        mode_s_list.append(mode_start-1)
        mode_e_list.append(mode_end)
        if config['flag_verbose']:
          print("Reconstruction mode name:",mode_dmd_reconstruct[n][0], ", Slice:",mode_s_list[n],'--',mode_e_list[n])
    return mode_s_list, mode_e_list

  def write_reconstructed_vtk(self, num_component, kind_reconst, filename_list, reader_init, var_list):
    reader_tmp = self.delete_initialarry_variables_vtk(reader_init)
    for n in range(0,len(kind_reconst)):
# -------------------
      add_value  = dataset_adapter.WrapDataObject(reader_tmp.GetOutput())
      value_name = kind_reconst[n]
      var_tmp = var_list[n]
      NewDataSet = np.real( var_tmp[:,0:num_component] )
      add_value.PointData.append(NewDataSet, value_name)
# The procedure below is slow
#          value_name = kind_reconst[n]
#          add_value = vtk.vtkDoubleArray()
#          add_value.SetName(value_name)
#          add_value.SetNumberOfComponents(num_component)
#          add_value.SetNumberOfTuples(num_element)
##          for i in range(num_element):
##            var_tmp = var_list[n]
##            add_value.SetValue(i, np.real( var_tmp[i,0:num_component,m] ))
#          var_tmp = var_list[n]
#          add_value.SetVoidArray(np.real(var_tmp[:,0:num_component,m]).ravel(), num_element*num_component, 1) # ＜−ChatGPT
#          reader_tmp.GetOutput().GetPointData().AddArray(add_value)
# -------------------
      # --Write VTK
      filename_tmp = filename_list[n]
      filename_ext = self.get_extension(filename_tmp)
      writer = self.get_vtk_writer(filename_tmp,filename_ext)
# -------------------
      writer.SetInputData(add_value.VTKObject)
#        writer.SetInputData(reader_tmp.GetOutput())
# -------------------
      writer.Write()
    return

  def write_reconstructed_image(self, num_component, kind_reconst, filename_list, reader_init, var_list):
    from PIL import Image
    # Get width and height
    dimensions = reader_init.GetOutput().GetDimensions()
    width  = dimensions[0]
    height = dimensions[1]
    # Write image file
    for n in range(0,len(kind_reconst)):
      filename_png = filename_list[n]
      var_tmp = np.real( var_list[n] ).astype(np.uint8)
      var_tmp = var_tmp.reshape(height, width, num_component)
      image = Image.fromarray(var_tmp)
      # rotation_angle=90
      #image = image.rotate(rotation_angle, expand=True)
      image.save(filename_png)
    return

  def reconstruct_DMD_mode(self, config, time_dict, num_element, num_component, reader_init, spatial_mean, spatial_mode, time_evolution):

    print('Reconstruct flow field based on DMD modes')

    # Make file output directory
    self.make_directory( config['dir_reconstruct_dmd'] )

    # Inital settings
    step_start    = time_dict['step_start']
    step_end      = time_dict['step_end']
    step_interval = time_dict['step_interval']
    num_step      = time_dict['num_step']
    num_rank      = self.get_number_rank(config, time_dict)
    step_digit    = config['step_digit']

    # Reconstructed variables name
    value_name    = self.get_variable_name(config)
    kind_reconst = self.get_kind_reconst(config, num_rank, value_name)
    num_reconst = len(kind_reconst)

    # Set Mode ID
    mode_s_list, mode_e_list = self.get_mode_list(config, num_rank)

    #dimensions = reader_init.GetOutput().GetDimensions()
    #width = dimensions[0]
    #height = dimensions[1]

    # Reconstruction from spatial mode and time evoluation data
#    filename_base = config['dir_reconstruct_dmd']+'/'+ self.split_file(config['file_dmd_reconstruct'],'_'+config['kind_flowvar'],'.')
    filename_base = config['dir_reconstruct_dmd']+'/'+ config['file_dmd_reconstruct']
    file_name_ext = self.get_extension(filename_base)
    for n in range(0,num_step):
      print('--Writing reconstructed VTK file (base file name): ', n, filename_base)
      var_list = []
      filename_list = []
      for m in range(0,num_reconst):
        # File name
        filename_tmp = self.split_file(filename_base,'_'+kind_reconst[m],'.')
        filename_tmp = self.split_file(filename_tmp,'_'+str(n).zfill(step_digit),'.')
        filename_list.append(filename_tmp)
        # Reconstruct a flow field in space and time from mode data
        mode_start = mode_s_list[m]
        mode_end   = mode_e_list[m]
        var_tmp = np.dot(spatial_mode[:,:,mode_start:mode_end], time_evolution[mode_start:mode_end,n]) + config['flag_dmd_reconstruct_add_mean'] * spatial_mean[:,:]
        var_list.append(var_tmp)

      # Write data
      if file_name_ext == '.png' :
        self.write_reconstructed_image(num_component, kind_reconst, filename_list, reader_init, var_list)
      else:
        self.write_reconstructed_vtk(num_component, kind_reconst, filename_list, reader_init, var_list )

    return