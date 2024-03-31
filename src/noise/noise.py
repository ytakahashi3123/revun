#!/usr/bin/env python3

import numpy as np
import vtk as vtk
from vtk.util import numpy_support
from vtk.numpy_interface import dataset_adapter
from modal.modal import modal


class noise(modal):

  def __init__(self):
    print("Calling class: noise")


  def get_variable_name(self, config):
    value_name = config['kind_flowvar'] 
    if config['kind_flowvar_replace'][0]:
      value_name = config['kind_flowvar_replace'][1]
    return value_name
    

  def routine_add_noise(self, config, time_dict):

    print("Adding noise to target variables")

    # Inital settings
    step_start    = time_dict['step_start']
    step_end      = time_dict['step_end']
    step_interval = time_dict['step_interval']
    num_step      = time_dict['num_step']
    step_digit    = config['step_digit']

    # Target flow variable to run POD
    kind_flowvar         = config['kind_flowvar']
    flag_flowvar_replace = config['kind_flowvar_replace'][0]
    if flag_flowvar_replace:
      varname_replaced_top = config['kind_flowvar_replace'][1]
      varname_replaced_bot = config['kind_flowvar_replace'][2]

    # Specific coordinate point
    coord_probe = config['coordinate_feature_noise']

    # Make file output directory
    self.make_directory( config['dir_result_noise'] )

    # Read object, add noise, and then, write objetc in VTK file formant
    for n in range(0, num_step):
      nstep = step_start + n*step_interval 
      # Reading data
      addfile      = '_'+str(nstep).zfill(step_digit)
      filename_tmp = self.split_file(config['dir_input'] +'/'+ config['file_input_base'],addfile,'.')
      filename_ext = self.get_extension(filename_tmp)
      print('--Reading VTK/VTU file: ', filename_tmp)
      reader       = self.get_vtk_reader(filename_tmp,filename_ext)

      # Extract variable
      flowvar_tmp = numpy_support.vtk_to_numpy( reader.GetOutput().GetPointData().GetAbstractArray(kind_flowvar) )
      if flag_flowvar_replace:
        var_replace_bot = numpy_support.vtk_to_numpy( reader.GetOutput().GetPointData().GetAbstractArray(varname_replaced_bot) )
        flowvar_tmp = flowvar_tmp/var_replace_bot[:,None]

      # Initial settings
      if n==0:
        print('Array names: ', self.get_arrayname(reader) )
        # --Set array shape
        if flowvar_tmp.ndim == 1:
          num_element  = flowvar_tmp.shape[0]
          num_component = 1
        else :
          num_element   = flowvar_tmp.shape[0]
          num_component = flowvar_tmp.shape[1]
        print('Number of elements   :',num_element)
        print('Number of components :',num_component)
        history_probe = np.zeros(num_step*num_component).reshape(num_step,num_component)

      # Add noise to the target variables
      noise_magnitude = config['magnitude_noise']
      if len(noise_magnitude) != flowvar_tmp.shape[1]:
        print("Dimension of noise_magnitude is different from that of the target variable, check config file.")
        exit()
      for i in range(0, len(noise_magnitude)):
        noise = np.random.normal(0, noise_magnitude[i], flowvar_tmp.shape[0]) # Make white noise
        flowvar_tmp[:,i] = flowvar_tmp[:,i] + noise[i]

      # Write mean data
      reader_tmp = self.delete_initialarry_variables_vtk(reader)
      value_name = self.get_variable_name(config)
      add_value  = dataset_adapter.WrapDataObject(reader_tmp.GetOutput())
      NewDataSet = flowvar_tmp
      add_value.PointData.append(NewDataSet, value_name)
      filename_tmp = config['dir_result_noise']+'/'+ self.split_file(config['file_result_noise'],addfile,'.')
      filename_ext = self.get_extension(filename_tmp)
      print('--Writing VTK/VTU file: ', filename_tmp)
      writer = self.get_vtk_writer(filename_tmp,filename_ext)
      writer.SetInputData(add_value.VTKObject)
      writer.Write()


      # Probe data extraction
      # Load the vtkPolyData object
      geometry_filter = vtk.vtkGeometryFilter()
      geometry_filter.SetInputConnection(reader.GetOutputPort()) 
      geometry_filter.Update()
      geometrydata = geometry_filter.GetOutput()
      # Set the dataset for the point locator and build it
      point_locator = vtk.vtkPointLocator()
      point_locator.SetDataSet(geometrydata)
      point_locator.BuildLocator()
      # Find the closest point
      closest_point_index = point_locator.FindClosestPoint( coord_probe )
      history_probe[n] = flowvar_tmp[closest_point_index]

    # Output history
    time = time_dict['time']
    time = time + config['time_offset_noise']
    filename_tmp = config['dir_result_noise']+'/'+ config['file_hisotry_noise']
    print('Writing history: ', filename_tmp)
    file_output = open( filename_tmp, 'w')
    header_tmp  = 'variables=Time[s]'
    for i in range(0, len(coord_probe)):
      header_tmp = header_tmp + ', Var_'+str(i+1)+'_'+str(coord_probe[i])
    header_tmp =  header_tmp + '\n'
    file_output.write( header_tmp )
    text_tmp = ''
    for n in range(0,num_step):
      text_tmp = text_tmp + str( time[n] )
      for i in range(0, len(coord_probe)):
        text_tmp = text_tmp + ', ' + str(history_probe[n][i])
      text_tmp = text_tmp + '\n'
    file_output.write( text_tmp )
    file_output.close()

    return
