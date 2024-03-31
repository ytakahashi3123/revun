#!/usr/bin/env python3

import numpy as np
import vtk as vtk
from vtk.util import numpy_support
from vtk.numpy_interface import dataset_adapter
from modal.modal import modal


class filter_data(modal):

  def __init__(self):
    print("Calling class: filter_data")
    
  def check_extension(self, filename):
    flagname_ext  = self.get_extension(filename)
    if flagname_ext != '.png' :
      print('The file format is not supported and only the PNG file is available:', flagname_ext)
      print('Program stopped.')
      exit()
    return

  def check_filesize(self, reader):
    dimensions = reader.GetOutput().GetDimensions()
    width  = dimensions[0]
    height = dimensions[1]
    print('--Width and height:', width, height)
    return

  def filter_png(self, config, time_dict):

    print("Making an image file coarse...")

    # Inital settings
    step_start    = time_dict['step_start']
    step_end      = time_dict['step_end']
    step_interval = time_dict['step_interval']
    num_step      = time_dict['num_step']
    step_digit    = config['step_digit']

    # Check file format; Only the PNG file is allowed
    filename_input_base = config['dir_input'] +'/'+ config['file_input_base']
    self.check_extension(filename_input_base)
    
    filename_result_base = config['dir_result_filterimages'] + '/'+ config['file_result_filterimages']
    self.check_extension(filename_result_base)

    # Make file output directory
    self.make_directory( config['dir_result_filterimages'] )

    # Set number of interval
    num_interval = config['num_interval_filterimages']

    # Read object and run filtering process
    for n in range(0, num_step):
      nstep = step_start + n*step_interval
      # Reading data
      addfile      = '_'+str(nstep).zfill(step_digit)
      filename_tmp = self.split_file(filename_input_base,addfile,'.')
      flagname_ext = self.get_extension(filename_tmp)
      reader       = self.get_vtk_reader(filename_tmp,flagname_ext)
      print('--Reading file: ', filename_tmp)

      # 画像データの取得
      image_data = reader.GetOutput()

      # 画像データの次元を取得
      image_dims = image_data.GetDimensions()
      new_dims = (image_dims[0]//num_interval, image_dims[1]//num_interval)

      print('--Width and height (Original):', image_dims[0], image_dims[1])
      print('--Width and height (filtered):', new_dims[0], new_dims[1])

      # 新しい画像データを作成
      new_image_data = vtk.vtkImageData()
      new_image_data.SetDimensions(new_dims[0], new_dims[1], 1)  # 2次元画像なので3つ目の次元は1

      new_image_data.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 3)  # 3チャンネルのカラー画像

      # 新しい画像データにデータをコピーして設定
      for j in range(new_dims[1]):
        for i in range(new_dims[0]):
          for k in range(0,3):
            pixel = image_data.GetScalarComponentAsFloat(i*num_interval,j*num_interval,0,k)
            new_image_data.SetScalarComponentFromFloat(i,j,0,k,pixel)

      # 新しい画像データをPNGとして出力
      filename_tmp = self.split_file(filename_result_base,addfile,'.')
      flagname_ext = self.get_extension(filename_tmp)

      png_writer = vtk.vtkPNGWriter()
      png_writer.SetInputData(new_image_data)
      png_writer.SetFileName(filename_tmp)
      png_writer.Write()
      print('--Writing file: ', filename_tmp)

    return
