#!/usr/bin/env python3

import numpy as np
import vtk as vtk
from vtk.util import numpy_support
from modal.modal import modal

class feature(modal):

  def __init__(self):
    print("Calling class: feature")


  def check_reconstruction(self, config):
    if not config['flag_dmd_reconstruct']:
      print('Flag of reconstruction is not valid.', config['flag_dmd_reconstruct'])
      print('Check flag_dmd_reconstruct in config. Program stopped.')
      exit()
    return

  def get_variable_name(self, config):
    value_name = config['kind_flowvar'] 
    if config['kind_flowvar_replace'][0]:
      value_name = config['kind_flowvar_replace'][1]
    return value_name

  def get_kind_reconst(self, config, num_rank, value_name):
    kind_reconst=[]
    #if config['flag_allmodes_reconstruct']:
    #  for n in range(0, num_rank):
    #    kind_reconst.append( value_name + '_R' + str(n+1))
    #else:
    #  mode_dmd_reconstruct = config['mode_dmd_reconstruct']
    #  for n in range(0, len(mode_dmd_reconstruct)):
    #    kind_reconst.append( value_name + '_' + str(mode_dmd_reconstruct[n][0]))
    for n in range(0, num_rank):
      kind_reconst.append( value_name + '_R' + str(n+1))
    return kind_reconst


  def get_closest_point(self, config, reader, coord_ref):
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
    closest_point_index = point_locator.FindClosestPoint( coord_ref )
    # Get the coordinates of the closest point
    points = geometrydata.GetPoints()
    if config['flag_verbose']:
      closest_point = points.GetPoint(closest_point_index)
      print('Closest point coordinate and index:',closest_point, closest_point_index)
    # Get the point data
    point_data = geometrydata.GetPointData()
    return closest_point_index, point_data


  def get_probe_data(self, var_name, closest_point_index, point_data ):
    # Get the vtkAbstractArray
    abstract_array = point_data.GetAbstractArray(var_name)
    # Cast vtkAbstractArray to vtkDataArray
    data_array = vtk.vtkDataArray.SafeDownCast(abstract_array)
    # Get the probe value at the closest point
    probe_data = np.array( data_array.GetTuple(closest_point_index) )
    return probe_data


  def extract_feature_DMD(self, config, time_dict, num_component, reader, spatial_mean, spatial_mode, time_evolution):

    print('Extract features from DMD spatiotemporal data')

    # Make file output directory
    self.make_directory( config['dir_feature'] )

    # Inital settings
    num_step   = time_dict['num_step']
    num_rank   = self.get_number_rank(config, time_dict)

    # Reconstructed variables name
    value_name = self.get_variable_name(config)
    kind_reconst = self.get_kind_reconst(config, num_rank, value_name)
    num_reconst = len(kind_reconst)

    # Case: target DMD variables:

    # Specific coordinate point
    # Connectivityは維持されるであろうという前提で、Initial or Sample　VTKからClosest pointを同定する
    # 問題があれば後日修正したい, 2023/06/30
    coord_ref = config['coordinate_feature']
    closest_point_index, point_data = self.get_closest_point(config, reader, coord_ref)

    # Number of array in feature extracting. Normally, this is the same as the DMD target variable's array (for example, n=3 for velocity), that is, num_component.
    #num_feature = spatial_mean.shape[1]
    num_feature = num_component

    # Extracting feature valus for DMD data
    # VTKの読み込みオブジェクト(reader)におけるインデックス構造と、それから読みこんでnumpy_support.vtk_to_numpy処理したあとの変数のインデクスは同じなのか？
    feature_var = np.zeros(num_step*num_reconst*num_feature).reshape(num_step,num_reconst,num_feature)
    feature_var_total = 0.0
    for n in range(0,num_step):
      for m in range(0, num_reconst):
        var_tmp = np.dot(spatial_mode[:,:,m:m+1], time_evolution[m:m+1,n]) + config['flag_dmd_reconstruct_add_mean'] * spatial_mean[:,:]
        feature_var[n,m,:] = np.real( var_tmp[closest_point_index,:] )
      feature_var_total = np.sum(feature_var[n,:,:], axis=0)
      print("Feature value", n, feature_var_total )

    feature_var_mean = np.zeros(num_step*num_feature).reshape(num_step,num_feature)
    if config['flag_subtraction_mean']:
      feature_var_mean[:,0:num_feature+1] = spatial_mean[closest_point_index,0:num_feature+1]
      print("Feature value (time-mean)", feature_var_mean[0])

    return num_feature, feature_var, feature_var_mean


  def get_aerdynamic_force_hisotry(self, config, time_dict):

    print('Get aerodynamic forces history from reconstruction data')

    # Make file output directory
    self.make_directory( config['dir_feature'] )

    # Inital settings
    num_step   = time_dict['num_step']
    num_rank   = self.get_number_rank(config, time_dict)
    step_digit = config['step_digit']

    # Reconstructed variables name
    value_name = self.get_variable_name(config)
    kind_reconst = self.get_kind_reconst(config, num_rank, value_name)
    num_reconst = len(kind_reconst)

    # Number of array in feature extracting
    num_feature = 6

    # Non-dimension parameters (output: aerodynamic forces or aerodynamic coefficients)
    nondimension_force  = 1.0
    nondimension_moment = 1.0
    if config['flag_aerodynamic_coefficient']:
      nondimension_force  = 0.5*config['aerodynamic_reference_density']*config['aerodynamic_reference_velocity']**2*config['aerodynamic_reference_area']
      nondimension_moment = nondimension_force*config['aerodynamic_reference_length']
    
    # Read object and calculate aerodynamic forces
    aerodynamic_force = np.zeros(num_step*num_reconst*num_feature).reshape(num_reconst,num_step,num_feature)
    filename_base = config['dir_reconstruct_dmd']+'/'+ config['file_dmd_reconstruct']
    filename_ext = self.get_extension(filename_base)
    for n in range(0,num_step):
      for m in range(0,num_reconst):
        filename_tmp = self.split_file(filename_base,'_'+kind_reconst[m],'.')
        filename_tmp = self.split_file(filename_tmp,'_'+str(n).zfill(step_digit),'.')
        reader       = self.get_vtk_reader(filename_tmp,filename_ext)
      
        force_list = self.calculate_aerodynamic_force(config, kind_reconst[m], reader)
        aerodynamic_force[m,n,0:3] = force_list[0:3]/nondimension_force
        aerodynamic_force[m,n,3:6] = force_list[3:6]/nondimension_moment

      print("aerodynamic force value", n, aerodynamic_force[0,n,:])

    aerodynamic_force_mean = np.zeros(num_step*num_feature).reshape(num_step,num_feature)
    if config['flag_subtraction_mean']:
      filename_tmp  = config['dir_result_mean']+'/'+ self.split_file(config['file_result_mean'],'_'+config['kind_flowvar'],'.')
      filename_ext  = self.get_extension(filename_tmp)
      reader        = self.get_vtk_reader(filename_tmp,filename_ext)

      force_list = self.calculate_aerodynamic_force(config, value_name, reader)
      aerodynamic_force_mean[:,0:3] = force_list[0:3]/nondimension_force
      aerodynamic_force_mean[:,3:6] = force_list[3:6]/nondimension_moment

      print("aerodynamic force value (time-mean)", aerodynamic_force_mean[0,:])

    return num_feature, aerodynamic_force, aerodynamic_force_mean


  def calculate_aerodynamic_force(self, config, kind_reconst, reader):
    # Note: This function is written with chatGPT's help...
    #
    # vtkUnstructuredGridReaderで読んだデータは一度vtkGeometryFilterでポリゴンデータ（三角形？）に変換する必要がある
    # （vtkUnstructuredGridReaderで読んだデータはPoint属性であり、Cell属性を有さない（SU2の場合に限る？））
    # その上でセル法線ベクトル、面積、流れ場のポイントデータをセルデータに変換し、各セルの空気力を計算・積分する必要がある。
    
    # Load the vtkPolyData object
    geometry_filter = vtk.vtkGeometryFilter()
    geometry_filter.SetInputConnection(reader.GetOutputPort()) 
    geometry_filter.Update()
    geometrydata = geometry_filter.GetOutput()

    # Get normal vector of cell surface
    cellnormal_filter = vtk.vtkPolyDataNormals()
    cellnormal_filter.SetInputData(geometrydata)
    cellnormal_filter.ComputeCellNormalsOn()
    cellnormal_filter.FlipNormalsOn()
    cellnormal_filter.Update()
    num_cell = cellnormal_filter.GetOutput().GetCellData().GetNormals().GetNumberOfTuples()
    cellnormal = numpy_support.vtk_to_numpy( cellnormal_filter.GetOutput().GetCellData().GetArray('Normals') )

    # Get area of cell surface
    area_filter = vtk.vtkCellSizeFilter()
    area_filter.SetInputData(geometrydata)
    area_filter.ComputeVertexCountOff()
    area_filter.ComputeLengthOff()
    area_filter.ComputeAreaOn()
    area_filter.Update()
    area = numpy_support.vtk_to_numpy( area_filter.GetOutput().GetCellData().GetArray("Area") )

    # Get cell center coordinate
    cellcenter_filter = vtk.vtkCellCenters()
    cellcenter_filter.SetInputData(geometrydata)
    cellcenter_filter.Update()
    points = cellcenter_filter.GetOutput().GetPoints()
    cellcenter = numpy_support.vtk_to_numpy( cellcenter_filter.GetOutput().GetPoints().GetData() )

    # Get variable oon cell surface by converting Point to Cell
    pointtocell_filter = vtk.vtkPointDataToCellData()
    pointtocell_filter.SetInputData(geometrydata)
    pointtocell_filter.Update()
    num_var = pointtocell_filter.GetOutput().GetCellData().GetNumberOfArrays()
    #for i in range(0, num_var):
    #  array = pointtocell_filter.GetOutput().GetCellData().GetArray(i)
    #  print(array.GetName())

    coord_ref = config['coordinate_feature']

    #aerodynamics_list = []
    force_tmp = np.zeros(num_cell*6).reshape(num_cell,6)
    aerosynamics_total = np.zeros(6)
    #for n in range(0, len(kind_reconst)):
    pressure = numpy_support.vtk_to_numpy( pointtocell_filter.GetOutput().GetCellData().GetArray(kind_reconst) )
    force_tmp[:,0:3] = (pressure*area).reshape(-1,1) * cellnormal
    force_tmp[:,3]  = (cellcenter[:,1]-coord_ref[1])*force_tmp[:,2] - (cellcenter[:,2]-coord_ref[2])*force_tmp[:,1]
    force_tmp[:,4]  = (cellcenter[:,2]-coord_ref[2])*force_tmp[:,0] - (cellcenter[:,0]-coord_ref[0])*force_tmp[:,2]
    force_tmp[:,5]  = (cellcenter[:,0]-coord_ref[0])*force_tmp[:,1] - (cellcenter[:,1]-coord_ref[1])*force_tmp[:,0]
    aerosynamics_total = np.sum(force_tmp, axis=0)

    aerodynamics_list = aerosynamics_total

    #aerodynamics_list.append(aerosynamics_total)

    return aerodynamics_list


  def write_feature_history(self, config, time_dict, num_feature, feature, feature_mean):
    # Inital settings
    num_step = time_dict['num_step']
    time     = time_dict['time']
    num_rank = self.get_number_rank(config, time_dict)
    delimiter_tmp = ','
    coord_ref = config['coordinate_feature']
    
    # Set kind of mode reconsructed
    value_name   = self.get_variable_name(config)
    kind_reconst = self.get_kind_reconst(config, num_rank, value_name)
    num_reconst  = len(kind_reconst)

    # Output aerodynamic data
    filename_tmp = config['dir_feature']+'/'+self.split_file(config['file_dmd_feature'],'_'+value_name,'.')
    print('Writing features:',filename_tmp)
    file_output = open( filename_tmp, 'w')
    header_tmp = 'variables=Time[s]'
    for k in range(0,num_feature):
      header_tmp = header_tmp + ', Feature'+str(k+1)
    header_tmp = header_tmp +'\n'
    file_output.write( header_tmp )
    for n in range(0,num_reconst):
      text_tmp='zone T="Reconstruction_'+str(kind_reconst[n])+'", i='+str(num_step)+' f=point'+'\n'
      for m in range(0,num_step):
        text_force = ''
        for k in range(0,num_feature):
          text_force = text_force + str( feature[m,n,k] ) + delimiter_tmp
        text_tmp = text_tmp + str(time[m]) + delimiter_tmp + text_force.rstrip(delimiter_tmp) + '\n'
      file_output.write( text_tmp )
    file_output.close()

    # TIme mean value
    if config['flag_subtraction_mean']:
      filename_tmp = config['dir_feature']+'/'+self.split_file(config['file_dmd_feature'],'_'+value_name+'_Mean','.')
      print('Writing time-mean features:',filename_tmp)
      file_output = open( filename_tmp, 'w')
      file_output.write( header_tmp )
      text_tmp='zone T="Time_mean", i='+str(num_step)+' f=point'+'\n'
      for m in range(0,num_step):
        text_force = ''
        for k in range(0,num_feature):
          text_force = text_force + str(feature_mean[m,k] ) + delimiter_tmp
        text_tmp = text_tmp + str(time[m]) + delimiter_tmp + text_force.rstrip(delimiter_tmp) + '\n'
      file_output.write( text_tmp )
      file_output.close()

    return

  def read_feature_history(self, config, time_dict, num_component):
    # Inital settings
    num_step = time_dict['num_step']
    num_rank = self.get_number_rank(config, time_dict)

    # Set kind of mode reconsructed
    value_name   = self.get_variable_name(config)
    kind_reconst = self.get_kind_reconst(config, num_rank, value_name)
    num_reconst  = len(kind_reconst)

    # Number of array in feature extracting. Normally, this is the same as the DMD target variable's array (for example, n=3 for velocity), that is, num_component.
    num_feature = num_component

    # Feature value hisotry array
    feature = np.zeros(num_reconst*num_step*num_feature).reshape(num_reconst,num_step,num_feature)

    # File name
    filename_tmp = config['dir_feature']+'/'+self.split_file(config['file_dmd_feature'],'_'+value_name,'.')
    print('Reading features history:',filename_tmp)
    file_input = open( filename_tmp, 'r')
    buffer = file_input.readlines()
    n_count = 1
    for n in range(0,num_reconst):
      for m in range(0,num_step):
        n_count=n_count+1
        buffr_tmp = buffer[n_count].split(',')
        feature[n,m,0:num_feature] = np.array(buffr_tmp[1:num_feature+1]).astype(float)
      n_count=n_count+1
    file_input.close()

    # Time mean value
    feature_mean = np.zeros(num_step*num_feature).reshape(num_step,num_feature)
    if config['flag_subtraction_mean']:
      filename_tmp = config['dir_feature']+'/'+self.split_file(config['file_dmd_feature'],'_'+value_name+'_Mean','.')
      print('Reading time-mean features:',filename_tmp)
      file_input = open( filename_tmp, 'r')
      buffer = file_input.readlines()
      n_count = 1
      for m in range(0,num_step):
        n_count = n_count+1
        buffr_tmp = buffer[n_count].split(',')
        feature_mean[m,0:num_feature] = np.array(buffr_tmp[1:num_feature+1]).astype(float)
      file_input.close()

    return num_feature, feature, feature_mean


  def mode_sensing_greedy(self, config, time_dict, num_feature, displacement, displacement_mean, frequency, amplitude):
    # 貪欲法を用いて再構築変位から重要なモードを抽出する。
    # 最初に最も貢献の大きなモードを抽出してプライマリーモードを定義、つぎにプライマリーモードに1つ1つ各モード構築値を足して貢献度を計測、最も成績の良いものを決定する、これを繰り返す。
    # 計算量はO(n^2) (or O(nlogn)?)となるはずだが、そもそもnは特異値分解が成立するくらいの量なので問題はあまりないと思う。

    print('Identifying the important modes based on reconstructed displacement by greedy algorithm')

    # Inital settings
    num_step = time_dict['num_step']
    num_rank = self.get_number_rank(config, time_dict)
    value_name   = self.get_variable_name(config)
    kind_reconst = self.get_kind_reconst(config, num_rank, value_name)
    num_reconst  = len(kind_reconst)

    # Search pair mode
    mode_pair = np.zeros(num_rank).reshape(num_rank).astype(int)
    for m in range(0,num_rank):
      mode_pair[m] = m
      idx = np.abs(np.asarray(frequency) - (-1.0*frequency[mode_pair[m]]) ).argmin()
      if abs(amplitude[m]- amplitude[idx])  <= 1.e-10:
        mode_pair[m] = idx

    # Modify reconstruction data according to the time-mean value subtraction
    # --Case for time mean subtraction
    if config['flag_subtraction_mean'] :
      # --Case for reconstuction data including time-mean data
      if config['flag_dmd_reconstruct_add_mean']: 
        displacement = displacement - displacement_mean[np.newaxis,:,:]

    # All modes considered
    displacement_total = np.sum(displacement,axis=0)

    # Identify the important modes using greedy algorithm
    num_rank_consider = config['num_considered_feature']
    displacement_reconst_accum = np.zeros(num_rank_consider*num_step*num_feature).reshape(num_rank_consider,num_step,num_feature)
    displacement_reconst       = np.zeros(num_step*num_feature).reshape(num_step,num_feature)
    contribution = np.zeros(num_reconst).reshape(num_reconst)
#    displacement_list = [[] for _ in range(num_reconst)]
    for i in range(0,num_feature):
      print('-'*100)
      print('Spatial index (for example, x):', str(i+1))
      flag_inclusion = [True]*num_reconst
      for n in range(0, num_rank_consider):
        #contribution =np.zeros(num_reconst)
        contribution.fill(1e50)
        for m in range(0,num_reconst):
          if flag_inclusion[m] :
            m_pair = mode_pair[m]
            # ユークリッド距離 or 相関係数で類似度を判断
            X = displacement_total[:,i]
            if m != m_pair :
              Y = displacement[m,:,i] + displacement[m_pair,:,i] + displacement_reconst[:,i]
            else :
              Y = displacement[m,:,i] + displacement_reconst[:,i]
            #contribution[m] = np.corrcoef(X, Y)[0, 1]
            contribution[m] = np.linalg.norm(X - Y)/(np.linalg.norm(X)+1.e-50)
    
        # Get the index which contribution value accumulated is hihgest
        max_index      = np.argmin(contribution, axis=0)
        max_index_pair = mode_pair[max_index]
        if max_index != max_index_pair :
          displacement_reconst[:,i] = displacement_reconst[:,i] + displacement[max_index:max_index+1,:,i] + displacement[max_index_pair:max_index_pair+1,:,i]
        else : 
          displacement_reconst[:,i] = displacement_reconst[:,i] + displacement[max_index:max_index+1,:,i]
        flag_inclusion[max_index] = False
        flag_inclusion[max_index_pair] = False
#        displacement_list[n].append(displacement_reconst[:,i])
        displacement_reconst_accum[n,:,i] = displacement_reconst[:,i]
        print('| Rank:',n+1,
              '| Mode IDs:',max_index+1, max_index_pair+1,
              '| Freq.','{:.5e}'.format(frequency[max_index]), 
              '| Amp.', '{:.5e}'.format(amplitude[max_index]), 
              '| Cont.','{:.5e}'.format(1.0-contribution[max_index])
              )
    return displacement_reconst_accum


  def write_feature_history_reconstructed(self, config, time_dict, num_feature, displacement, displacement_mean, displacement_reconst_accum):

    # Inital settings
    num_step = time_dict['num_step']
    num_rank = self.get_number_rank(config, time_dict)
    value_name   = self.get_variable_name(config)
    kind_reconst = self.get_kind_reconst(config, num_rank, value_name)
    num_reconst  = len(kind_reconst)
    num_rank_consider = config['num_considered_feature']

    # Output data to file
    time = time_dict['time'] + config['time_offset_feature']
    filename_tmp = config['dir_feature']+'/'+config['file_sensing_feature']
    delimiter_tmp = ','
    print('Writing feature history:',filename_tmp)
   
    # Modify reconstruction data according to the time-mean value subtraction
    # --Case for time mean subtraction
    if config['flag_subtraction_mean'] :
      # --Case for reconstuction data including time-mean data
      if config['flag_dmd_reconstruct_add_mean']: 
        displacement = displacement - displacement_mean[np.newaxis :, :]

    # Time mean value subtraction case
    displacement_total = np.sum(displacement,axis=0)
    displacement_reconst_mod = displacement_reconst_accum
    if config['flag_subtraction_mean'] :
      displacement_total = displacement_total + displacement_mean
      displacement_reconst_mod = displacement_reconst_mod + displacement_mean

    # File output
    file_input = open( filename_tmp, 'w')
    header_tmp = 'variables=Time[s]'
    for k in range(0,num_feature):
      header_tmp = header_tmp + ', Feature'+str(k+1)
    header_tmp = header_tmp +'\n'
    file_input.write( header_tmp )

    # All modes
    text_tmp='zone T="Allmodes", i='+str(num_step)+' f=point'+'\n'
    for m in range(0,num_step):
      text_force = ''
      for k in range(0,num_feature):
        text_force = text_force + str( displacement_total[m,k] ) + delimiter_tmp
      text_tmp = text_tmp + str(time[m]) + delimiter_tmp + text_force.rstrip(delimiter_tmp) + '\n'
    file_input.write( text_tmp )

    # Reconstructed modes
    for n in range(0,num_rank_consider):
      text_tmp='zone T="Reconstructed_numaccum_'+str(n+1)+'", i='+str(num_step)+' f=point'+'\n'
      for m in range(0,num_step):
        text_force = ''
        for k in range(0,num_feature):
          text_force = text_force + str( displacement_reconst_mod[n,m,k] ) + delimiter_tmp
        text_tmp = text_tmp + str(time[m]) + delimiter_tmp + text_force.rstrip(delimiter_tmp) + '\n'
      file_input.write( text_tmp )
    file_input.close()

    return
