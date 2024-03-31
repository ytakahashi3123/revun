#!/usr/bin/env python3

import numpy as np
from scipy.linalg import svd
import os as os
import shutil as shutil
import yaml as yaml
import sys as sys
import vtk as vtk
from vtk.util import numpy_support
from vtk.numpy_interface import dataset_adapter
from general.general import general


class modal(general):

  def __init__(self):
    print("Calling class: modal")

    self.KIND_MODE_SUPPORT = ['snapshotpod', 'dmd']

  def check_config(self, config):
    print('Check config file')
    if not config['flag_subtraction_mean']:
      print('--Mean value is not subtracted.')
      if config['flag_dmd_reconstruct_add_mean']:
        print('--Option adding mean value to reconstucted data in reconstruction proces is valid. But this is not reasonable.')
        print('--Check "flag_subtraction_mean" and "flag_dmd_reconstruct_add_mean" in config file.')
        print('--Program stopped.')
        exit()
    return

  def check_mode(self, config):
    # Mode decomposition which is performed in this routine
    kind_mode_performed = config['kind_mode_performed']
    flag_mode={}
    flag_tmp=False
    for n in range(0,len(self.KIND_MODE_SUPPORT)):
      flag_mode[ self.KIND_MODE_SUPPORT[n].lower() ] = False 
      for m in range(0,len(kind_mode_performed)):
        if kind_mode_performed[m].lower() == self.KIND_MODE_SUPPORT[n].lower() :
          flag_mode[ self.KIND_MODE_SUPPORT[n].lower() ] = True
          flag_tmp = True
          break
    if not flag_tmp : 
      print( 'Any modes are not performed. Check kind_mode_performed in configuration file.' )
      print( 'Program stopped.' )
      exit()
    return flag_mode

  def check_xml(self, file_name_base):
    # Check file is XML or not by extension name (.vtk or .vtu)
    filename_ext=os.path.splitext(file_name_base)[-1]
    if filename_ext=='.vtu' :
      flag_xml = True
    elif filename_ext=='.vtk':
      flag_xml = False
    else:
      print('File type does not match VTK or VTU.')
      print('Program stopped.')
      exit()
    return flag_xml

  def get_extension(self, file_name_base):
    filename_ext = os.path.splitext(file_name_base)[-1]
    return filename_ext

  def set_time_variables(self, config):
    step_start      = config['step_start']
    step_end        = config['step_end']
    step_interval   = config['step_interval']
    num_step        = (step_end-step_start)//step_interval
    time_step       = config['time_step']
    time_step_modal = time_step*float(step_interval)
    time            = np.linspace(0, time_step_modal*(num_step-1), num_step)

    time_dict     = {'step_start':step_start,
                    'step_end': step_end,
                    'step_interval': step_interval,
                    'num_step': num_step,
                    'time_step': time_step_modal,
                    'time': time
                    }

    return time_dict

  def get_number_rank(self, config, time_dict):
    num_step               = time_dict['num_step']
    mode_dmd_index_maximum = config['mode_dmd_index_maximum']
    num_rank               = min(num_step, mode_dmd_index_maximum)-1
    return num_rank

  def get_vtk_reader(self,filename,filename_ext):
    if filename_ext == '.vtu' :
      reader = vtk.vtkXMLUnstructuredGridReader()
    elif filename_ext == '.vtk' :
      reader = vtk.vtkUnstructuredGridReader()
      reader.ReadAllScalarsOn()
      reader.ReadAllVectorsOn()
      reader.ReadAllTensorsOn()
    elif filename_ext == '.vti' :
      reader = vtk.vtkXMLImageDataReader()
    elif filename_ext == '.png' :
      reader = vtk.vtkPNGReader()
    else:
      print('This file formant is not supported:',filename_ext)
      print('Program stipped.')
      exit()
    reader.SetFileName(filename)
    reader.Update()
    return reader

  def get_vtk_writer(self,filename,filename_ext):
    if filename_ext == '.vtu' :
      writer = vtk.vtkXMLUnstructuredGridWriter()
    elif filename_ext == '.vtk' :
      writer = vtk.vtkUnstructuredGridWriter()
    elif filename_ext == '.vti' :
      writer = vtk.vtkXMLImageDataWriter()
    else:
      print('This file formant is not supported:',filename_ext)
      print('Program stipped.')
      exit()
    writer.SetFileName(filename)
    return writer

  def read_vtk(self, config,time_dict):
    # Read VTK(U) files
    print('Start reading VTK(U) files')

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

    file_name_base  = config['dir_input'] +'/'+ config['file_input_base']
    #flag_xml = self.check_xml(file_name_base)
    index_inc = 0

#    for n in range(step_start, step_end, step_interval):
    for n in range(0, num_step):
      nstep        = step_start + n*step_interval
      filename_tmp = self.split_file(file_name_base,'_'+str(nstep).zfill(step_digit),'.')
      filename_ext = self.get_extension(filename_tmp)
      reader       = self.get_vtk_reader(filename_tmp,filename_ext)
      
      # Number of cells
      #num_cell=reader.GetNumberOfCells()
      # Number of Points
      #num_node=reader.GetNumberOfPoints()

      # Extract flow variable
      flowvar_tmp = numpy_support.vtk_to_numpy( reader.GetOutput().GetPointData().GetAbstractArray(kind_flowvar) )
      if flag_flowvar_replace:
        # Velocity=Momentum/Density
        var_replace_bot = numpy_support.vtk_to_numpy( reader.GetOutput().GetPointData().GetAbstractArray(varname_replaced_bot) )
        flowvar_tmp = flowvar_tmp/var_replace_bot[:,None]

      # Initial settings
      if n==0:
        # --Get sample VTK file
        reader_init = reader
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
        # --Set initial array
        flowvar_t = np.zeros(num_step*num_element*num_component).reshape(num_step,num_element*num_component)
        # --Mean value
        #flowvar_mean_t = np.zeros(num_element*num_component).reshape(num_element*num_component)

      # When using flowvar_tmp.ravel() for multi-component such velocities (u, v, w), elements are stored to one-dimension array like u0, v0, w0, u1, u2, w2, u3..., 
      # because "u", "v", and "w" are stored to flowvar_tmp[:][0],flowvar_tmp[:][1] and flowvar_tmp[:][2].
      # To set u0, u1, u2.., v1, v2, v3..., w0 ..., transpose methos is used
      flowvar_trans          = flowvar_tmp.transpose()
      flowvar_t[n,:] = flowvar_trans.ravel()[:]

      # Calculate time-mean value by on-line algorithm
      #flowvar_mean_t = flowvar_mean_t + (flowvar_t[index_inc,:]-flowvar_mean_t)/float(index_inc+1)

      #index_inc = index_inc+1
      print('Read input file: ',filename_tmp)

      del reader

    print('Finish reading files')

    return flowvar_t, reader_init, num_element, num_component

  def read_vtk_initial(self, config, time_dict):
    # Used in reconstruction
    print('Start reading VTK(U) files')

    step_start    = time_dict['step_start']
    step_digit    = config['step_digit']

    # Target flow variable to run POD
    kind_flowvar         = config['kind_flowvar']
    #flag_flowvar_replace = config['kind_flowvar_replace'][0]
    #if flag_flowvar_replace:
    #  varname_replaced_top = config['kind_flowvar_replace'][1]
    #  varname_replaced_bot = config['kind_flowvar_replace'][2]

    file_name_base = config['dir_input'] +'/'+ config['file_input_base']
    filename_tmp   = self.split_file(file_name_base,'_'+str(step_start).zfill(step_digit),'.')
    filename_ext   = self.get_extension(filename_tmp)
    reader         = self.get_vtk_reader(filename_tmp,filename_ext)

      # Extract flow variable
    flowvar_tmp = numpy_support.vtk_to_numpy( reader.GetOutput().GetPointData().GetAbstractArray(kind_flowvar) )
    #if flag_flowvar_replace:
    #  # Velocity=Momentum/Density
    #  var_replace_bot = numpy_support.vtk_to_numpy( reader.GetOutput().GetPointData().GetAbstractArray(varname_replaced_bot) )
    #  flowvar_tmp = flowvar_tmp/var_replace_bot[:,None]

    # --Get sample VTK file
    #reader_init = reader
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

    return reader, num_element, num_component

  def get_arrayname(self, reader):
    num_of_array_vtk = reader.GetOutput().GetPointData().GetNumberOfArrays()
    arrayname=[]
    for n in range(0,num_of_array_vtk):
      arrayname.append(reader.GetOutput().GetPointData().GetArrayName(n))
    return arrayname

  def get_vtk_sample(self, config):
    # Read sample VTK(U) file 
    filename_tmp=config['dir_read_samplevtk'] +'/'+config['file_read_samplevtk']
    flag_xml = self.check_xml(filename_tmp)
    print('Reading input sample file: ',filename_tmp)
    reader = self.get_vtk_reader(filename_tmp,flag_xml)
    return reader

  def delete_initialarry_variables_vtk(self, reader):
    # Delete initial data-array varibales of sample VTK file
    name_list=self.get_arrayname(reader)
    for n in range(0,len(name_list)):
      arrayname=name_list[n]
      reader.GetOutput().GetPointData().RemoveArray(arrayname)
    return reader

  def get_mean(self, flowvar_t):
    # Calculate time-mean value
    print('Calculate mean values')
    flowvar_mean_t = np.mean(flowvar_t,axis=0)
    return flowvar_mean_t

  def subtract_mean(self, time_dict, flowvar_t, flowvar_mean_t):
    print('Subtract mean values')
    num_step = time_dict['num_step']
    for n in range(0,num_step):
      flowvar_t[n,:] = flowvar_t[n,:]-flowvar_mean_t[:]
    return flowvar_t

  def get_eigen(self, flowvar_t):
    # Calculate eigen value and vector

    print('Start to calculate eigen velues and vectors')

    flowvar = flowvar_t.T
    product = np.dot(flowvar_t, flowvar)

    # Eigen decomposition
    eigen_value_tmp, eigen_vector_tmp = np.linalg.eigh(product)

    # Sort
    index_sort = eigen_value_tmp.argsort()[::-1]
    eigen_value_sort  = eigen_value_tmp[index_sort]  #eigenvalue vector
    eigen_vector_sort = eigen_vector_tmp[:,index_sort]

    eigen_value   = eigen_value_sort
    lamda_sqrtinv = np.diag( eigen_value_sort**(-1/2) )
    lamda_sqrt    = np.diag( eigen_value_sort**( 1/2) )

    # Compute POD modes and temporary coefficient
    eigen_vector  = np.dot( np.dot(flowvar, eigen_vector_sort), lamda_sqrtinv)
    temporal_coef = np.dot(lamda_sqrt, eigen_vector_sort.T)
    #print(lamda_sqrt.shape, eigen_vector_tmp.T.shape)

    return eigen_value, eigen_vector, temporal_coef 

  def extract_POD_mode(self, config,time_dict,num_element,num_component,eigen_vector):
    num_step = time_dict['num_step']
    flowvar_mode = np.zeros(num_element*num_component*num_step).reshape(num_element,num_component,num_step)
    for i in range(0,num_step):
      for j in range(0,num_component):
        flowvar_mode[0:num_element,j,i] = np.matrix(eigen_vector[num_element*j:num_element*(j+1),i])
    return flowvar_mode

  def write_POD_mode(self, config, time_dict, num_component, reader_init, eigen_value, temporal_coef, flowvar_mode):

    # Make file output directory
    dir_result_path = config['dir_result_pod']
    if not os.path.exists(dir_result_path):
      os.mkdir(dir_result_path)

    # Initial settings
    step_start    = time_dict['step_start']
    step_end      = time_dict['step_end']
    step_interval = time_dict['step_interval']
    num_step      = time_dict['num_step']
    time          = time_dict['time']

    delimiter_tmp = ','

    # Set index
    mode_pod_index_maximum = config['mode_pod_index_maximum']
    mode_max               = min(num_step, mode_pod_index_maximum)
    mode_display_digit     = len( str( mode_max ) )

    # Eigen value
    if config['flag_result_pod_eigenvalue'] :
      filename_tmp = config['dir_result_pod']+'/'+self.split_file(config['file_result_pod_eigenvalue'],'_'+config['kind_flowvar'],'.')
      print('Writing eigen values: ', filename_tmp)
      file_output = open( filename_tmp, 'w')
      text_tmp  = 'variables=Mode, Eigen_value, Eigen_value_ratio'+'\n'
      eigen_value_sum = 1e-15
      for n in range(0,num_step):
        eigen_value_sum = eigen_value_sum + eigen_value[n]
      eigen_value_deno = 0.0
      for n in range(0,num_step):
        eigen_value_deno = eigen_value_deno + eigen_value[n]
        text_tmp = text_tmp + str(n+1) + delimiter_tmp + str(eigen_value[n]) + delimiter_tmp + str(eigen_value_deno/eigen_value_sum) + '\n'
      file_output.write( text_tmp )
      file_output.close()

    # Temporal coefficient
    if config['flag_result_pod_temporal'] :
      filename_tmp = config['dir_result_pod']+'/'+self.split_file(config['file_result_pod_temporal'],'_'+config['kind_flowvar'],'.')
      print('Writing temporal coefficients: ', filename_tmp)
      file_output = open( filename_tmp, 'w')
      header_tmp  = 'variables=Time[s],Temporal_coefficient'+'\n'
      file_output.write( header_tmp )
      for n in range(0,mode_max):
        text_tmp='zone T="Mode'+str(n+1).zfill(mode_display_digit)+'", i='+str(num_step)+' f=point'+'\n'
        for m in range(0,num_step):
          text_tmp = text_tmp + str(time[m]) + delimiter_tmp + str(temporal_coef[n][m])+'\n'
        file_output.write( text_tmp )
      file_output.close()

    # Spatial mode
    if config['flag_result_pod_spatial']:
      # Add new data array
      value_name_base = config['kind_flowvar'] 
      if config['kind_flowvar_replace'][0]:
        value_name_base = config['kind_flowvar_replace'][1]

      add_value = dataset_adapter.WrapDataObject(reader_init.GetOutput())
      for n in range(0,mode_max):
        value_name = value_name_base + "_mode_" +str(n+1).zfill(mode_display_digit)
        NewDataSet = flowvar_mode[:,0:num_component,n]
        add_value.PointData.append(NewDataSet, value_name)

      # Write VTK
      filename_tmp = config['dir_result_pod']+'/'+ self.split_file(config['file_result_pod_spatial'],'_'+config['kind_flowvar'],'.')
      filename_ext = self.get_extension(filename_tmp)
      print('Writing mode file: ', filename_tmp)
      writer = self.get_vtk_writer(filename_tmp,filename_ext)
      #if flag_xml :
      #  writer = vtk.vtkXMLUnstructuredGridWriter()
      #else :
      #  writer = vtk.vtkUnstructuredGridWriter()
      #writer.SetFileName(filename_tmp)
      writer.SetInputData(add_value.VTKObject)
      writer.Write()

    return

  def get_DMD_mode(self, config, time_dict, flowvar_t):
    # Extract DMD modes

    print('Start to build DMD modes')

    step_start    = time_dict['step_start']
    step_end      = time_dict['step_end']
    step_interval = time_dict['step_interval']
    num_step      = time_dict['num_step']
    dt            = time_dict['time_step']
    time          = time_dict['time']

    # Make data matrix
    D = flowvar_t.T

    # Make DMD input/output matix (Y=AX)
    X = D[:,:-1] # D要素2つ目[0,1,...,jmax-1]
    Y = D[:,1:]  # D要素2つ目[1,2,...,jmax]

    # Singular value decomposition for input matrix X
    # X = U \Sigma V*
    print('--Run singular value decomposition for input matrix')
    U2,S2,Vh2 = svd(X, False)
    # U2:  Unitary matrix (m, n): U^{*}U=UU^{*}=I (U^{*}:Adjoint matrix)
    # Vh2: Adjoint matrix of unitary matrix (n, n)
    # Sig2: Singular values
    singular_value_full = S2

    # Truncate rank
    mode_dmd_index_maximum = config['mode_dmd_index_maximum']
#    num_rank               = min(num_step-1, mode_dmd_index_maximum)
    num_rank               = min(num_step, mode_dmd_index_maximum)-1
    print('--Maximum modes considered:', num_rank+1)

    unitary        = U2[:,:num_rank]
    singular_value = np.diag(S2)[:num_rank,:num_rank]
    adjoint        = Vh2.conj().T[:,:num_rank]

    # Reconstruct A~ which is the projection of the full matrix 
    # A~ = U^{*} A U = U^{*} Y V Sigma^{-1}
    # Here, the eigenvalues for A are equivalent to those of A~.
    print('--Reconstruct eigen value and eigen vector for projection of the full matrix (A~) and the full matrix (A).') 
    print('--The eigenvalues for A are equivalent to those of A~')
    A_tilde = np.dot(np.dot(np.dot(unitary.conj().T, Y), adjoint), np.linalg.inv(singular_value))

    # Computed eigen value and eigen vector of A~ (eigen_value: mu, eigen_vec)
    eigen_value, eigen_vector_tilde = np.linalg.eig(A_tilde)

    # Build DMD mode
    # Reconstruct the eigendecomposition of A from W and \Lambda.
    # The eigenvectors of "A" are given by the columns of "\Phi". Here, "A \Phi = \Phi \Lambda" --> "\Phi=Y V \Sigma^{-1} W"
    eigen_vector = np.dot( np.dot( np.dot(Y, adjoint), np.linalg.inv(singular_value)), eigen_vector_tilde)

    # Compute time evolution
    print('--Reconstruct the system time evolution')
    # Reconstruct a matrix "\Psi" corresponding to the system’s time evolution in continuous time
    amplitude_coef = np.dot(np.linalg.pinv(eigen_vector), X[:,0])
    time_evolution = np.zeros([num_rank, len(time)], dtype='complex')
    for i,_t in enumerate(time):
      time_evolution[:,i] = np.multiply(np.power(eigen_value, _t/dt), amplitude_coef)

    # Sort (not recommended)
    if config['flag_dmd_sort']:
      print('--Sort mode ID based on amplitude')
      index_sort          = np.abs(amplitude_coef).argsort()[::-1]
      # When num_rank < num_step, indexes are sorted in num_rank, and then, the indexes of the remaining ranks are added.
      l = np.linspace(0, num_step-2, num_step-1, dtype='int')
      index_sort_tmp = np.append(index_sort, l).tolist()
      index_sort_tmp = sorted(set(index_sort_tmp), key=index_sort_tmp.index)
      index_sort_alt = np.array(index_sort_tmp)
      # Sort
      singular_value      = singular_value[index_sort] 
      singular_value_full = singular_value_full[index_sort_alt]
      eigen_value         = eigen_value[index_sort] 
      eigen_vector        = eigen_vector[:,index_sort]
      amplitude_coef      = amplitude_coef[index_sort]
      time_evolution      = time_evolution[index_sort,:]

    # DMD再構成の計算
    if config['flag_verbose']:
      print('--Check whether the reconstructed data reproduces the original data set')
      D2 = np.dot(eigen_vector, time_evolution)
      #print(np.dot(eigen_vector[:,0:1], time_evolution[0:1,:]))
      print( '--',np.allclose(D, D2) ) # True

    print('Finish to build DMD modes')

    return singular_value, singular_value_full, eigen_value, eigen_vector, amplitude_coef, time_evolution

  def extract_DMD_mode(self, config,time_dict,num_element,num_component,eigen_vector):
    print('Extract DMD mode (spatial)')
    num_step               = time_dict['num_step']
    mode_dmd_index_maximum = config['mode_dmd_index_maximum']
    #num_rank               = min(num_step-1, mode_dmd_index_maximum)
    num_rank               = min(num_step, mode_dmd_index_maximum)-1
    flowvar_mode = np.zeros(num_element*num_component*num_rank, dtype='complex').reshape(num_element,num_component,num_rank)
    for i in range(0,num_rank):
      for j in range(0,num_component):
        #flowvar_mode[0:num_element,j,i] = np.matrix(np.real( eigen_vector[num_element*j:num_element*(j+1),i]))
        flowvar_mode[0:num_element,j,i] = eigen_vector[num_element*j:num_element*(j+1),i]
    return flowvar_mode

  def write_DMD_mode(self, config, time_dict, num_element, num_component, reader_init, flowvar_mean_t, singular_value, singular_value_full, eigen_value, spatial_mode, amplitude_coef, time_evolution):
    # Output results

    # Make file output directory
    self.make_directory( config['dir_result_dmd'] )

    # Inital settings
    step_start    = time_dict['step_start']
    step_end      = time_dict['step_end']
    step_interval = time_dict['step_interval']
    num_step      = time_dict['num_step']
    dt            = time_dict['time_step']
    t             = time_dict['time']

    delimiter_tmp = ','

    mode_dmd_index_maximum = config['mode_dmd_index_maximum']
    num_rank               = min(num_step, mode_dmd_index_maximum)-1

    # Set index
    mode_display_digit = len( str( num_rank ) )

    # DMD variables 
    frequency = np.imag( np.log( eigen_value )/dt)/(2.0*np.pi)
    amplitude = np.abs(amplitude_coef)
    index_sort_amp = amplitude.argsort()[::-1]
    index_mode_tmp = np.linspace(0, num_rank-1, num_rank).astype(int)

    # Display DMD variables at each mode
    print('Displaying DMD variables at each mode...')
    print('Rank, Frequency, Amplitude, Abs.Eigenvalue, Rank before sorted by Amplitude')
    space_tmp = '  '
    for n in range(0,num_rank):
#      m = index_mode_tmp[index_sort_amp[n]]
      m = index_sort_amp[n]
      print( space_tmp+str(n+1)
            +space_tmp+str('{:e}'.format(frequency[m]))
            +space_tmp+str('{:e}'.format(amplitude[m]))
            +space_tmp+str('{:e}'.format(np.abs(eigen_value[m]))) 
            +space_tmp+str(m+1)
            #+space_tmp+str((index_sort_amp[n])) 
            )

    # Sigunlar values distribution
    if config['flag_result_dmd_singularvalue'] :
      filename_tmp = config['dir_result_dmd']+'/'+self.split_file(config['file_result_dmd_singularvalue'],'_'+config['kind_flowvar'],'.')
      print('Writing singular values: ', filename_tmp)
      file_output = open( filename_tmp, 'w')
      text_tmp  = 'variables=Mode, Singular_value, Singular_value_ratio'+'\n'
      singular_value_sum = 1e-15
      for n in range(0,num_step-1):
        singular_value_sum = singular_value_sum + singular_value_full[n]
      singular_value_deno = 0.0
      for n in range(0,num_step-1):
        singular_value_deno = singular_value_deno + singular_value_full[n]
        text_tmp=text_tmp + str(n+1) + delimiter_tmp + str(singular_value_full[n]) + delimiter_tmp + str(singular_value_deno/singular_value_sum) + '\n'
      file_output.write( text_tmp )
      file_output.close()

    # Eigen value distribution
    if config['flag_result_dmd_eigenvalue'] :
      filename_tmp = config['dir_result_dmd']+'/'+self.split_file(config['file_result_dmd_eigenvalue'],'_'+config['kind_flowvar'],'.')
      print('Writing eigen values: ', filename_tmp)
      file_output = open( filename_tmp, 'w')
      header_tmp  = 'variables=Real, Imaginary, Absolute, Argument, Mode'+'\n'
      file_output.write( header_tmp )
      # Unit circle
      degree_start =-180
      degree_end   = 180
      degree_div   = 360
      c_x,c_y,c_abs,c_arg = [],[],[],[]
      for _x in np.linspace(degree_start,degree_end,degree_div):
        c_x.append(np.sin(np.radians(_x)))
        c_y.append(np.cos(np.radians(_x)))
        c_abs.append(1.0)
        c_arg.append(_x)
      text_tmp ='zone T="Unit_circle", i='+str(degree_div)+' f=point'+'\n'
      for n in range(0,degree_div):
        text_tmp = text_tmp + str( c_x[n] )   + delimiter_tmp \
                            + str( c_y[n] )   + delimiter_tmp \
                            + str( c_abs[n] ) + delimiter_tmp \
                            + str( c_arg[n] ) + delimiter_tmp \
                            + str(1) + '\n' # Mode: Dummy
      file_output.write( text_tmp )
      # Eigen value--Discrete spectrum
      text_tmp = 'zone T="Discrete_Spectrum_Mode", i='+str(num_rank)+' f=point'+'\n'
      for n in range(0,num_rank):
        abs_tmp = np.abs(eigen_value[n])
        arg_tmp = np.degrees( np.arctan2(np.imag(eigen_value[n]), np.real(eigen_value[n])) )
        text_tmp = text_tmp + str( np.real(eigen_value[n]) ) + delimiter_tmp \
                            + str( np.imag(eigen_value[n]) ) + delimiter_tmp \
                            + str(abs_tmp) + delimiter_tmp \
                            + str(arg_tmp) + delimiter_tmp \
                            + str(n+1) + '\n'
      file_output.write( text_tmp )
      # Eigen value--Continuous spectrum
      text_tmp = 'zone T="Continuous_Spectrum_Mode", i='+str(num_rank)+' f=point'+'\n'
      for n in range(0,num_rank):
        eigen_value_cont = np.log(eigen_value[n])/dt
        abs_tmp = np.abs( eigen_value_cont)
        arg_tmp = np.degrees( np.arctan2(np.imag( eigen_value_cont), np.real( eigen_value_cont)) )
        text_tmp = text_tmp + str( np.real( eigen_value_cont ) ) + delimiter_tmp \
                            + str( np.imag( eigen_value_cont ) ) + delimiter_tmp \
                            + str(abs_tmp) + delimiter_tmp \
                            + str(arg_tmp) + delimiter_tmp \
                            + str(n+1) + '\n'
      file_output.write( text_tmp )
      # File closing
      file_output.close()

    # Osillation: Amplitude and frequency
    if config['flag_result_dmd_oscillation']:
      filename_tmp = config['dir_result_dmd']+'/'+self.split_file(config['file_result_dmd_oscillation'],'_'+config['kind_flowvar'],'.')
      print('Writing amplitudes and frequencies:',filename_tmp)
      file_output = open( filename_tmp, 'w')
      header_tmp  = 'variables=Mode, Frequency[Hz], Amplitude'+'\n'
      file_output.write( header_tmp )
      text_tmp = ''
      for n in range(0,num_rank):
        text_tmp = text_tmp + str(n+1) + delimiter_tmp + str( frequency[n] ) + delimiter_tmp + str( amplitude[n] ) +'\n'
      file_output.write( text_tmp )
      file_output.close()

    # Temporal mode
    if config['flag_result_dmd_temporal'] :
      filename_tmp = config['dir_result_dmd']+'/'+self.split_file(config['file_result_dmd_temporal'],'_'+config['kind_flowvar'],'.')
      print('Writing temporal mode:',filename_tmp)
      file_output = open( filename_tmp, 'w')
      header_tmp  = 'variables=Time[s], Real, Imaginary'+'\n'
      file_output.write( header_tmp )
      for n in range(0,num_rank):
        text_tmp='zone T="Mode'+str(n+1).zfill(mode_display_digit)+'", i='+str(num_step)+' f=point'+'\n'
        for m in range(0,num_step):
          text_tmp=text_tmp + str(t[m]) + delimiter_tmp + str(np.real(time_evolution[n][m]))+ delimiter_tmp +str(np.imag(time_evolution[n][m]))+'\n'
        file_output.write( text_tmp )
      file_output.close()

    # Spatial mode
    if config['flag_result_dmd_spatial']:
      reader_tmp = self.delete_initialarry_variables_vtk(reader_init)
      # --Add new data array
      value_name_base = config['kind_flowvar'] 
      if config['kind_flowvar_replace'][0]:
        value_name_base = config['kind_flowvar_replace'][1]
      add_value = dataset_adapter.WrapDataObject(reader_tmp.GetOutput())
      for n in range(0,num_rank):
        # Real part
        value_name = value_name_base + "_Real_mode_" +str(n+1).zfill(mode_display_digit)
        NewDataSet = np.real( spatial_mode[:,0:num_component,n] )
        add_value.PointData.append(NewDataSet, value_name)
        # Imaginary part
        value_name = value_name_base + "_Imag_mode_" +str(n+1).zfill(mode_display_digit)
        NewDataSet = np.imag( spatial_mode[:,0:num_component,n] )
        # Append
        add_value.PointData.append(NewDataSet, value_name)
      # --Write VTK
      filename_tmp = config['dir_result_dmd']+'/'+ self.split_file(config['file_result_dmd_spatial'],'_'+config['kind_flowvar'],'.')
      filename_ext = self.get_extension(filename_tmp)
      print('Writing spatial mode: ', filename_tmp)
      writer = self.get_vtk_writer(filename_tmp,filename_ext)
      #if flag_xml :
      #  writer = vtk.vtkXMLUnstructuredGridWriter()
      #else :
      #  writer = vtk.vtkUnstructuredGridWriter()
      #writer.SetFileName(filename_tmp)
      writer.SetInputData(add_value.VTKObject)
      writer.Write()
    return


  def read_flowvar_mean(self, config, time_dict, num_element, num_component):
    print('Read time-mean flow varable data')
    num_step = time_dict['num_step']
    # Set array
    spatial_mean = np.zeros(num_element*num_component).reshape(num_element,num_component)
    if config['flag_subtraction_mean']:
      # Set variable name
      value_name = config['kind_flowvar'] 
      if config['kind_flowvar_replace'][0]:
        value_name = config['kind_flowvar_replace'][1]
      # Read mean data
      filename_tmp = config['dir_result_mean']+'/'+ self.split_file(config['file_result_mean'],'_'+config['kind_flowvar'],'.')
      filename_ext = self.get_extension(filename_tmp)
      reader       = self.get_vtk_reader(filename_tmp,filename_ext)
      if num_component == 1:
        spatial_mean[:,0] = numpy_support.vtk_to_numpy( reader.GetOutput().GetPointData().GetAbstractArray(value_name) )
      else :
        spatial_mean[:,:] = numpy_support.vtk_to_numpy( reader.GetOutput().GetPointData().GetAbstractArray(value_name) )
    return spatial_mean


  def write_flowvar_mean(self, config, num_element, num_component, reader_init, flowvar_mean_t):
    # Make file output directory
    self.make_directory( config['dir_result_mean'] )

    spatial_mean = np.zeros(num_element*num_component).reshape(num_element,num_component)
    for j in range(0,num_component):
      spatial_mean[0:num_element,j] = flowvar_mean_t[num_element*j:num_element*(j+1)]
    # Write mean data
    reader_tmp = self.delete_initialarry_variables_vtk(reader_init)
    value_name = config['kind_flowvar'] 
    if config['kind_flowvar_replace'][0]:
      value_name = config['kind_flowvar_replace'][1]
    add_value  = dataset_adapter.WrapDataObject(reader_tmp.GetOutput())
    NewDataSet = spatial_mean
    add_value.PointData.append(NewDataSet, value_name)
    filename_tmp = config['dir_result_mean']+'/'+ self.split_file(config['file_result_mean'],'_'+config['kind_flowvar'],'.')
    filename_ext = self.get_extension(filename_tmp)
    print('Writing time-mean flow varable: ', filename_tmp)
    writer = self.get_vtk_writer(filename_tmp,filename_ext)
    writer.SetInputData(add_value.VTKObject)
    writer.Write()
    return


  def read_DMD_mode(self, config, time_dict, num_element, num_component):

    print('Read DMD model results')

    # Inital settings
    step_start    = time_dict['step_start']
    step_end      = time_dict['step_end']
    step_interval = time_dict['step_interval']
    num_step      = time_dict['num_step']
    dt            = time_dict['time_step']
    time          = time_dict['time']

    delimiter_tmp = ','

    # Set index
    mode_dmd_index_maximum = config['mode_dmd_index_maximum']
    num_rank               = min(num_step, mode_dmd_index_maximum)-1
    mode_display_digit = len( str( num_rank ) )

    # Make array
    time_evolution = np.zeros([num_rank, num_step], dtype='complex')
    spatial_mode = np.zeros(num_element*num_component*num_rank, dtype='complex').reshape(num_element,num_component,num_rank)

    # Spatial mode
    if config['flag_result_dmd_spatial']:
      # Set variables name
      value_name_base = config['kind_flowvar'] 
      if config['kind_flowvar_replace'][0]:
        value_name_base = config['kind_flowvar_replace'][1]
      # Read spatial mode file
      filename_tmp = config['dir_result_dmd']+'/'+ self.split_file(config['file_result_dmd_spatial'],'_'+config['kind_flowvar'],'.')
      filename_ext = self.get_extension(filename_tmp)
      reader   = self.get_vtk_reader(filename_tmp,filename_ext)
      # Store DMD spatial mode variables 
      for n in range(0,num_rank):
        value_name = value_name_base + "_Real_mode_" +str(n+1).zfill(mode_display_digit)
        value_real = numpy_support.vtk_to_numpy( reader.GetOutput().GetPointData().GetAbstractArray(value_name) )
        value_name = value_name_base + "_Imag_mode_" +str(n+1).zfill(mode_display_digit)
        value_imag = numpy_support.vtk_to_numpy( reader.GetOutput().GetPointData().GetAbstractArray(value_name) )
        if num_component == 1:
          spatial_mode[:,0,n] = value_real.astype(float)+1j*value_imag.astype(float)
        else:
          spatial_mode[:,:,n] = value_real.astype(float)+1j*value_imag.astype(float)

    # Temporal mode
    if config['flag_result_dmd_temporal'] :
      filename_tmp = config['dir_result_dmd']+'/'+self.split_file(config['file_result_dmd_temporal'],'_'+config['kind_flowvar'],'.')
      print('--Reading temporal mode:',filename_tmp)
      file_output = open( filename_tmp, 'r')
      buffer = file_output.readlines()
      n_count  = 1
      for n in range(0,num_rank):
        for m in range(0,num_step):
           n_count=n_count+1
           buffr_tmp = buffer[n_count].split(',')
#           time_evolution[n][m] = complex( float(buffr_tmp[1]), float(buffr_tmp[2]))
           time_evolution[n][m] = float(buffr_tmp[1])+1j*float(buffr_tmp[2])
        n_count=n_count+1
      file_output.close()

    return spatial_mode, time_evolution


  def read_DMD_oscillation(self, config, time_dict):

    # Inital settings
    num_step      = time_dict['num_step']
    num_rank = self.get_number_rank(config, time_dict)

    frequency = np.zeros(num_rank)
    amplitude = np.zeros(num_rank)

    if config['flag_result_dmd_oscillation']:
      filename_tmp = config['dir_result_dmd']+'/'+self.split_file(config['file_result_dmd_oscillation'],'_'+config['kind_flowvar'],'.')
      print('Reading amplitudes and frequencies:',filename_tmp)
      file_input = open( filename_tmp, 'r')
      buffer = file_input.readlines()
      n_count  = 0
      for n in range(0,num_rank):
        n_count   = n_count+1
        buffr_tmp = buffer[n_count].split(',')
        frequency[n] = np.array(buffr_tmp[1]).astype(float)
        amplitude[n] = np.array(buffr_tmp[2]).astype(float)
      file_input.close()

    return frequency, amplitude