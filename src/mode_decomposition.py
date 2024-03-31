#!/usr/bin/env python3

# REVUN: Reconstruction script for VTK/VTU data format based on numerical mode decomposition
# Mode decomposition routine for VTK unstructured grid
# Version: 0.80
# Author: Yusuke Takahashi, Hokkaido University
# Contact: ytakahashi@eng.hokudai.ac.jp
# Date: 2023/08/21

import os as os
from modal.modal import modal

version = "0.80"

def main():

  # Read controlfile
  file_control = "config.yml"
  config = modal.read_config_yaml(file_control)

  # Check control file
  modal.check_config(config)

  # OpenMP number of theads in calculating eigen value and vector
  #os.environ["OMP_NUM_THREADS"] = str(config['omp_num_threads'])

  # Select the mode decomposition to run: snapshotPOD or DMD
  flag_mode = modal.check_mode(config)
  
  # Set time variables
  time_dict = modal.set_time_variables(config)

  # VTK/VTU files reading
  flowvar_t, reader_init, num_element, num_component = modal.read_vtk(config, time_dict)

  # Get time-mean value
  flowvar_mean_t = modal.get_mean(flowvar_t)

  # Subtract time-mean value
  if config['flag_subtraction_mean']:
    flowvar_t = modal.subtract_mean(time_dict, flowvar_t, flowvar_mean_t)

  # Sample VTK file (reader_input is updated)
  if config['flag_read_samplevtk']:
    reader_init = modal.get_vtk_sample(config)

  # Delete initial array of variables in VTK
  reader_init = modal.delete_initialarry_variables_vtk(reader_init)

  # Write time mean flow variables
  if config['flag_subtraction_mean']:
    modal.write_flowvar_mean(config, num_element, num_component, reader_init, flowvar_mean_t)

  # Snapshot proper orthogonal decomposition (POD)
  if flag_mode['snapshotpod']:
    # Get eigen value and vector
    eigen_value, eigen_vector, temporal_coef = modal.get_eigen(flowvar_t)
    # Extract modes
    flowvar_mode = modal.extract_POD_mode(config, time_dict, num_element, num_component, eigen_vector)
    # Write modes
    modal.write_POD_mode(config, time_dict, num_component, reader_init, eigen_value, temporal_coef, flowvar_mode)

  # Dynamic mode decompositin (DMD)
  if flag_mode['dmd']:
    # Build DMD modes
    singular_value, singular_value_full, eigen_value, eigen_vector, amplitude_coef, time_evolution = modal.get_DMD_mode(config, time_dict, flowvar_t)
    # Extract modes
    spatial_mode = modal.extract_DMD_mode(config, time_dict, num_element, num_component, eigen_vector)
    # Write modes
    modal.write_DMD_mode(config, time_dict, num_element, num_component, reader_init, flowvar_mean_t, singular_value, singular_value_full, eigen_value, spatial_mode, amplitude_coef, time_evolution)

  print('Mode decomposition finished.')

  return


if __name__ == '__main__':

  print('Mode decomposition routine for VTK/VTU file in unstructure grid format')
  print('--Version: ', version)

  print('Initializing modal analysis')

  # Class orbital
  modal = modal()

  # main 
  main()

  print('Finalizing modal analysis')

  exit()
