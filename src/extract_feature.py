#!/usr/bin/env python3

# Routine to extract feature on DMD spatiotemporal data

from modal.modal import modal
from feature.feature import feature

def main():

  # Read controlfile
  file_control = "config.yml"
  config = modal.read_config_yaml(file_control)

  # Check control file
  modal.check_config(config)
  
  # Check reconstruction flag
  feature.check_reconstruction(config)

  # Set time variables
  time_dict = modal.set_time_variables(config)

  # VTK/VTU files reading
  #flowvar_t, reader_init, num_element, num_component = modal.read_vtk(config, time_dict)
  reader_init, num_element, num_component = modal.read_vtk_initial(config, time_dict)

  # Sample VTK file (reader_input is updated)
  if config['flag_read_samplevtk']:
    reader_init = modal.get_vtk_sample(config)

  # Read time-mean value
  spatial_mean = modal.read_flowvar_mean(config, time_dict, num_element, num_component)

  # Read DMD result
  spatial_mode, time_evolution = modal.read_DMD_mode(config, time_dict, num_element, num_component)
  
  # Extract feature values
  num_feature, feature_var, feature_var_mean = feature.extract_feature_DMD(config, time_dict, num_component, reader_init, spatial_mean, spatial_mode, time_evolution)
  
  # Write feature value
  feature.write_feature_history(config, time_dict, num_feature, feature_var, feature_var_mean)

  return


if __name__ == '__main__':

  print('Initializing feature-extracting routine based on DMD spatiotemporal data')

  # Call Class
  modal = modal()
  feature = feature()

  # main 
  main()

  print('Finalizing feature-extracting routine based on DMD spatiotemporal data')

  exit()
