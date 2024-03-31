#!/usr/bin/env python3

# Routine for DMD mode sensing based on featured values by greedy algorithm

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
  reader_init, num_element, num_component = modal.read_vtk_initial(config, time_dict)
  #flowvar_t, reader_init, num_element, num_component = modal.read_vtk(config, time_dict)

  # Sample VTK file (reader_input is updated)
  if config['flag_read_samplevtk']:
    reader_init = modal.get_vtk_sample(config)

  # Read oscillation variables of modes
  frequency, amplitude = modal.read_DMD_oscillation(config, time_dict)
  
  # Calculate aerodynamic forces
  num_feature, feature_var, feature_var_mean = feature.read_feature_history(config, time_dict, num_component)

  # Run mode sensing approach
  feature_reconst = feature.mode_sensing_greedy( config, time_dict, num_feature, feature_var, feature_var_mean, frequency, amplitude)

  # Output
  feature.write_feature_history_reconstructed(config, time_dict, num_feature,feature_var, feature_var_mean, feature_reconst)

  return


if __name__ == '__main__':

  print('Initializing mode sensing routine based on feature value')

  # Call Class
  modal = modal()
  feature = feature()

  # main 
  main()

  print('Finalizing mode sensing routine based on feature value')

  exit()
