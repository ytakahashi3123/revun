#!/usr/bin/env python3

# Routine to find contributions of modes based on reconstructed aerodynamic forces

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

  # Read oscillation variables of modes
  frequency, amplitude = modal.read_DMD_oscillation(config, time_dict)
  
  # Calculate aerodynamic forces
  num_feature, aerodynamic_force, aerodynamic_force_mean = feature.get_aerdynamic_force_hisotry(config, time_dict)

  # Run mode sensing approach
  feature_reconst = feature.mode_sensing_greedy( config, time_dict, num_feature, aerodynamic_force, aerodynamic_force_mean, frequency, amplitude)

  # Output
  feature.write_feature_history_reconstructed(config, time_dict, num_feature, aerodynamic_force, aerodynamic_force_mean, feature_reconst)

  return


if __name__ == '__main__':

  print('Initializing mode sensing routine based on aerodynamic forces reconstructed')

  # Call Class
  modal = modal()
  feature = feature()

  # main 
  main()

  print('Finalizing mode sensing routine based on aerodynamic forces reconstructed')

  exit()
