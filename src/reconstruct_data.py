#!/usr/bin/env python3

# Reconstruction routine based on DMD results

from modal.modal import modal
from reconstruction.reconstruction import reconstruction

def main():

  # Read controlfile
  file_control = "config.yml"
  config = modal.read_config_yaml(file_control)

  # Check control file
  modal.check_config(config)

  # Select the mode decomposition to run: snapshotPOD or DMD
  flag_mode = modal.check_mode(config)
  
  # Set time variables
  time_dict = modal.set_time_variables(config)

  # VTK/VTU files reading
  reader_init, num_element, num_component = modal.read_vtk_initial(config, time_dict)

  # Sample VTK file (reader_input is updated)
  if config['flag_read_samplevtk']:
    reader_init = modal.get_vtk_sample(config)

  # Read time-mean value
  spatial_mean = modal.read_flowvar_mean(config, time_dict, num_element, num_component)

  # Read DMD result
  spatial_mode, time_evolution = modal.read_DMD_mode(config, time_dict, num_element, num_component)

  # Reconstruct flow field from DMD modes
  reconstruction.reconstruct_DMD_mode(config, time_dict, num_element, num_component, reader_init, spatial_mean, spatial_mode, time_evolution)

  return


if __name__ == '__main__':

  print('Reconstruction routine based on DMD results')

  # call classes
  modal = modal()
  reconstruction = reconstruction()

  # main 
  main()

  exit()
