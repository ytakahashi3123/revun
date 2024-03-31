#!/usr/bin/env python3

# Noise test

from modal.modal import modal
from noise.noise import noise

def main():

  # Read controlfile
  file_control = "config.yml"
  config = modal.read_config_yaml(file_control)

  # Check control file
  modal.check_config(config)
  
  # Set time variables
  time_dict = modal.set_time_variables(config)

  # VTK/VTU files reading
  noise.routine_add_noise(config, time_dict)

  return


if __name__ == '__main__':

  # call classes
  modal = modal()
  noise = noise()

  # main 
  main()

  exit()
