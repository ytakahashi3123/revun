#!/usr/bin/env python3

# Routine to filter image file in coarse graining

from modal.modal import modal
from filter_data.filter_data import filter_data

def main():

  # Read controlfile
  file_control = "config_support.yml"
  config = modal.read_config_yaml(file_control)
  
  # Set time variables
  time_dict = modal.set_time_variables(config)

  # Make image file coarse
  filter_data.filter_png(config, time_dict)

  return


if __name__ == '__main__':

  print('Initializing routine to make an image file coarse')

  # call classes
  modal = modal()
  filter_data = filter_data()

  # main 
  main()

  print('Finalizing routine to make an image file coarse')

  exit()
