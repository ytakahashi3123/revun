#!/usr/bin/env python3

# Probe data

from modal.modal import modal
from probe.probe import probe

def main():

  # Read controlfile
  file_control = "config.yml"
  config = modal.read_config_yaml(file_control)

  # Check control file
  modal.check_config(config)
  
  # Set time variables
  time_dict = modal.set_time_variables(config)

  # Get probe data at points
  probe.get_probe_data(config, time_dict)

  return


if __name__ == '__main__':

  # call classes
  modal = modal()
  probe = probe()

  # main 
  main()

  exit()
