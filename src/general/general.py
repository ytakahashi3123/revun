#!/usr/bin/env python3

class general:

  def __init__(self):
    print("Calling class: general")

  def argument(self, filename_default):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-file', action='store', type=str, default=filename_default)
    args = parser.parse_args()
    return args

  def read_config_yaml(self, file_control):
    import yaml as yaml
    import sys as sys
    import pprint as pprint
    print("Reading control file...:", file_control)
    try:
      with open(file_control) as file:
        config = yaml.safe_load(file)
#        pprint.pprint(config)
    except Exception as e:
      print('Exception occurred while loading YAML...', file=sys.stderr)
      print(e, file=sys.stderr)
      sys.exit(1)
    return config

  def make_directory(self, dir_path):
    import os as os
    import shutil as shutil
    if not os.path.exists(dir_path):
      os.mkdir(dir_path)
    #else:
    #  shutil.rmtree(dir_path)
    #  os.mkdir(dir_path)
    return

  def split_file(self, filename,addfile,splitchar):
    # 'filename'の特定の文字'splitchar'の前に'addfile'を加える.'splitchar'が２つ以上ある場合は末尾優先（のはず）
    splitchar_tmp   = splitchar
    filename_split  = filename.rsplit(splitchar_tmp, 1)
    filename_result = filename_split[0]+addfile+splitchar+filename_split[1]
    return filename_result