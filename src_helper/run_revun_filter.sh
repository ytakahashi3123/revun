#!/bin/bash

. run_revun.h

LD=${DIR_LD}/filter_images.py
LOG=log_filter_images

# Modal analysis
python3 $LD 2>&1 | tee ${LOG}
