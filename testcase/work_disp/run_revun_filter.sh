#!/bin/bash

. run_revun.h

LD=${DIR_LD}/filter_images.py
LOG=log_filter_images

# Modal analysis
$PYTHON $LD 2>&1 | tee ${LOG}
