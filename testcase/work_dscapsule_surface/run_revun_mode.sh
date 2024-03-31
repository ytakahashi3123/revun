#!/bin/bash

. run_revun.h

LD=${DIR_LD}/mode_decomposition.py
LOG=log_mode

# Modal analysis
$PYTHON $LD 2>&1 | tee ${LOG}
