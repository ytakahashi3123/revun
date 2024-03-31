#!/bin/bash

. run_revun.h

LD_SEN=${DIR_LD}/extract_feature.py
LOG_SEN=log_extract

# Mode sensing from reconstruction data
$PYTHON ${LD_SEN} 2>&1 | tee ${LOG_SEN}
