#!/bin/bash

. run_revun.h

LD_SEN=${DIR_LD}/sensing_aerodynamics.py
LOG_SEN=log_sensing

# Mode sensing from reconstruction data
python3 ${LD_SEN} 2>&1 | tee ${LOG_SEN}
