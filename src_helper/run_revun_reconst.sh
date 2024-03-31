#!/bin/bash

. run_revun.h

LD=${DIR_LD}/reconstruct_data.py
LOG=log_reconstruct

# Reconstruct dataset based on modal analysis results
python3 ${LD} 2>&1 | tee ${LOG}
