#!/bin/bash

DIR_LD=../src
ncpu=8

export OMP_NUM_THREADS=$ncpu

PYTHON=python3.9

# note
#tar zcvf result_dmd.tar.gz result_dmd 
