#!/bin/bash

FILE_INPUT=$1  # example-->hayabusa.%04d
FILE_OUTPUT=$2 # hayabusa
#num_frame=48
#num_frame2=96
num_frame=30
num_frame2=60
START_STEP=0
#FILE_NAME=${FILE_INPUT}

#ffmpeg -r 60 -i sma_anim.%04d.png -vcodec libx264 -pix_fmt yuv420p out.mp4
ffmpeg \
   -r $num_frame \
   -start_number ${START_STEP} \
   -i ${FILE_INPUT}.png \
   -vcodec libx264 \
   -pix_fmt yuv420p \
   -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
   -r $num_frame2 \
   ${FILE_OUTPUT}.mp4

