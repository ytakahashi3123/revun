#

FILE_INPUT=$1  # movie
FILE_OUTPUT=$2 # image_%04d
FORMAT=png

ffmpeg -i ${FILE_INPUT}.mp4 -vcodec ${FORMAT} ${FILE_OUTPUT}.${FORMAT}
