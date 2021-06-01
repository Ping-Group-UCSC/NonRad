#!/usr/bin/env bash
module load numpy
dir="/home/tjsmart/Programs/Ping-Group/NonRad/code"
python2 -u $dir/calc_auto.py | tee std.out
