#!/usr/bin/env bash

for img in ../../demonstrator/output/*.h5; do
    echo $img
    ./PlotSedov.py -i $img -o output/ -p 1 -a
done

ffmpeg -framerate 30 -pattern_type glob -i "output/*.png" -c:v libx264 -pix_fmt yuv420p sedov.mp4
