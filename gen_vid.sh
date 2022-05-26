#!/usr/bin/env bash
# generate a video "output.mp4" using the captured images in "output" folder
ffmpeg -i 'output/image%06d.png' -r 10 -pix_fmt yuv420p inverse_dynamics.mp4
