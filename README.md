# rgb_animation_julia

Run with command line julia rgb_animation.jl output_folder --args

details of args can be seen via julia rgb_animation.jl --help

This will create a new directory output_folder (in the same directory as rgb_animation), then create a series of images based on the RGB animation algorithm.  To create an mp4 (24 frames/s) in a unix system, you can cd to the folder and run: --  

ffmpeg -f image2 -r 24 -pattern_type glob -i '*.png' animation.mp4

Blog post about this program can be found at ::TODO::
