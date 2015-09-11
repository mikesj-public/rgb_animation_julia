# rgb_animation_julia

Project to create animations which have every 24-bit RGB colour in them once (or a set number of times each)

Run with command line julia rgb_animation.jl output_folder --args  
Details of args can be seen via julia rgb_animation.jl --help

This will create a new directory output_folder (in the same directory as rgb_animation), then create a series of images in the folder based on the RGB animation algorithm.  To create an mp4 (24 frames/s) in a unix system, you can cd to the folder and run 
ffmpeg -f image2 -r 24 -pattern_type glob -i '*.png' animation.mp4

Blog post about this program can be found at https://swarbrickjones.wordpress.com/2015/09/11/all-24-bit-rgb-colours-in-one-animation/
