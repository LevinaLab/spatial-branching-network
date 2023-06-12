#!/usr/bin/python

import os

#Ask the OS to return absolute path to this file
path_2_this = os.path.dirname(__file__)
path_2_this = path_2_this if path_2_this != "" else "."

#Paths where the cpp lives (relative to this file)
cpp = path_2_this + "/../../cpp/"
file_to_compile = cpp + "ibnm.cpp"

#Where to save the binary
bin_folder = cpp + "bin/"
output_bin = bin_folder + "ibnm.exe"

#Compiler flags
cpp_flags = "-std=c++11 -O3 -DMODE=DIAGRAM"

#Make sure bin folder exists and compile
os.system("mkdir -p " + bin_folder)
os.system("g++ {path} {flags} -o {binfile}".format(path=file_to_compile, flags=cpp_flags, binfile=output_bin))

#Define the parameters: system size, autoexcitation, neighbourhood, and branching 
L = 64
ps = 0.5
k = 1
m0,mf = 1.0, 1.3
number_m = 30

#Set output path and ensure it exists
output_folder = path_2_this + "/../../data/phase_diagram/"
os.system("mkdir -p " + output_folder)
output_path = output_folder + "diagram"

#Execute the program
params = {"L":L, "p_s":ps, "k":k, "m":[m0, mf, number_m]}
os.system("./{binfile} {L} {p_s} {k} {m[0]} {m[1]} {m[2]} {output}".format(**params, binfile=output_bin, output=output_path))
