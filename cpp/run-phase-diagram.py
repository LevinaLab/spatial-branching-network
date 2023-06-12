#!/usr/bin/python

#Launcher to generate several avalanche programs in Granada's supercomputation PROTEUS

import os

#Paths where the cpp lives (relative to this file)
file_to_compile = "lbnm.cpp"

#Where to save the binary
bin_folder = "bin/"
output_bin = bin_folder + "lbnm-diagram.exe"

#Compiler flags
cpp_flags = "-std=c++11 -O3 -DMODE=DIAGRAM" 

#Make sure bin folder exists and compile 
os.system("mkdir -p " + bin_folder)
os.system("g++ {path} {flags} -o {binfile}".format(path=file_to_compile, flags=cpp_flags, binfile=output_bin))

#Parameters to obtain phase diagram
Ls = [8] # 32, 64, 128] 
ps = 0.5
ks = [1,3]
m0,mf = 1.0, 1.2
number_m = 10

#To divide work among different programs 
number_programs = 1 
delta_program = (mf-m0)/number_programs

params = {"p_s":ps, "p_rewire":0.1}

#Set output path and ensure it exists
output_folder = "../results/phase-diagram-cpp/"
os.system("mkdir -p " + output_folder)

#Divide the work among different programs, then launch them all
for k in ks:
    params["k"] = k
    for L in Ls:
        params["L"] = L
        for j in range(number_programs):
            m0_program = m0 + delta_program*j 
            mf_program = m0 + delta_program*(j+1)

            params["m"] = [m0_program, mf_program, number_m]

            output_path = output_folder + "diagram_k{0}_L{1}_part{2}".format(k, L, j)

            #Execute the program
            os.system("./{binfile} {L} {p_s} {k} {m[0]} {m[1]} {m[2]} {p_rewire} {output}".format(**params, binfile=output_bin, output=output_path))