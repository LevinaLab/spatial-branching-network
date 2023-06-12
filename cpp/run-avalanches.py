#!/usr/bin/python

#Launcher to generate several avalanche programs in Granada's supercomputation PROTEUS

import os
import numpy as np


#Paths where the cpp lives (relative to this file)
file_to_compile = "lbnm.cpp"

#Where to save the binary
bin_folder = "bin/"
output_bin = bin_folder + "lbnm_aval.exe"

is_manhattan = False
is_mf = True

#Compiler flags
cpp_flags = "-std=c++11 -O3 -DMODE=AVAL "
if is_manhattan:
    cpp_flags += "-DLATTICE=MANHATTAN "
if is_mf:
    cpp_flags += "-DMF=TRUE"


#Make sure bin folder exists and compile 
os.system("mkdir -p " + bin_folder)
os.system("g++ {path} {flags} -o {binfile}".format(path=file_to_compile, flags=cpp_flags, binfile=output_bin))


#Define the parameters: system size, autoexcitation, neighbourhood, branching ratio 
L_list = [16]#, 32, 64, 128, 256, 512]#[8,32,128] 
k_list = [1]#[1,2,3]
p_list = [0.001, 0.1, 0.5]
ps = 0.5
naval = 10 #1e7

is_subcritical = False
use_biggest_size = False 
rewiring = False

#Subcritical: set all the dictionary with unit branching
if is_subcritical:
    branching = {}
    if not rewiring:
        for k in k_list:
            branching[k] = {}
            for L in L_list:
                branching[k][L] = 1.0
    else:
        for p in p_list:
            branching[p] = {}
            for L in L_list:
                branching[p][L] = 1.0
else:
    #In any other case, get the values we have determined
    if not is_manhattan:
        branching = {1: {8: 1.13, 16:1.1396, 32:1.11, 64:1.1152, 128:1.109, 256:1.109, 512:1.109},
                    2: {8: 1.10,  16:1.108,  32:1.062, 64:1.0628, 128:1.05},
                    3: {8: 1.08, 16:1.088, 32:1.045,64:1.0436, 128:1.03},
                    0.001: {8:1.13, 32:1.112, 128:1.104},
                    0.1:   {8:1.10, 32:1.07, 128: 1.0548},
                    0.5:   {8:1.08, 32:1.05, 128:1.022}
                    }

    else:
        branching = {1: {16:1.1644, 32:1.1422, 64:1.1284, 128:1.1248},
                    3: {16:1.1032, 32:1.0726, 64:1.0552, 128:1.048},
                    5: {16:1.09, 32:1.054, 64:1.0354, 128:1.0258}}


nfiles = {8:10, 16:10, 32:10, 64:10, 128:100, 256:200, 512:300}

#Set output path and ensure it exists
output_folder = "../results/aval-mf/"
os.system("mkdir -p " + output_folder)

#Set parameters and run simulation
if not rewiring:
    params = {"p_s":ps, "naval":naval, "p_rewire": 0.0} 

    if not is_mf:
        for i,k in enumerate(k_list):
            params["k"] = k

            for j,L in enumerate(L_list):

                #Get the critical branching for this k,size. Allows also to get the largest size for everyone
                if use_biggest_size:
                    m = branching[k][L_list[-1]]
                else:
                    m = branching[k][L]


                #The C++ program needs p_n instead of branching, compute it
                num_neighs = 4*k*(k+1) if not is_manhattan else 2*k*(k+1)
                pr = (m - ps) / num_neighs

                params["L"] = L
                params["p_r"] = pr
                repetitions = nfiles[L]
                for r in range(repetitions):
                    output_path = output_folder + "aval_{0}_{1}_part{2}".format(k, L, r)
                    os.system("./{binfile} {L} {p_s} {p_r} {k} {naval} {p_rewire} {output} {index}".format(**params, binfile=output_bin, output=output_path, index=r))
    else:
        pr = 0.393469 #Critical: 0.393469  
        k=1 #Needed by some formaters and so on, but does not matter
        params["p_r"] = pr
        params["k"] = k
        for j,L in enumerate(L_list):

            params["L"] = L
            repetitions = nfiles[L]
            for r in range(repetitions):
                output_path = output_folder + "aval_{0}_{1}_part{2}".format(k, L, r)
                os.system("./{binfile} {L} {p_s} {p_r} {k} {naval} {p_rewire} {output} {index}".format(**params, binfile=output_bin, output=output_path, index=r))

else:
    params = {"p_s":ps, "naval":naval, "k": 1} 
    k = 1
    for i,p in enumerate(p_list):
        for j,L in enumerate(L_list):

            #Get the critical branching for this k,size. Allows also to get the largest size for everyone
            if use_biggest_size:
                m = branching[p][L_list[-1]]
            else:
                m = branching[p][L]

            #The C++ program needs p_n instead of branching, compute it
            num_neighs = 4*k*(k+1) if not is_manhattan else 2*k*(k+1)
            pr = (m - ps) / num_neighs

            params["L"] = L
            params["p_r"] = pr
            repetitions = nfiles[L]
            for r in range(repetitions):
                output_path = output_folder + "aval_{0:.3f}_{1}_part{2}".format(p, L, r)
                os.system("./{binfile} {L} {p_s} {p_r} {k} {naval} {p_rewire} {output} {index}".format(**params, binfile=output_bin, output=output_path, index=r))