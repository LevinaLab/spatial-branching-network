#!/usr/bin/env python

import avalanche_generator as av
import numpy as np
import sys
from os import system

if len(sys.argv) == 5:
    L = int(sys.argv[1])
    sigma_ = float(sys.argv[2])
    k = int(sys.argv[3])
    index = int(sys.argv[4])
else:
    print("Bad number of arguments where given. Call this as")
    print("'python program.py L sigma k folder index'")
    print("L: lenght of side, N=L*L")
    print("sigma: branching parameter")
    print("k: neighbourhood (int > 1)")
    print("folder: the folder where the output will be saved (relative path!)")
    print("index: a number to avoid all programs writing to the same file")
    exit()


# lattice size
m = L
n = L

# branching ratio
sigma = sigma_

# external input, put to 0 for avalanches
pext = 0

# self-excitaion (p_s in thesis)
ps = 0.5

# number of avalanches to generate, 
num_av = 1000
counter = num_av

# connectivity structure (put 'uniform' for structured lattice)
ptype = 'uniform'

# some default settings
self_exciteP = 1
self_excite_neigh = 1


#where to save the results from this simulation
system("mkdir ../results/avalanches/")
filename = "../results/avalanches/av-{size}-{index}".format(size=n*m, index=index)

av.av_gen(m,n,sigma,pext,ps,k,num_av,ptype,self_exciteP,self_excite_neigh,filename)


