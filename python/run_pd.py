#!/usr/bin/env python

import avalanche_generator as av
import numpy as np
import sys
from os import system

# lattice size
m = 64 
n = 64 

# branching ratio
sigma_vec = np.linspace(1.0, 1.3, 2)

# external input, put to 0 for avalanches
pext = 0

# self-excitaion (p_s in thesis)
ps = 0.5

# neighborhood radius (k in thesis)
k = 1

# number of iterations
num_its = 300
num_meas = 10
relax_its = 70

# connectivity structure (put 'uniform' for structured lattice)
ptype = 'uniform'

# some default settings
self_exciteP = 1
self_excite_neigh = 1

# dirctories to save avlanches' size and duration
system("mkdir ../results/phase-diagram/")
filename = '../results/phase-diagram/phase-diagram'


#Single phase diagram
av.phase_diagram(m,n,sigma_vec,pext,ps,k,num_its,num_meas,relax_its,ptype,self_exciteP,self_excite_neigh,filename)



