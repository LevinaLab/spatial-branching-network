#!/usr/bin/env python

import avalanche_generator as av
import numpy as np
import sys
from os import system

m = 8
n = m
sigma = 1.
k = 1
# sigma = 1.054
# num_neigh = 2
# sigma = 1.0354
# num_neigh = 3
pext = 0
ps = 0.5
num_av = 10
ptype = 'uniform'
self_exciteP = 1
self_excite_neigh = 1
counter = num_av



system('mkdir ../results/coalescence/')
direc = '../results/coalescence/'

act, mcoal, avsize = av.av_gen_coalescence(m,n,sigma,pext,ps,k,num_av,ptype,self_exciteP,self_excite_neigh,counter)

# compute average and std coalesence conditioned on number of active units    
act = np.asarray(act)
mcoal = np.asarray(mcoal)
act_uniq = np.unique(act)

avg_mcoal_vsAct = []
std_mcoal_vsAct = []
num_mcoal_vsAct = []
for i in range(len(act_uniq)):
    avg_mcoal_vsAct.append(np.mean(mcoal[act==act_uniq[i]]))
    std_mcoal_vsAct.append(np.std(mcoal[act==act_uniq[i]]))
    num_mcoal_vsAct.append(len(mcoal[act==act_uniq[i]]))
    
avg_mcoal_vsAct = np.asarray(avg_mcoal_vsAct)
std_mcoal_vsAct = np.asarray(std_mcoal_vsAct)


np.save(direc + 'mcoal_' + ptype +'_k'+str(k) +'_m'+str(sigma)+ '_ps' + str(ps),\
     [avg_mcoal_vsAct, std_mcoal_vsAct, act_uniq, num_mcoal_vsAct, np.mean(mcoal), np.std(mcoal), avsize])