import numpy as np
import scipy.stats as stats
import sys

num_max = 10**10
num_avalanches = 50
ps = 0.5
L = 8
raw_sigma = 1.0
k = 1

n = (2*k+1)**2 - 1 # number of connected neighbors


path_coal = '../results/coalescence/'
av_path = '../results/coalescence/'

num_neigh = 1
pext = 0
num_av = 50
ptype = 'uniform'
self_exciteP = 1
self_excite_neigh = 1
counter = num_av


def bp_avSize_adaptive(conv_list, raw_sigma, n, num_max, num_avalanches):
    """
    Function to simulate an adaptive BR process and return the avalanche size distribution
    num_max: size of the largest avalanche (equal to the network size)
    num_avalanches: total number of avlanches
    raw_sigma: raw branching parameter in the network
    n: number of neighbors in the network
    """
    avSize = []
    for t in range(num_avalanches):
        num_descend = 1
        track_descend = []
        num_people = 1
        while num_people <= num_max :

            track_descend.append(num_descend)
            descend = 0
            
            #find the adaptive sigma
            satLevel = len(conv_list)
            if num_descend < satLevel:
                sigma = raw_sigma - conv_list[num_descend-1]
                p = sigma / n
                pk = []
                for k in range(n+1):
                    pk.append(stats.binom.pmf(k,n,p))
            else:
                sigma = raw_sigma - conv_list[-1]
                p = sigma / n
                pk = []
                for k in range(n+1):
                    pk.append(stats.binom.pmf(k,n,p))
            
            pk_cum = np.cumsum(pk)                            
            #generate decendents
            for i in range(num_descend):
                d = np.random.rand()
                descend = descend + np.min(np.where(pk_cum >= d)[0])


            num_descend = descend
            num_people = num_people + num_descend

            if num_descend == 0:
                break

        if num_people > num_max:
            num_people = num_max + 1

        avSize.append(num_people)
    return avSize   



[avg_mcoal_vsAct, std_mcoal_vsAct, act_uniq, nums, mean_br, st_br, avs] = np.load(path_coal + 'mcoal_' + ptype +'_k'+str(k) +'_m'+str(raw_sigma)+ '_ps' + str(ps) + ".npy",  allow_pickle=True)
avs = 0 # empty not needed avsize

coal_list = np.zeros((int(np.max(act_uniq))))
for i in range(len(act_uniq)):
    coal_list[int(act_uniq[i])-1] = avg_mcoal_vsAct[i]         
    
    
for i in range(1,len(coal_list)):
    if coal_list[i] == 0:
        coal_list[i] = coal_list[i-1]+0
        
        

avSize = bp_avSize_adaptive(coal_list, raw_sigma, n, num_max, num_avalanches)

np.save(av_path +'avSize_adaptBR_'+ ptype +'_k'\
        +str(k) +'_sig'+str(raw_sigma)+ '_ps' + str(ps)+'_size'+str(L),avSize)


