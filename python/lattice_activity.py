
import numpy as np

def onedim_to_twodim(k,m,n):
    i = k // m + 1 - 1
    j = k % m 
    return i,j


def twodim_to_onedim (i,j,m,n):
    i = i + 1
    j = j +1
    k = (i-1) * n + j -1 
    return k


def find_nth_neigh_general(k,m,n,nth):
    i,j = onedim_to_twodim(k,m,n)
    
    i_up_all = []
    i_down_all = []
    for ct in (np.arange(nth)+1):
        i_up = int(i-ct >= 0) * (i-ct) + (m-(ct-i)) * int(i-ct < 0)
        i_down = int(i+ct <= m-1) * (i+ct) + (ct - ((m-1)-i)-1) * int(i+ct > m-1)
        i_up_all.append(i_up)
        i_down_all.append(i_down)
        
    j_left_all = []
    j_right_all = []
    for ct in (np.arange(nth)+1):
        j_left = int(j-ct >= 0) * (j-ct) + (n-(ct-j)) * int(j-ct < 0)
        j_right = int(j+ct <= n-1) * (j+ct) + (ct - ((n-1)-j)-1) * int(j+ct > n-1)
        j_left_all.append(j_left)
        j_right_all.append(j_right)
        
    x = [i_up_all[-1]]*(2*nth+1)
    y = [i_down_all[-1]]*(2*nth+1)
    z = i_up_all[:-1] + [i] + i_down_all[:-1] 
    NB_i = np.array(x + y + z + z)

    
    xx = [j_right_all[-1]]*(2*nth-1)
    yy = [j_left_all[-1]]*(2*nth-1)
    zz = j_left_all + [j] + j_right_all 
    NB_j = np.array(zz + zz + xx + yy)
    NB = twodim_to_onedim (NB_i,NB_j,m,n)

    return NB

#Manhattan distance with periodic initial conditions
def manh_distance(x1,y1,x2,y2,m,n):
    dx = np.minimum(np.abs(x2 - x1) , n - np.abs(x2 - x1))
    dy = np.minimum(np.abs(y2 - y1) , m - np.abs(y2 - y1))
    return dx + dy

def find_allneigh_manhattan(n,m,k):
    # n and m are lattice dimensions
    # num_neigh is the maximum number of nearest neighbors(first,second,etc.)
    k = k + 1 #adding each cell as the strongest neighbor to itself
    num_cell = n*m
    neigh_all = []
    for i in range(num_cell):
        temp_NB = []
        for j in range(k):
            if j == 0:
                NB = np.array([i])
            else:
                NB = find_nth_neigh_general(i,m,n,j)
            temp_NB.append(NB)

        #Eliminate from the neighbourhood those which are further away than a Manhattan distance
        #(Absolutely inefficient but allows me to do this without changing the format of the neighbourhood list)
        temp_NB_flatten = np.concatenate(temp_NB).ravel() #convert to 1D list
        distance_to_element = np.empty(temp_NB_flatten.size) 

        #Manhattan distance to each one of the neighbours
        x1,y1 = onedim_to_twodim(i,m,n)
        for j,nb in enumerate(temp_NB_flatten):
            x2,y2 = onedim_to_twodim(nb, m, n)
            d = manh_distance(x1,y1,x2,y2,m,n)
            distance_to_element[j] = d
        
        #Re-arrange in the same format via masking the array
        temp_NB = [] #Recycling variables is a bad practice, but I was a bit lazy :)
        for j in range(k):
            temp_NB.append(temp_NB_flatten[distance_to_element == j])

        neigh_all.append(temp_NB)
    return neigh_all

def find_allneigh(n,m,k):
    # n and m are lattice dimensions
    # num_neigh is the maximum number of nearest neighbors(first,second,etc.)
    num_cell = n*m
    neigh_all = []
    for i in range(num_cell):
        temp_NB = []
        for j in range(k):
            NB = find_nth_neigh_general(i,m,n,j+1)
            temp_NB.append(NB)
        neigh_all.append(temp_NB)
    return neigh_all

def find_allneigh_strongSelf(n,m,k): # STRONG SELF_CONNECTIVITY
    # n and m are lattice dimensions
    # num_neigh is the maximum number of nearest neighbors(first,second,etc.)
    k = k + 1 #adding each cell as the strongest neighbor to itself
    num_cell = n*m
    neigh_all = []
    for i in range(num_cell):
        temp_NB = []
        for j in range(k):
            if j == 0:
                NB = np.array([i])
            else:
                NB = find_nth_neigh_general(i,m,n,j)
            temp_NB.append(NB)
        neigh_all.append(temp_NB)
    return neigh_all


def rewire_uniform(n, m, k, neigh_all, prw):
    # rewire the connections in structured network with probability pr
    n_units = m*n # total number of units
    for i in range(n_units):
        for j in range(1, k+1): # go through each radius, exclude self-excite
            n_neighbors = 8*j # find n_neighbors in that radius
            for k in range(n_neighbors):
                d = np.random.rand()
                if d < prw: # if probability of rewiring replace with a random unit
                    neigh_all[i][j][k] = random.randint(0, n_units-1)
    return neigh_all


def compute_p_hyp(k,pext,ps,sigma,self_excite):
    # Finding normalization constant
    if self_excite:
        c = sigma/(k*8+1)
    else:
        c = sigma/(k*8)
    # Computing p vector
    p = [pext]
    for i in range(1,k+1):
        if k==1:
            if self_excite:
                pi = sigma/(i*8+1)
            else:
                pi = sigma/(i*8)
        else:
            pi = c * (1/i)
        p.append(pi)
    return p

def compute_p_uniform(k,pext,ps,sigma,self_excite, is_manhattan=False):
    # Compute total number of neighbors
    if not is_manhattan:
        total_num_nei = (8 * k * (k+1))/2
    else:
        total_num_nei = (4 * k * (k+1))/2
    
    # Computing p vector
    p = [pext]
    if self_excite:
        sigma = sigma - ps
        p.append(ps)
        
    for i in range(1,k+1):
        pi = sigma/total_num_nei
        p.append(pi)
    return p
        
def update_network_state(k,m,n,s,s_ref,p,neigh_all,ref):  
    p_notActive = 1
    for i in range(len(p)):
        if i==0:
            p_notActive = p_notActive * (1-p[0])
        else:             
            NB = neigh_all[k][i-1]
            active_NB = sum(s[NB]>0)
            p_notActive = p_notActive * (1-p[i])**active_NB

    d = np.random.rand()
    s_new = int(s_ref[k] == 0) * int(d < (1-p_notActive)) + int(s_ref[k]>0) * (s[k])
    return s_new    


def update_network_states_ref0 (k,m,n,s,p,neigh_all):
    p_notActive = 1
    #print("index " + str(k) + " " + str(s[k-1:k+2]))
    for i in range(len(p)):
        if i==0:
            p_notActive = p_notActive * (1-p[0])
            #print(p[0], p_notActive)
        else:             
            NB = neigh_all[k][i-1]  
            active_NB = sum(s[NB]>0)
            p_notActive = p_notActive * (1-p[i])**active_NB
            #print(NB)
            #print(i,p[i],active_NB,p_notActive)
    
    d = np.random.rand()
    s_new =  int(d < (1-p_notActive))
    #print(s_new)
    #print("")
    return s_new


def update_network_states_ref0_conv (k,m,n,s,p,neigh_all):
    pot_s = 0
    #p_notActive = 1
    for i in range(len(p)):
        if i==0:
            #p_notActive = p_notActive * (1-p[0])
            d = np.random.rand()
            pot_s = pot_s + int(d < p[0])
        else:             
            NB = neigh_all[k][i-1]  
            active_NB = sum(s[NB]>0)
            #p_notActive = p_notActive * (1-p[i])**active_NB

            for j  in range(active_NB):
                d = np.random.rand()
                pot_s = pot_s + int(d < p[i])
    
    s_new =  int(pot_s>0)
    #d = np.random.rand()
    #s_new = (d < (1 - p_notActive))
    return s_new, pot_s



