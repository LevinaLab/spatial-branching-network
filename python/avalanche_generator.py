

import numpy as np
from lattice_activity import *

def av_gen(m,n,sigma,pext,ps,k,num_av,ptype,filename, prw=None, self_exciteP=1):
    """
    Computes a number num_av of avalanches and stores its size and duration into the file.

    Parameters
    ===========

    - m, n : int
        The lattice size
    - sigma : float
        Branching parameter
    - pext : float
        Spontaneous firing probability, between 0 and 1
    - ps : float
        Probability of self-excitation
    - k : int
        Connectivity radious
    - num_av : int
        Number of avalanches to simulate
    - ptype : string
        Should be either uniform or rewired_uniform
    - filename : string
        Where to store the result
    - prw : float
        Rewiring probability, only if ptype is rewired_uniform. Defaults to None.
    - self_excitedP : legacy code, should be set to 1
    """
         
    neigh_all = find_allneigh_strongSelf(n,m,k+1) 
   
    if ptype == 'uniform':
        p = compute_p_uniform(k,pext,ps,sigma,self_exciteP)
        print('p=',p)
    elif ptype == 'rewired_uniform':
        p = compute_p_uniform(k,pext,ps,sigma,self_exciteP)
        neigh_all = rewire_uniform(n, m, k, neigh_all, prw)
        print('p=',p)
    else:
        print("Wrong ptype")
        raise ValueError

    center = int(np.ceil(m/2)*n + np.ceil(n/2))

    with open(filename, 'w') as output:
        for i in range(num_av):
            #save_counter = save_counter+1
            s = np.zeros(n*m)
            s[center] = 1
            total_active = 1
            t = 0
            while(sum(s>0) > 0):
                #lattice_activity.append(s)
                s = np.array([update_network_states_ref0(k,m,n,s,p,neigh_all) for k in range(0,len(s))])     

                t = t + 1
                total_active = total_active + sum(s>0)

            output.write("{size} {duration}\n".format(size = total_active, duration=t))




def av_gen_coalescence(m,n,sigma,pext,ps,k,num_av,ptype, prw = None,self_exciteP=1):
    """
    Computes a number num_av of avalanches and stores its size and duration into the file.

    Parameters
    ===========

    - m, n : int
        The lattice size
    - sigma : float
        Branching parameter
    - pext : float
        Spontaneous firing probability, between 0 and 1
    - ps : float
        Probability of self-excitation
    - k : int
        Connectivity radious
    - num_av : int
        Number of avalanches to simulate
    - ptype : string
        Should be either uniform or rewired_uniform
    - filename : string
        Where to store the result
    - prw : float
        Rewiring probability, only if ptype is rewired_uniform. Defaults to None.
    - self_excitedP : legacy code, should be set to 1

    Returns
    =======
    - numActive : list of ndarray
        Each element of this list contains the number of active units in each timestep. Each element corresponds to a different avalanche. 
    - m_coal : list of ndarray
        Each element of this list contains the amount of coalescence at each timestep. Each element corresponds to a different avalanche. 
    - avalanche_size: ndarray
        Array with the size of each avalanche. It has num_av elements.
    """
    
    neigh_all = find_allneigh_strongSelf(n,m,k+1) 
   
    if ptype == 'uniform':
        p = compute_p_uniform(k,pext,ps,sigma,self_exciteP)
        print('p=',p)
    elif ptype == 'rewired_uniform':
        p = compute_p_uniform(k,pext,ps,sigma,self_exciteP)
        neigh_all = rewire_uniform(n, m, k, neigh_all, prw)
        print('p=',p)
    else:
        print("Wrong ptype")
        raise ValueError

    center = int(np.ceil(m/2)*n + np.ceil(n/2))
    avalanche_size = []
    avalanche_lifetime = []
    lattice_activity = []
    m_coal = []
    numActive_all = []
    m_real = []
    save_counter = 0
    for i in range(num_av):
        save_counter = save_counter+1
        s = np.zeros(n*m)
        s[center] = 1
        total_active = 1
        t = 0
        while(np.sum(s>0) > 0):
#             lattice_activity.append(s)
            s_old = s
            s = np.zeros(n*m)
            pot_s = np.zeros(n*m) # potential number of active units without coalesence
            for k in range(0,len(s_old)):
                s[k],pot_s[k] = update_network_states_ref0_conv(k,m,n,s_old,p,neigh_all)     
            
            #t = t + 1
            
            if np.sum(s) > 0:
                m_coal.append((np.sum(pot_s) - np.sum(s))/np.sum(s_old))
                numActive_all.append(np.sum(s_old))
                total_active = total_active + np.sum(s>0)
                
                
            if total_active > 10**15: # for killing large avalanches in case of supercriticality
                s = np.zeros(n*m)

        avalanche_size.append(total_active)
        #avalanche_lifetime.append(t)

    return numActive_all, m_coal, avalanche_size






def phase_diagram(m,n,sigma_vec,pext,ps,k,num_its,num_meas,relax_its,ptype,filename,prw=None, self_exciteP=1):
    """
    Computes a number num_av of avalanches and stores its size and duration into the file.

    Parameters
    ===========

    - m, n : int
        The lattice size
    - sigma : float
        Branching parameter
    - pext : float
        Spontaneous firing probability, between 0 and 1
    - ps : float
        Probability of self-excitation
    - k : int
        Connectivity radious
    - num_its : int
        Number of iterations to measure on
    - num_meas : int
        Measures are taken each num_meas steps
    - relax_its : string
        Leave the system relax to stationary state for relax_its iterations prior to start measurements
    - ptype : string
        Should be either uniform or rewired_uniform
    - filename : string
        Where to store the result
    - prw : float
        Rewiring probability, only if ptype is rewired_uniform. Defaults to None.
    - self_excitedP : legacy code, should be set to 1
    """
    

    neigh_all = find_allneigh_manhattan(n,m,k+1) 

    #Initialize vectors to store phase diagram results
    n_points = sigma_vec.size

    output = open(filename, "w")

    #Phase diagram!
    for j,sigma in enumerate(sigma_vec):
        if ptype == 'uniform':
            p = compute_p_uniform(k,pext,ps,sigma,self_exciteP, True)
        elif ptype == 'rewired_uniform':
            p = compute_p_uniform(k,pext,ps,sigma,self_exciteP)
            neigh_all = rewire_uniform(n, m, k, neigh_all, prw)
            print('p=',p)
        else:
            print("Wrong ptype")
            raise ValueError


        #Random initial conditions with density ~0.5
        s = np.random.randint(low=0, high=2, size=n*m)
        active = s.sum()

        #To store the activity at each time
        current_activity = np.zeros(num_its // num_meas)

        #Thermalization steps
        t = 0
        while (t < relax_its and active > 0):
            s = np.array([update_network_states_ref0(k,m,n,s,p,neigh_all) for k in range(0,len(s))])     
            active = s.sum()
            t += 1
        
        #Measurement steps

        if active > 0: 
            t = 0
            while (t < num_its and active > 0):
                #Update lattice
                s = np.array([update_network_states_ref0(k,m,n,s,p,neigh_all) for k in range(0,len(s))])     
                
                #Measure activity at this time
                active = s.sum()
                if t % num_meas == 0:
                    current_activity[t // num_meas] = active / (n*m)
                
                #Next time
                t += 1

            #If activity survived enough time, we are either supercritical or either approaching criticality
            #rho_s[j,0] = current_activity.mean()
            #rho_s[j,1] = current_activity.var() 
            output.write("{s} {rho} {v}\n".format(s=sigma, rho=current_activity.mean(), v=current_activity.var()))

        else:
            #If we fell during thermalization in absorbing state, then without any doubt...
            #rho_s[j,0] = 0.0
            #rho_s[j,1] = 0.0
            output.write("{s} {rho} {v}\n".format(s=sigma, rho=0.0, v=0.0))

    #np.save("phase_diagram", rho_s)
    output.close()


def phase_diagram_parallel(m,n,sigma,pext,ps,k,num_its,num_meas,relax_its,ptype,filename,prw=None,self_exciteP=1):
    """
    Computes a number num_av of avalanches and stores its size and duration into the file.

    Parameters
    ===========

    - m, n : int
        The lattice size
    - sigma : float
        Branching parameter
    - pext : float
        Spontaneous firing probability, between 0 and 1
    - ps : float
        Probability of self-excitation
    - k : int
        Connectivity radious
    - num_its : int
        Number of iterations to measure on
    - num_meas : int
        Measures are taken each num_meas steps
    - relax_its : string
        Leave the system relax to stationary state for relax_its iterations prior to start measurements
    - ptype : string
        Should be either uniform or rewired_uniform
    - filename : string
        Where to store the result
    - prw : float
        Rewiring probability, only if ptype is rewired_uniform. Defaults to None.
    - self_excitedP : legacy code, should be set to 1
    """
    
    neigh_all = find_allneigh_strongSelf(n,m,k+1) 
   
    #Initialize vectors to store phase diagram results
    rho_s = np.empty(2)
    
    if ptype == 'uniform':
        p = compute_p_uniform(k,pext,ps,sigma,self_exciteP, True)
        #print('p=',p)
    elif ptype == 'rewired_uniform':
        p = compute_p_uniform(k,pext,ps,sigma,self_exciteP)
        neigh_all = rewire_uniform(n, m, k, neigh_all, prw)
        print('p=',p)
    else:
        print("Wrong ptype")
        raise ValueError


    #Random initial conditions with density ~0.5
    s = np.random.randint(low=0, high=2, size=n*m)
    active = s.sum()

    #To store the activity at each time
    #Zero so if simulation is finished by reaching absorbing, averages make sense
    current_activity = np.zeros(num_its // num_meas)

    #Thermalization steps
    t = 0
    while (t < relax_its and active > 0):
        s = np.array([update_network_states_ref0(k,m,n,s,p,neigh_all) for k in range(0,len(s))])     
        
        active = s.sum()
        t += 1
    

    #Measurement steps

    if active > 0: 
        t = 0
        while (t < num_its and active > 0):
            #Update lattice
            s = np.array([update_network_states_ref0(k,m,n,s,p,neigh_all) for k in range(0,len(s))])     
            
            #Measure activity at this time
            active = s.sum()
            if t % num_meas == 0:
                current_activity[t // num_meas] = active / (n*m)
            
            #Next time
            t += 1

        #If activity survived enough time, we are either supercritical, or subcritical approaching criticality
        rho_s[0] = current_activity.mean()
        rho_s[1] = current_activity.var() 

    else:
        #If we fell during thermalization in absorbing state, then without any doubt...
        rho_s[0] = 0.0
        rho_s[1] = 0.0

    #Write result
    with open(filename, 'w') as the_file:
        the_file.write("{sigma} {rho} {varrho}".format(sigma=sigma, rho=rho_s[0], varrho=rho_s[1]))
 