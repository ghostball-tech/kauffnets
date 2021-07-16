# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 20:10:52 2021

Metabolic Stability and Epigenesis in Randomly Constructed Genetic Nets by Stuart 
Kauffman, Journal of Theoretical Biology (1969) vol 22, 437-467


Random (uniform) values of p and K for each iteration.

Take mean over all agent states at T = 20.

save results as a data frame (panda) and write to csv to be assimilated in R

@author: bruce
"""



import numpy as np
import pandas as pd 

import time

start_time = time.time()

PLOT_GATE = True

######################################
######################################
# fixed population parameters
######################################

N = 1000 # number of "genes" or agents or nodes
T = 40#time steps

# free parameters ######################################

# two parameters: k and p
# chosen uniformly at random from following intervals
K = [1, 8] # [K_min, K_max]
# p uniform on [0,1)

# expected number of 1's in initial conditions
# x_0 random, defined inside parameter loop

# number of iterations / samples
nIter = 30
#####################################

#Got seed ? 
#np.random.seed(718281828)

#####################################################
# ABM function
######################################################

def kaufNetABM(N,x_0,K,p,T):
    stateMean = np.array([])
    #choose influencers for each node
    edgeList = []
    for agent in range(N):
        k = K
        #k = np.random.poisson(K-1)+1
        edgeList = edgeList + [np.random.choice(N,k,replace=False)]
        
    #state vector
    X = np.zeros((T+1,N))
    #set initial conditions
    initialConditions = np.random.binomial(1,x_0,N)
    X[0,] = initialConditions   
        
    stateMean = np.append(stateMean,X[0,].mean())
        
    #random logic for each agent
    # each agent assigned a random map, agentMap (dictionary)
    # maps are indexed in a (list), mapList
    # there are 2**K possible inputs, and a p-coin is flipped to determine the value at each input
    mapList = []
    for agent in range(N):
        agentMap = {}
        for key in range(2**K):
            agentMap[key] = np.random.binomial(1,p)
        mapList = mapList + [agentMap]
        
    #dynamics
    for t in range(T):
        for agent in range(N):
       #base 2 representation of inputs, to match encodeing in dictionary 
            mapKey = 0
            for k in range(K):     
                mapKey = mapKey + (2**k)*X[t,edgeList[agent][k]]
            X[t+1,agent] = mapList[agent][mapKey]
        
        #record stats
        stateMean = np.append(stateMean, X[t+1,].mean())
    return(stateMean[T])

######################################################### end ABM function



# Data ####################
output = []
###############################


################################
# Iterating ABM over parameters
################################


for iter in range(nIter):
    x_0 = np.random.rand()
    k = np.random.randint(K[0],K[1])
    p = np.random.rand()
    stateMean = kaufNetABM(N,x_0,k,p,T)
    output.append([p,k,stateMean])

            
print("---  ABM sims: %s seconds ---" % (time.time() - start_time))


#############################################################
#%% write to csv
output = np.array(output)
output_df = pd.DataFrame(data = output, columns = ['p','k','state_mean'])

#output_df.to_csv('kauffman_abm_output.csv',index = False)


