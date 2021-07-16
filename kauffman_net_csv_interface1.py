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

######################################
######################################
# fixed population parameters
######################################

N = 1000 # number of "genes" or agents or nodes
T = 40#time steps

# free parameters ######################################

# two parameters: k and p
# chosen uniformly at random from following intervals
K = [1.5, 6.5] # [K_min, K_max]
P = [0.01,0.5]

# expected number of 1's in initial conditions
# x_0 random, defined inside parameter loop

# number of iterations / samples
nIter = 50
nSample = 4
#####################################

#Got seed ? 
#np.random.seed(718281828)


#%% Helper functions
def HamDist(x,y):
    return(sum((x-y)**2))

#####################################################
#%% ABM function
######################################################

def kaufNetABM(N,x_0,K,p,T):

    #choose influencers for each node
    #Construct Graph
    edgeList = []
    edge_prob = K / N
    degreeList = []
    for agent in range(N):
        # choose how many neighbors each agent will have
        agent_degree = np.random.binomial(N,edge_prob)
        degreeList = degreeList + [agent_degree]
        # choose exactly which neighbors they are
        edgeList = edgeList + [np.random.choice(N,agent_degree,replace=False)]
        
    #state vector
    X = np.zeros((T+1,N))
    #set initial conditions
    initialConditions = np.random.binomial(1,x_0,N)
    X[0,] = initialConditions 
    #Perturbation
    # random reassignment of each element of X, iid p = 0.05 (Np = expected number of perturbed sites)
    Xb = X.copy()
    Xb[0,] = (Xb[0,] + np.random.binomial(1,0.01,N))%2
        
    #random logic for each agent
    # each agent assigned a random map, agentMap (dictionary)
    # maps are indexed in a (list), mapList
    # there are 2**K possible inputs, and a p-coin is flipped to determine the value at each input
    # here K is the (in)degree of each agent, degreeList[agent]
    mapList = []
    for agent in range(N):
        agentMap = {}
        coinflipArray = np.random.binomial(1,p,2**degreeList[agent])
        for key in range(2**degreeList[agent]):
            agentMap[key] = coinflipArray[key]
        mapList = mapList + [agentMap]
        
    #dynamics
    for t in range(T):
        for agent in range(N):
       #base 2 representation of inputs, to match encodeing in dictionary 
            mapKey = 0
            mapKeyB = 0
            for k in range(degreeList[agent]):     
                mapKey = mapKey + (2**k)*X[t,edgeList[agent][k]]
                mapKeyB = mapKeyB + (2**k)*Xb[t,edgeList[agent][k]]
            X[t+1,agent] = mapList[agent][mapKey]
            Xb[t+1,agent] = mapList[agent][mapKeyB]
        #record stats
    stateMean = X[t+1,].mean()
    hamDist = HamDist(X[t+1,],Xb[t+1,]) / N
    return([stateMean,hamDist])


#%%
# Data ####################
output = []
###############################


################################
# Iterating ABM over parameters
################################
start_time = time.time()

for iter in range(nIter):
    x_0 = np.random.rand()
    k = np.random.uniform(low = K[0],high = K[1])
    p = np.random.uniform(low = P[0], high = P[1])
    for exper in range(nSample):
        stateMean, hamDist = kaufNetABM(N,x_0,k,p,T)
        output.append([p,k,stateMean, hamDist])

            
print("---  ABM sims: %s seconds ---" % (time.time() - start_time))


#############################################################
#%% write to csv
output = np.array(output)
output_df = pd.DataFrame(data = output, columns = ['p','k','state_mean','hamming_dist'])

#example 2: K = [1.5, 6.5]  P = [0.01,0.5]
#output_df.to_csv('kauffman_abm_example2.csv',index = False)


