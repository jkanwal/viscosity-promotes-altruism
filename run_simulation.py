#Load libraries
import numpy as np
import random
import copy
import sys 
import time
from datetime import datetime

#Set parameter values
n = int(sys.argv[1]) #breeders per patch
N = int(50000/n) #number of patches (set so that all simulation runs have equal total breeder population, regardless of breeders per patch)
c = float(sys.argv[2]) #cost of dispersal
k = 100 #offspring per preeder
m = 0.01 #mutation rate

z_star = 2/(1 + 2*c*n + np.sqrt(1+4*n*(n-1)*(c**2)))
b1 = z_star
b2 = z_star

#Generate filename using parameters
now = datetime.now()
timestamp = str(now.date()) + '-' + str(now.time().hour) + '.' + str(now.time().minute)
run_name = "n_" + str(n) + "_c_" + str(c) + "_" + timestamp
filename = 'scratch/' + run_name + '.csv'
disp_filename = 'scratch/' + run_name + '_disp.csv'
ind_filename = 'scratch/' + run_name + '_ind.csv'

# Run simulation for G generations
G = 20000

#Clone function (with mutation)
#Dispersal genes mutate incrementally, cooperation gene mutates to a random value between 0 and 1
def Clone(ind, m):
    offspring = copy.deepcopy(ind)
    #gene b1
    r0 = random.random()
    if r0 <= m/2:
        offspring[0] = ind[0] - 0.05
    elif r0 <= m:
        offspring[0] = ind[0] + 0.05
    #gene b2  
    r1 = random.random()
    if r1 <= m/2:
        offspring[1] = ind[1] - 0.05
    elif r1 <= m:
        offspring[1] = ind[1] + 0.05
    #gene c
    r2 = random.random()
    if r2 <= m:
        offspring[2] = random.random()
    return offspring

#Reproduction function
#Replace each individual in previous gen with k clones, with some chance of mutation each time
def Reproduce(Gen, N, n, k, m):
    new = np.empty((N,n*k,3))
    n = 0
    for i in Gen:
        o_list = []
        for j in i:
            for kk in range(k):
                offspring = Clone(j, m)
                o_list.append(offspring)
        new[n] = o_list
        n += 1
    return new

#Survival function
#This function takes a patch as input, and selects the offspring that survive
def Survival(patch, n, k):
    survivors = []
    avg_c = np.sum(patch, axis=0)[2]/(n*k)
    #Add some extrinsic patch variation - each patch has extrinsic survival rate between 0.5 and 1
    SE = random.uniform(0.5,1)
    for i in patch:
        S = 0.5 + 0.5*(1-i[2])*avg_c
        r = random.random()
        if r <= S:
            r2 = random.random()
            if r2 <= SE:
                survivors.append(i.tolist())
    return survivors

#Dispersal functions
def DisperseOut(patch, c, n, k):
    P = len(patch)/n #per-breeder patch productivity
    d_count = 0
    dispersers = []
    remainers = []
    for i in patch:
        D = ((P-0.5*k)*(i[1]-i[0])/(0.5*k)) + i[0]
        r = random.random()
        if r <= D:
            d_count += 1
            if r <= D*(1-c):
                dispersers.append(i)
        else:
            remainers.append(i)
    return [d_count, dispersers, remainers]

def DisperseIn(dispersers, patches):
    pn = len(patches)-1
    for i in dispersers:
        r = random.randint(0,pn)
        patches[r].append(i)
    return patches

#Select new breeders
def SelectBreeders(patches, N, n):
    newGen = np.empty((N,n,3))
    for i in range(N):
        if len(patches[i]) == 0:
            for j in range(n):
                newGen[i,j] = [b1, b2, random.random()]               
        elif len(patches[i]) < n:
            for j in range(n):
                rand = random.randint(0,len(patches[i])-1)
                newGen[i,j] = patches[i][rand]
        else:
            np.random.shuffle(patches[i])
            newGen[i] = patches[i][0:n]
    return newGen

#Generation advancer:
def Step(Gen, N, n, k, m, c, step):
    offspr = Reproduce(Gen, N, n, k, m)
    patches = []
    for patch in offspr:
        patches.append(Survival(patch, n, k))
    dispersers = []
    for i in range(len(patches)):
        disp = DisperseOut(patches[i], c, n, k)
        if step >= G-1000:     #only write disp data for final 1000 gens
            with open(disp_filename, 'a+') as g:
                g.write(str(step) + ', ' + str(len(patches[i])) + ', ' + str(disp[0]) + '\n')
        dispersers += disp[1]
        patches[i] = disp[2]        
    patches = DisperseIn(dispersers, patches)
    newGen = SelectBreeders(patches, N, n)
    return newGen

#Function to write data to a file - write medians, then means
def WriteAvg(Gen, filename):
    pop_b1 = Gen[:,:,0].flatten()
    pop_b2 = Gen[:,:,1].flatten()
    pop_c = Gen[:,:,2].flatten()
    with open(filename, 'a+') as f:
        f.write(str(np.median(pop_b1)) + ',' + str(np.median(pop_b2)) + ',' +  str(np.median(pop_c)) + ',' + str(np.mean(pop_b1)) + ',' + str(np.mean(pop_b2)) + ',' +  str(np.mean(pop_c)) + '\n')

#Function to write EVERY INDIVIDUAL'S data to a file for a single generation
def WriteInd(Gen, filename):
    with open(filename, 'a+') as f:
        patch = 0
        for i in Gen0:
            for j in i:
                f.write( str(patch) + ',' + str(j).strip(']').strip('[').replace(' ', ',').strip(',,') + '\n')
            patch += 1

#Generate first generation
#Each individual is an array of 3 gene values. Each patch is an array of these arrays. A generation is an array of THESE arrays.
Gen0 = np.empty((N,n,3))
for i in range(N):
    for j in range(n):
        Gen0[i,j] = [b1, b2, random.random()]

#Run simulation
print("start time: " + str(datetime.now().time()))
start_time = time.time()

WriteAvg(Gen0, filename)
step = 0
Gen = Gen0
for i in range(G):
    Gen_new = Step(Gen, N, n, k, m, c, step)
    Gen = Gen_new
    WriteAvg(Gen, filename)
    step += 1
WriteInd(Gen, ind_filename) #write individual data for final gen only

print("time taken: %s minutes" % ((time.time() - start_time)/60))   