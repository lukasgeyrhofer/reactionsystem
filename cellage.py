#!/usr/bin/env python

import numpy as np
import argparse
import math,sys

def age(x):
    return 1.*np.mean(x)

parser = argparse.ArgumentParser()
parser.add_argument("-n","--initialcellnumbers",type=int,default=1)
parser.add_argument("-G","--generations",type=int,default=10)
args = parser.parse_args()

pop = [np.zeros(2,dtype=np.int) for i in range(args.initialcellnumbers)]

for g in range(args.generations):
    ps = len(pop)
    mothercells = np.random.randint(low = 0, high = ps, size = ps)
    for i in mothercells:
        pop.append(np.array([pop[i][1],g],dtype=np.int))
        pop[i][1] = g

agehisto = np.zeros(shape = (args.generations,args.generations))
for age in pop:
    agehisto[age[0],age[1]] += 1

for x in range(args.generations):
    for y in range(args.generations):
        print x,y,agehisto[x,y]
    print
