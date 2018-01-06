#!/usr/bin/env python

import numpy as np
import argparse
import sys,math


parser = argparse.ArgumentParser()
parser.add_argument("-N","--popsize_final",default=100000,type=int)
parser.add_argument("-n","--popsize_initial",default=10,type=int)
parser.add_argument("-a","--growthrate_mean",default=1,type=float)
parser.add_argument("-s","--growthrate_std",default=.2,type=float)
parser.add_argument("-t","--timescale_decorrelation",default=3,type=float)
parser.add_argument("-o","--outfile",default=None)
args = parser.parse_args()

pop = list(np.random.normal(args.growthrate_mean,args.growthrate_std,args.popsize_initial))
pop[pop < 0] = 0

a_heritable = np.exp(-1./args.timescale_decorrelation)
a_random    = 1. - a_heritable

for generation in range(args.popsize_initial,args.popsize_final):
    p = pop/np.sum(pop)
    new_gr = pop[np.random.choice(generation,1,p=p)[0]]
    pop.append(a_heritable * new_gr + a_random * np.random.normal(args.growthrate_mean,args.growthrate_std))

print '{:.6f} {:.6f} {:.6f}'.format(np.mean(pop),np.std(pop),np.median(pop))

if not args.outfile is None:
    h,b = np.histogram(pop,range = (args.growthrate_mean - 3*args.growthrate_std, args.growthrate_mean + 3*args.growthrate_std),bins = 100, density = True)
    b = b[:-1] + np.diff(b)
    np.savetxt(args.outfile,np.transpose([b,h]))



