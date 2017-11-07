#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import reactionsystem as rs

def output(time,pops):
        print "{:8.3f}".format(time),
        pops = r.get_populations(allpops)
        for p in pops:
            print " {:5d}".format(p),
        print
    
parser = argparse.ArgumentParser()
parser.add_argument("-P","--populations",type=int,default=4)
parser.add_argument("-N","--initialcond_firstpop",type=int,default=25)
parser.add_argument("-n","--initialcond_otherpop",type=int,default=0)
parser.add_argument("-S","--substrate",type=int,default=10000)

parser.add_argument("-r","--repetitions",type=int,default=10)
parser.add_argument("-m","--mu",type=float,default=1e-2)
parser.add_argument("-a","--alpha",type=float,default=1.)
parser.add_argument("-o","--outputsteps",type=int,default=100)
args = parser.parse_args()

assert 2 <= args.populations <= 26,"populations indexed by letters in alphabet..."

r       = rs.reactionsystem(indexset = "Aa")
prevn   = "A"
allpops = "A"
r.add_reaction("Aa", "AA", rate = args.alpha, coefficients = "A")



for i in range(66,65+args.populations):
    n = chr(i)     # microbial population
    s = chr(i+32)  # consumed substrate
    r.add_reaction(n+s, n+n, rate = args.alpha, coefficients = n,     permissive = True)  # growth
    r.add_reaction(n, prevn, rate = args.mu,    coefficients = n,     permissive = True)  # migration to previous deme
    r.add_reaction(prevn, n, rate = args.mu,    coefficients = prevn, permissive = True)  # migration from previous deme
    prevn    = n
    allpops += n # keep whole string of all growing populations (for output)
    
for rep in range(args.repetitions):
    r.set_population("A",args.initialcond_firstpop)
    r.set_population("a",args.substrate)
    for i in range(66,65+args.populations):
        n = chr(i)
        s = chr(i+32)
        r.set_population(n,args.initialcond_otherpop)
        r.set_population(s,args.substrate)

    r.set_time(0)
    o=0
    while r.is_present("a"):
        if o%args.outputsteps == 0:
            output(r.get_time(),r.get_populations(allpops))
        o = r.step()
    output(r.get_time(),r.get_populations(allpops))
    print
    
