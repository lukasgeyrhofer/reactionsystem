#!/usr/bin/env python

import numpy as np
import argparse
import sys,math


def KullbackLeibler(dist1,dist2, dx = None):
    tmp1 = dist1[dist2 > 0]
    tmp2 = dist2[dist2 > 0]

    tmp2 = tmp2[tmp1 > 0]
    tmp1 = tmp1[tmp1 > 0]
    
    if np.sum(tmp1) > 0 and np.sum(tmp2) > 0:
        step = dx
        if dx is None:
           step = 1 
        tmp1 = tmp1/(np.sum(tmp1) * step)
        tmp2 = tmp2/(np.sum(tmp2) * step)

        return np.sum(tmp1*np.log(tmp1/tmp2))
    else:
        return None


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs = "*")
args = parser.parse_args()


allKL = list()

for filename in args.infiles:
    try:
        data = np.genfromtxt(filename)
    except:
        continue
        
    y       = data[:,0]
    dy      = y[1] - y[0]
    ONdistr = data[:,2]
    
    for i in range(3,len(data[0])):
        allKL.append(KullbackLeibler(data[:,i],ONdistr,dx = dy) * dy)
        print allKL[-1]

allKL = np.array(allKL)
allKL = allKL[allKL != np.array(None)]

print >> sys.stderr,np.mean(allKL),np.var(allKL)

h,b = np.histogram(allKL,range = (0,.4),bins = 100)
b = b[:-1] + 0.5* np.diff(b)

for x,y in zip(b,h):
    print x,y

