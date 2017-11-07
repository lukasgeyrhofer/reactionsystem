#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==================================================================== #
#                                                                      #
#  Simulation program to investigate growing populations of cells,     #
#  that have a decorrelation timescale tau for inheriting yield        #
#  from their parent                                                   #
#                                                                      #
#  First, an ON (overnight) culture is created, that starts from       #
#  'ON_seedingsize' cells for 'ON_generations' generations. This       #
#  ON culture is used to seed multiple droplets with 'seedingsize'     #
#  cells initially, this time for 'generation' generations.            #
#                                                                      #
#  Inheritance follows a stochastic process,                           #
#    Y(n+1) = a Y(n) + (1-a) xi,     xi ~ Uniform(YieldMin,YieldMax),  #
#  where 1/tau = log(1/a)                                              #
#                                                                      #
#  Lukas Geyrhofer, 2017                                               #
#  l.geyrhofer@technion.ac.il                                          #
#                                                                      #
# ==================================================================== #


import numpy as np
import argparse
import sys



class inoculumeffect(object):
    def __init__(self,**kwargs):
        # parse cmdline arguments
        self.__PoissonSeeding       = kwargs.get("PoissonSeeding",False)
        
        self.__generations          = kwargs.get("generations",8)
        self.__correlation          = kwargs.get("correlation",8)
        self.__seedingsize          = kwargs.get("seedingsize",25)
        if not self.__PoissonSeeding:
            self.__seedingsize      = int(self.__seedingsize)
        
        self.__ONinitialcorrelation = kwargs.get("ON_initialcorrelation",8)
        self.__ONgenerations        = kwargs.get("ON_generations",8)
        self.__ONseedingsize        = int(kwargs.get("ON_seedingsize",25))
        
        self.__yieldinterval        = np.array([kwargs.get("yieldmin",.5),kwargs.get("yieldmax",1.5)])
        self.__verbose              = kwargs.get("verbose",False)
        self.__onlymeanhisto        = kwargs.get("onlymeanhisto",False)
        self.__outputgenerationstep = kwargs.get("outputgenerationstep",1)
        if not self.__outputgenerationstep is None:
            self.__intermediateoutput = True
        else:
            self.__intermediateoutput = False   
        
        # coefficients for faster reference instead of computing them every step
        self.__coefficient          = np.array([np.exp(-1./self.__correlation),1. - np.exp(-1./self.__correlation)])
        
        # statistics, analysis
        self.__histogrambins       = 20
        self.__histograms          = list()
        self.__finalpopulationsize = list()
        
        # have startingconditions?
        self.__haveovernightculture = False
    
    
    
    def rng(self):
        return np.random.uniform(low = self.__yieldinterval[0], high = self.__yieldinterval[1])
            
    def newyield(self,xn):
        return self.__coefficient[0] * xn + self.__coefficient[1] * self.rng()
    
    
    # seed an ON culture to generate starting conditions for single droplets
    def run_overnightculture(self,seedingsize = None, generations = None, initialcorrelation = None):
        if  seedingsize        is None:
            seedingsize        = self.__ONseedingsize
        if  initialcorrelation is None:
            initialcorrelation = self.__ONinitialcorrelation
        if  generations        is None:
            generations        = self.__ONgenerations
        
        self.__overnightculture = list()
        
        # make the initial seeding for the overnight culture
        x = self.rng()
        for i in range(seedingsize):
            self.__overnightculture.append(x)
            # wait 'initialcorrelation' generations before adding a new value, this is only a rough estimate of this distribution
            for j in range(int(initialcorrelation)):
                x = self.newyield(x)
        
        # from these initial seedings, run on average g generations
        self.__currentsubstrate = np.power(2.,generations) * seedingsize / np.mean(self.__overnightculture)
        
        # add more cells
        while self.add(population = "overnightculture"):
            continue
        
        # starting substrate chosen such that the ON culture would take on average g generations to use up all nutrients
        self.__ONyieldmean_inv = 1./ np.mean(self.__overnightculture)
        self.__startingsubstrate = np.power(2.,self.__generations) * self.__seedingsize * self.__ONyieldmean_inv
        
        # we're done here
        self.__haveovernightculture = True
    
    # seed a droplet and let cells grow until substrate is depleted
    def run(self,seedingsize = None, generations = None):
        # need overnightculture for seeding
        if not self.__haveovernightculture:
            self.run_overnightculture()
        
            
        # use default values from object creation if no argument given here
        if  seedingsize is None:
            seedingsize = self.__seedingsize
        if  generations is None:
            generations = self.__generations
            self.__currentsubstrate = self.__startingsubstrate
        else:
            self.__currentsubstrate = np.power(2,generations) * seedingsize * self.__ONyieldmean_inv

        if self.__intermediateoutput:
            outgen  = np.arange(start = 0, stop = generations, step = self.__outputgenerationstep)
            outsize = np.power(2,outgen) * seedingsize * self.__ONyieldmean_inv
            current_outsize_index = 0

        
        if self.__PoissonSeeding:
            seedingsize = np.random.poisson(seedingsize)
        
        if seedingsize > 0:
            # set initial conditions
            self.__population = list(np.random.choice(self.__overnightculture,size = seedingsize))

            # run until nutrients are out
            while self.add():
                if self.__intermediateoutput:
                    if len(self.__population) >= outsize[current_outsize_index]:
                        self.verbose("# gen: {} size: {}".format(outgen[current_outsize_index],outsize[current_outsize_index]))
                        current_outsize_index += 1
        else:
            # empty droplet due to Poisson seeding
            self.__population = list()
        
        
        # do statistics on run
        self.__histograms.append(np.histogram(self.__population,range = self.__yieldinterval, bins = self.__histogrambins))
        fps = len(self.__population)
        self.__finalpopulationsize.append(fps)
        return fps

    # add a single cell to the population, return False if not enough substrate anymore
    def add(self,population = "population"):
        # use dict representation of self to chose either "self.__population" or "self.__overnightculture"
        x  = self.newyield(np.random.choice(self.__dict__["_inoculumeffect__{:s}".format(population)]))
        xi = 1./x
        if self.__currentsubstrate > xi:
            self.__currentsubstrate -= xi
            self.__dict__["_inoculumeffect__{:s}".format(population)].append(x)
            return True
        else:
            return False
    
    # output funneled through this method
    def verbose(self,msg = "", handle = None, flush = False):
        if self.__verbose:
            if handle is None:
                handle = sys.stdout
            print >>handle, msg
            if flush:
                handle.flush()
    

    def __getattr__(self,key):
        # reading out those attributes resets them to empty!
        if key == "finalpopulationsize":
            fps = self.__finalpopulationsize
            self.__finalpopulationsize = list()
            return fps
        elif key == "histograms":
            if len(self.__histograms) > 0:
                bins  = self.__histograms[0][1][:-1] + 0.5 * np.diff(self.__histograms[0][1])
                h     = np.concatenate([[x[0] for x in self.__histograms]],axis=0)
                meanh = np.mean(h,axis=0)
                
                # ON culture yield histo
                ONh,ONb = np.histogram(self.__overnightculture,range = self.__yieldinterval, bins = self.__histogrambins)
                

                r     = np.transpose([bins,meanh,ONh])
                if not self.__onlymeanhisto:
                    r = np.concatenate([r,h.T],axis=1)
                self.__histograms = list()
                return r
                
            else:
                raise ValueError("no histograms found. run the populations")
        elif "substraterange":
            return np.sort(np.power(2,self.__generations) * self.__seedingsize / self.__yieldinterval)







def main():
    # run multiple ON cultures with many 'droplets' (=second populations) each
    
    # add cmdline arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-g","--generations",           type = float, default = 8.)
    parser.add_argument("-t","--correlation",           type = float, default = 8.)
    parser.add_argument("-n","--seedingsize",           type = float, default = 25) # can assume float for PoissonSeeding
    parser.add_argument("-T","--ON_initialcorrelation", type = float, default = 8.)
    parser.add_argument("-G","--ON_generations",        type = float, default = 8.)
    parser.add_argument("-N","--ON_seedingsize",        type = int,   default = 25)
    
    parser.add_argument("-k","--droplets",              type = int, default = 1000)
    parser.add_argument("-O","--overnightculturecount", type = int, default = 3)

    parser.add_argument("-y","--yieldmin", type = float, default = 0.5)
    parser.add_argument("-Y","--yieldmax", type = float, default = 1.5)
    
    parser.add_argument("-P","--PoissonSeeding", default = False, action = "store_true")
    
    parser.add_argument("-H","--onlymeanhisto",                      default = False, action = "store_true")
    parser.add_argument("-v","--verbose",                            default = False, action = "store_true")
    parser.add_argument("-L","--logfile",                            default = None)
    parser.add_argument("-o","--outfilebasename",                    default = "out")
    parser.add_argument("-s","--outputgenerationstep", type = float, default = None)
    args = parser.parse_args()


    # create logfile, otherwise default later to sys.stdout
    try:
        logfile = open(args.logfile,"w")
        if args.verbose:
            print "opening logfile '{}'".format(args.logfile)
        if not args.verbose:
            args.verbose = True
    except:
        logfile = None


    # initialize object and datastructure
    ie = inoculumeffect(**vars(args))
    
    # loop over different ON cultures
    for i in range(args.overnightculturecount):
        ie.verbose("# starting overnight culture ({:4d}/{:4d})".format(i+1,args.overnightculturecount), handle = logfile, flush = True)
        ie.run_overnightculture()
    
        # seed droplets from this ON culture
        for j in range(args.droplets):
            current_fps = ie.run()
            ie.verbose("#   droplet ({:4d}/{:4d}) FPS {:d}".format(j+1,args.droplets,current_fps),handle = logfile)

        # reading destroys the data, so only read once
        fps         = ie.finalpopulationsize
        histo_yield = ie.histograms
        
        # compute moments of FPS
        fps_mean = np.mean(fps)
        fps_var  = np.var(fps)
        ie.verbose("# MomentsFPS:   {:.4e} {:.4e}".format(fps_mean,fps_var), handle = logfile)
        
        
        # make histogram for population sizes
        ps,psbin = np.histogram(fps,range = ie.substraterange,bins = 200)
        histo_fps = np.transpose(np.array([psbin[:-1] + 0.5 * np.diff(psbin),ps]))
        
        # save histograms to files
        np.savetxt("{}_N{:04d}".format(args.outfilebasename,i),histo_fps)
        np.savetxt("{}_Y{:04d}".format(args.outfilebasename,i),histo_yield)



if __name__ == "__main__":
    main()




