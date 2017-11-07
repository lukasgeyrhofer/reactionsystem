#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
from scipy import stats

class reactionsystem:
    def __init__(self,indexset = ""):
        
        # define populations and set initial conditions =0 for all of them
        self.__indexset = indexset.replace("0","")
        self.__numpops  = len(self.__indexset)
        self.__n        = dict()
        self.set_population(self.__indexset,0,permissive = True)

        # reactions are stored in these arrays
        # a single first reaction is already stored: "0" -> "0" with rate 0.
        self.__reactionrates = np.array((0.),dtype = float)
        self.__reactants     = np.array(("0"),dtype = str)
        self.__products      = np.array(("0"),dtype = str)
        self.__coefficients  = np.array(("0"),dtype = str)
        self.__numreactions  = 1

        # internal time tracking
        self.__time = 0.
        self.__steps = 0
    
    
    def load_populations_from_file(self,filename = None,permissive = False):
        try:
            data = np.genfromtxt(filename,dtype=(str,int))
        except:
            raise IOError("Could not load populations from file '%s'"%filename)
        for i in range(len(data[:,0])):
            p = data[i,0]
            n = data[i,1]
            self.set_population(p,n,permissive)
    
    
    def load_reactions_from_file(self,filename = None,permissive = False):
        try:
            fp = open(filename,"r")
        except:
            raise IOError("Could not load reactions from file '%s'"%filename)
        for reaction in fp.readline():
            e = reaction.split()
            if len(e) < 2:
                # skip lines with not at least two entries
                continue
            elif len(e) == 2:
                self.add_reaction(e[0],e[1],permissive = permissive)
            elif len(e) == 3:
                self.add_reaction(e[0],e[1],rate = float(e[2]),permissive = permissive)
            elif len(e) > 3:
                self.add_reaction(e[0],e[1],rate = float(e[2]),coefficients = e[3],permissive = permissive)
    
    
    def existing_populations(self,populations = "0"):
        # splits the string for the population in two parts, with those in the indexset and those which are not
        populations = str(populations).replace("0","")
        tmp_exist    = ""
        tmp_notexist = ""
            
        for p in populations:
            if (p in self.__indexset) and (not p in tmp_exist):
                tmp_exist += p
            if (not p in self.__indexset) and (not p in tmp_notexist):
                tmp_notexist += p
        
        return list([tmp_exist,tmp_notexist])

    
    def set_population(self,population = "0",value = 0,permissive = False):
        populations = self.existing_populations(population)
        
        for p in populations[0]:
            self.__n[p] = value
        if permissive:
            for p in populations[1]:
                self.__indexset += p
                self.__numpops  += 1
                self.__n[p]      = value

    
    def add_reaction(self,reactants,products,rate = 1.,coefficients = "0",permissive = False):
        correctreaction = True
        
        if rate <= 0:
            correctreaction = False
            
        if correctreaction:
            rpop = self.existing_populations(reactants)
            if not permissive and len(rpop[1]) > 0:
                correctreaction = False
            else:
                self.set_population(rpop[1],0,permissive = True)
            
            ppop = self.existing_populations(products)
            if not permissive and len(ppop[1]) > 0:
                correctreaction = False
            else:
                self.set_population(ppop[1],0,permissive = True)
            
            cpop = self.existing_populations(coefficients)
            if not permissive and len(cpop[1]) > 0:
                correctreaction = False
            else:
                self.set_population(cpop[1],0,permissive = True)
                
        if correctreaction:
            self.__reactants     = np.append(self.__reactants,reactants)
            self.__products      = np.append(self.__products,products)
            self.__reactionrates = np.append(self.__reactionrates,rate)
            self.__coefficients  = np.append(self.__coefficients,coefficients)
            self.__numreactions += 1
    
    
    def isavailable(self,populations = "0"):
        # check for any populations in the parameter string, if even one of them is 0 => return false
        a    = True
        pops = self.existing_populations(populations)
        
        if len(pops[0]) == 0 or len(pops[1]) > 0:
            a = False
        else:
            for p in pops[0]:
                if self.__n[p] == 0:
                    a = False
        return a
    

    def nextreaction(self):
        # need to build rates from reactionrate and coefficients
        # so far, only linear dependence implemented
        currentrates = np.zeros(self.__numreactions)
        for i in range(self.__numreactions):
            currentrates[i] = self.__reactionrates[i]
            for r in self.__coefficients[i].replace("0",""):
                currentrates[i] *= self.__n[r]
                    
        # pick next reaction
        totalrate = np.sum(currentrates)
        if totalrate > 0:
            currentrates /= totalrate
            nr = np.random.choice(np.arange(self.__numreactions),p = currentrates)
        else:
            nr = 0
        
        # return as tuple
        return nr,totalrate
        
    
    def step(self):
        # draw random numbers for next reaction until one is found with all reactants present
        nextreaction = 0
        totalrate = 0
        while not self.isavailable(self.__reactants[nextreaction]):
            nextreaction,totalrate = self.nextreaction()
        
        # only works, if something happens
        if totalrate > 0:
            for r in self.__reactants[nextreaction].replace("0",""):
                self.__n[r] -= 1
            for r in self.__products[nextreaction].replace("0",""):
                self.__n[r] += 1
            
            self.__steps += 1
            self.__time  += np.random.exponential(1./totalrate)
            return self.__steps
        else:
            return None


    def get_populations(self, populations = None):
        a = np.array((),dtype = int)
        if populations == None:
            listpops = self.__indexset
        else:
            listpops = populations.replace("0","")
        for r in listpops:
            a = np.append(a,self.__n[r])
        return a
    
    def get_time(self):
        return self.__time
    
    def set_time(self,time):
        self.__time = time
    
    def get_step(self):
        return self.__steps
    
    def print_reactions(self):
        if self.__numreactions > 1:
            print "# Reactants\tProducts\tRate\tCoefficients"
            print "# ============================================="
            for i in range(1,self.__numreactions):
                print "# %s\t->\t%s\t%e\t%s"%(self.__reactants[i],self.__products[i],self.__reactionrates[i],self.__coefficients[i])
        else:
            print "# No reactions defined"
    
    
    def is_present(self,reactant):
        a = True
        if isinstance(reactant,str):
            for r in reactant:
                if r in self.__indexset.replace("0",""):
                    if self.__n[r] == 0:
                        a = False
                else:
                    a = False
        else:
            a = False
        return a



def main():
    
    r = reactionsystem(indexset = "PNRAG")
    r.set_population("PNRA",100)
    r.add_reaction("PR","PP",1.)
    r.add_reaction("NR","NN",1.)
    r.add_reaction("P","PG",1.)
    r.add_reaction("GA","0",1.)
    r.add_reaction("PRA","R",1.)
    r.add_reaction("NRA","R",1.)
    r.print_reactions()
    
    while r.is_present("R"):
        s =  r.step()
        if s%10 == 0:
            print r.get_time(),r.get_populations()
    
if __name__ == "__main__":
    main()
    
    
    
    
    
        
