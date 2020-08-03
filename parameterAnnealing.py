import tellurium as te
import roadrunner
import antimony
import time
import copy
import random
import numpy as np
import math
#import scipy
#import pylab
#import seaborn as sns
import re

rateLaws = []                                                
rateLaws.append ("enzyme1*k1*(S0-S1/q0)/(1 + $S$/Ki$)")     
rateLaws.append ("enzyme2*k2*(S0-S2/q1)/(1 + $S$/Ki$)")    
rateLaws.append ("enzyme3*k3*(S2-S3/q2)/(1 + $S$/Ki$)")    
rateLaws.append ("enzyme4*k4*(S1-S0/q3)/(1 + $S$/Ki$)")     

modelWithRegulation = """
species $S0, S1, S2, $S3;


// Reactions:
J1: $S0 -> S1;   enzyme1*k1*(S0-S1/q0); 
J2: S0 -> S2;    enzyme2*k2*(S0-S2/q1)/(1 + S1/0.01); // Inhibited by S1
J3: S2 -> S3;    enzyme3*k3*(S2-S3/q2);
J4: S1 -> $S3;    enzyme4*k4*(S1-S3/q3)/(1 + S2/0.01);    // Inhibited by S2


factor = 1

enzyme1 = 1;   enzyme2 = 1;  enzyme3 = 1;  enzyme4 = 1; 

q0 = 1.5;  q1 = 1.5; q2 = 1.5;   q3 = 1.5;  q4 = 1.5;  

// Species initializations:
S0 = 10;
S1 = 0;
S2 = 0;
S3 = 5;


// Compartment initializations:
compartment_ = 1;

// Variable initializations:
k1 = 7
k2 = 0.8
k3 = 0.6
k4 = 3

"""

modelWithOutRegulation = """
species $S0, S1, S2, $S3;


// Reactions:
J1: $S0 -> S1;   enzyme1*k1*(S0-S1/q0); 
J2: S0 -> S2;    enzyme2*k2*(S0-S2/q1);      // Inhibited by S1
J3: S2 -> S3;    enzyme3*k3*(S2-S3/q2);
J4: S1 -> $S3;    enzyme4*k4*(S1-S3/q3);    // Inhibited by S2


factor = 1

enzyme1 = 1;   enzyme2 = 1;  enzyme3 = 1;  enzyme4 = 1; 

q0 = 1.5;  q1 = 1.5; q2 = 1.5;   q3 = 1.5;  q4 = 1.5;  

// Species initializations:
S0 = 10;
S1 = 0;
S2 = 0;
S3 = 5;


// Compartment initializations:
compartment_ = 1;

// Variable initializations:
k1 = 7
k2 = 0.8
k3 = 0.6
k4 = 3

"""

### Ground truth model
check = antimony.loadString(modelWithRegulation)
if check < 0:
    print(antimony.getLastError())
sbml_Reg = antimony.getSBMLString("__main")
r =  roadrunner.RoadRunner(sbml_Reg)


truth_K = [7, 0.8, 0.6, 3]
    
truth_S = [14.9777, 3.34025]
       

#truth_J = [ 6.66282856e+00,  3.50581162e-01,  2.56490726e-01,  9.32610851e-03,
       #3.32055336e+00,  1.14570471e-01,  7.67271858e-02,  5.69228603e-01,
       # 3.48804944e-02,  3.43545366e-01, -1.00035830e+00,  9.32610851e-03,
       # 4.47873886e-02,  3.19397972e-02,  8.85716966e+00,  3.48804944e-02,
       # 1.68458764e-01,  1.68458764e-01, -9.37229492e-01, -3.29089532e+00,
       # 4.92501417e-01,  2.56490726e-01, -3.41929268e-01,  6.66282856e+00,
       # 3.32055336e+00,  3.19397972e-02,  1.23013249e+00,  9.09384794e-01,
        #9.32610851e-03,  2.09691252e+00,  2.09691252e+00,  2.56490726e-01,
      # -2.51861735e+00,  9.32610851e-03,  7.67271858e-02, -3.27928745e+00,
        #4.58861287e+00,  9.42639826e+00]
    
fitUsingSensitivityMatrix = True
fitUsingConcentrations = True
fitUsingFlux_J1 = False
fitUsingFlux_All = False

r.steadyState()
truth_CS = r.getScaledConcentrationControlCoefficientMatrix()

antimony.loadString(modelWithOutRegulation)
sbml_noReg = antimony.getSBMLString("__main")
r =  roadrunner.RoadRunner(sbml_noReg)

parameters = ['k1', 'k2', 'k3', 'k4']

nReactions = r.getNumReactions()
reactionIds = r.getReactionIds()

def getFitness (rp):
    try:
        diffSqr_sum = 0
        rp.steadyState()
        if fitUsingSensitivityMatrix:
           diff = truth_CS - rp.getScaledConcentrationControlCoefficientMatrix()
           diffSqr = np.multiply (diff, diff)
           diffSqr_sum = np.sum (diffSqr)
        
        if fitUsingConcentrations:
           diff_S = truth_S - rp.getFloatingSpeciesConcentrations()
           diffSqr_S = np.multiply (diff_S, diff_S)
           diffSqr_sum = diffSqr_sum + np.sum (diffSqr_S)
    
    #    if fitUsingFlux_J1:
    #       diff_J = r.J1 - truth_J1
    #       diffSqr_J = diff_J*diff_J
    #       diffSqr_sum = diffSqr_sum + diffSqr_J
           
        #if fitUsingFlux_All:
           #for index, flux in enumerate (reactionIds):
               #diff_J = rp.getValue (flux) - truth_J[index]
               #diffSqr_J = diff_J*diff_J
               #diffSqr_sum = diffSqr_sum + diffSqr_J          
           
        return math.sqrt (diffSqr_sum)
    except:
        return "This didn't work 10000"
    


nSpecies = r.getNumFloatingSpecies()
speciesNames = r.getFloatingSpeciesIds()

currentFitness = 10000

# Mutate a model
def mutate (rp):
    # Pick a reaction
    p = random.randint (0,nReactions-1)
    ratelaw = rateLaws[p-1]
    ri = reactionIds[p-1]
    
    px = -1
    while px == -1:
       # pick a species that will regulate
       px = random.randint (0, nSpecies-1)
       si = speciesNames[px]
       # We are not allowed to pick a regulator that is also 
       # a reactant/product of the chosen reaction.
       if ratelaw.find (si) != -1:
          p = random.randint (0,nReactions-1)
          ratelaw = rateLaws[p-1]
          ri = reactionIds[p-1]          
          px = -1
       else:
          px = 0
      
    #print (ri, si, ratelaw)
    p = random.random()
    # Swap in a regulated step
    if ratelaw.find('$S$') != -1: 
    # Empty rate law
        if p>0.5:
            ratelaw = ratelaw.replace ('$S$', si)
            ratelaw = ratelaw.replace ('Ki$', '0.01')
        else:  
            ratelaw = ratelaw.replace ('/(1 + $S$/Ki$)', '')
    elif ratelaw.find(('S\d\d/\d'),ratelaw) != None:
    # Populated rate law, double digit reg species
        pattern = re.search(('S\d\d/\d'),ratelaw)
        if p>0.5:
            newLaw = si + '/0'
            ratelaw = ratelaw.replace(pattern.group(),newLaw)
        else:
            ratelaw.replace(pattern.group(), '')
    elif ratelaw.find(('S\d/\d'),ratelaw) != None:
    # Populated rate law, single digit reg species
        pattern = re.search(('S\d/\d'),ratelaw)
        if p>0.5:
            newLaw = si + '/0'
            ratelaw = ratelaw.replace(pattern.group(),newLaw)
        else:
            ratelaw.replace(pattern.group(), '')
    else:
    # Rate law completely absent
        if p>0.5:
            newLaw = '/(1 + ' + si +'/0.01)'
   
    #print (ratelaw)
    rp.setKineticLaw (ri, ratelaw, True)      
    

def byFitness (elem):
    return elem[1]

def copyModel (rx):
   rx.saveState ('c:\\Users\\user\\LogFiles\\parameterAnnealing.txt')
   original = roadrunner.RoadRunner()
   original.loadState ('c:\\Users\\user\\LogFiles\\parameterAnnealing.txt')
   return original
        
        
seed = 12379 #12341
random.seed (seed)
np.random.seed (seed)
nGenerations = 100
popSize = 15   
useTemperature = False
numberOfTrials = 1 
derivedModels = []

for trials in range (numberOfTrials):
    pop = []
    for i in range (popSize):
        rp = te.loada (modelWithOutRegulation)
        mutate (rp)
        pop.append ([rp, getFitness (rp)])
       
    pop.sort(key=byFitness)  
    
    fitArray = [] 
    nElite = 2
        
    timeStart = time.time()
    Temperature = 1.0; actualGenerationsRan = nGenerations
    for gen in range (nGenerations):
           
        
        fitArray.append (pop[0][1])
        print('Gen: ', gen, ' Fitness = ', pop[0][1])
        if pop[0][1] < 0.0001:
           actualGenerationsRan = nGenerations
           break
            
        newPop = []
        # Copy over the elite
        for i in range (nElite):
            newPop.append ([copyModel (pop[i][0]), pop[i][1]])
                
        count = 0
        for i in range (popSize-nElite):
            # Pick a random individual
            r1 = random.randint(0,popSize-1)
    
            original = copyModel (pop[r1][0])
            originalFitness = getFitness (original)
            
            mutatedModel = copyModel (original)
            mutate (mutatedModel)         
            mutatedFitness = getFitness (mutatedModel)
            
            if useTemperature:
               a = math.pow (math.e, (originalFitness-mutatedFitness)/Temperature)
               if a > 1:
                  newPop.append ([mutatedModel, mutatedFitness])
               else:
                  newPop.append ([original, originalFitness])
            else:
               if originalFitness > mutatedFitness:
                  # New one is better, add it to the new population
                  newPop.append ([copyModel (mutatedModel), mutatedFitness])
               else:
                 # Its worse, 0.8 means we will likley keep the worse model
                 if random.random() < 0.5:
                     # Keep the worse one
                     newPop.append ([copyModel (mutatedModel), mutatedFitness])
                 else:
                     # Keep the original
                     newPop.append ([copyModel (original), originalFitness])
        if gen % 100 == 0:
            if Temperature > 0.01:
               Temperature = Temperature*0.95;
        pop = newPop    
        
        pop.sort(key=byFitness)

        regulatedSteps = []
        modelStr = te.sbmlToAntimony (pop[0][0].getCurrentSBML())    
        s1 = modelStr.find('ns:')
        s2 = modelStr.find('// Species')
        # Cut out the reaction section
        modelStr = modelStr[s1+4:s2-3]
        modelStr = modelStr.splitlines()
        for count, s in enumerate (modelStr):
            if s.find ('0.01);') != -1:
               sf = s
               pattern = " \+ S(.*?)/0.01"
               spe = 'S' + re.search(pattern, sf).group(1)
               regulatedSteps.append ([count+1, spe])
    print ('Trial = ', trials)
    regulatedSteps.append (getFitness (pop[0][0]))# Tag on the fitness value
    derivedModels.append (regulatedSteps)
    print ('Time taken = ', time.time() - timeStart)
    print (derivedModels)
