
import tellurium as te
import roadrunner
import time
import copy
import random
import numpy as np
import math
import scipy
import pylab
import seaborn as sns
import re
​
​
rateLaws = []                                                
rateLaws.append ("enzyme1*k1*(S28-S9/q0)/(1 + $S$/Ki$)")     
rateLaws.append ("enzyme2*k2*(S24-S34/q1)/(1 + $S$/Ki$)")    
rateLaws.append ("enzyme3*k3*(S24-S30/q2)/(1 + $S$/Ki$)")    
rateLaws.append ("enzyme4*k4*(S14-S0/q3)/(1 + $S$/Ki$)")     
rateLaws.append ("enzyme5*k5*(S20-S23/q4)/(1 + $S$/Ki$)")    
rateLaws.append ("enzyme6*k6*(S32-S24/q5)/(1 + $S$/Ki$)")    
rateLaws.append ("enzyme7*k7*(S29-S6/q6)/(1 + $S$/Ki$)")     
rateLaws.append ("enzyme8*k8*(S19-S29/q7)/(1 + $S$/Ki$)")    
rateLaws.append ("enzyme9*k9*(S15-S27/q8)/(1 + $S$/Ki$)")    
rateLaws.append ("enzyme10*k10*(S34-S31/q9)/(1 + $S$/Ki$)") 
rateLaws.append ("enzyme11*k11*(S27-S17/q10)/(1 + $S$/Ki$)")
rateLaws.append ("enzyme12*k12*(S33-S20/q11)/(1 + $S$/Ki$)")
rateLaws.append ("enzyme13*k13*(S26-S31/q12)/(1 + $S$/Ki$)")
rateLaws.append ("enzyme14*k14*(S26-S8/q13)/(1 + $S$/Ki$)") 
rateLaws.append ("enzyme15*k15*(S19-S17/q14)/(1 + $S$/Ki$)")
rateLaws.append ("enzyme16*k16*(S34-S15/q15)/(1 + $S$/Ki$)")
rateLaws.append ("enzyme17*k17*(S4-S12/q16)/(1 + $S$/Ki$)") 
rateLaws.append ("enzyme18*k18*(S27-S4/q17)/(1 + $S$/Ki$)") 
rateLaws.append ("enzyme19*k19*(S34-S14/q18)/(1 + $S$/Ki$)")
rateLaws.append ("enzyme20*k20*(S31-S17/q19)/(1 + $S$/Ki$)")
rateLaws.append ("enzyme21*k21*(S29-S24/q20)/(1 + $S$/Ki$)")
rateLaws.append ("enzyme22*k22*(S3-S11/q21)/(1 + $S$/Ki$)") 
rateLaws.append ("enzyme23*k23*(S32-S14/q22)/(1 + $S$/Ki$)")
rateLaws.append ("enzyme24*k24*(S17-S28/q23)/(1 + $S$/Ki$)")
rateLaws.append ("enzyme25*k25*(S23-S2/q24)/(1 + $S$/Ki$)") 
rateLaws.append ("enzyme26*k26*(S8-S20/q25)/(1 + $S$/Ki$)") 
rateLaws.append ("enzyme27*k27*(S14-S27/q26)/(1 + $S$/Ki$)")
rateLaws.append ("enzyme28*k28*(S34-S31/q27)/(1 + $S$/Ki$)")
rateLaws.append ("enzyme29*k29*(S22-S33/q28)/(1 + $S$/Ki$)")
rateLaws.append ("enzyme30*k30*(S27-S7/q29)/(1 + $S$/Ki$)") 
rateLaws.append ("enzyme31*k31*(S7-S17/q30)/(1 + $S$/Ki$)") 
rateLaws.append ("enzyme32*k32*(S30-S3/q31)/(1 + $S$/Ki$)") 
rateLaws.append ("enzyme33*k33*(S14-S25/q32)/(1 + $S$/Ki$)")
rateLaws.append ("enzyme34*k34*(S0-S22/q33)/(1 + $S$/Ki$)") 
rateLaws.append ("enzyme35*k35*(S6-S26/q34)/(1 + $S$/Ki$)") 
rateLaws.append ("enzyme36*k36*(S20-S16/q35)/(1 + $S$/Ki$)")
rateLaws.append ("enzyme37*k37*(S31-S1/q36)/(1 + $S$/Ki$)") 
rateLaws.append ("enzyme38*k38*(S21-S19/q37)/(1 + $S$/Ki$)")
​
modelWithRegulation = """
​
species S0, $S1, S4, S6;
species S7, S8, S14, S15;
species S19, S20, S24, S26;
species S27, S28, S30, S31;
species $S2, S3, $S9, $S11;
species $S12, $S16, S17, $S21;
species S22, S23, $S25, S29;
species $S32, S33, S34;
​
// Reactions:
J1: S28 -> $S9;   enzyme1*k1*(S28-S9/q0); 
J2: S24 -> S34;   enzyme2*k2*(S24-S34/q1);
J3: S24 -> S30;   enzyme3*k3*(S24-S30/q2);
J4: S14 -> S0;    enzyme4*k4*(S14-S0/q3)/(1 + S20/0.01);    // Inhibited by S20
J5: S20 -> S23;   enzyme5*k5*(S20-S23/q4);
J6: $S32 -> S24;  enzyme6*k6*(S32-S24/q5);
J7: S29 -> S6;    enzyme7*k7*(S29-S6/q6);
J8: S19 -> S29;   enzyme8*k8*(S19-S29/q7);
J9: S15 -> S27;   enzyme9*k9*(S15-S27/q8);
J10: S34 -> S31;  enzyme10*k10*(S34-S31/q9);
J11: S27 -> S17;  enzyme11*k11*(S27-S17/q10);
J12: S33 -> S20;  enzyme12*k12*(S33-S20/q11);
J13: S26 -> S31;  enzyme13*k13*(S26-S31/q12)/( + S29/0.01);  // Inhibited by S29
J14: S26 -> S8;   enzyme14*k14*(S26-S8/q13);
J15: S19 -> S17;  enzyme15*k15*(S19-S17/q14);
J16: S34 -> S15;  enzyme16*k16*(S34-S15/q15);
J17: S4 -> $S12;  enzyme17*k17*(S4-S12/q16);
J18: S27 -> S4;   enzyme18*k18*(S27-S4/q17);
J19: S34 -> S14;  enzyme19*k19*(S34-S14/q18);
J20: S31 -> S17;  enzyme20*k20*(S31-S17/q19);
J21: S29 -> S24;  enzyme21*k21*(S29-S24/q20);
J22: S3 -> $S11;  enzyme22*k22*(S3-S11/q21);
J23: $S32 -> S14; enzyme23*k23*(S32-S14/q22);
J24: S17 -> S28;  enzyme24*k24*(S17-S28/q23);
J25: S23 -> $S2;  enzyme25*k25*(S23-S2/q24);
J26: S8 -> S20;   enzyme26*k26*(S8-S20/q25);
J27: S14 -> S27;  enzyme27*k27*(S14-S27/q26);
J28: S34 -> S31;  enzyme28*k28*(S34-S31/q27);
J29: S22 -> S33;  enzyme29*k29*(S22-S33/q28);
J30: S27 -> S7;   enzyme30*k30*(S27-S7/q29);
J31: S7 -> S17;   enzyme31*k31*(S7-S17/q30);
J32: S30 -> S3;   enzyme32*k32*(S30-S3/q31);
J33: S14 -> $S25; enzyme33*k33*(S14-S25/q32);
J34: S0 -> S22;   enzyme34*k34*(S0-S22/q33);
J35: S6 -> S26;   enzyme35*k35*(S6-S26/q34);
J36: S20 -> $S16; enzyme36*k36*(S20-S16/q35);
J37: S31 -> $S1;  enzyme37*k37*(S31-S1/q36);
J38: $S21 -> S19; enzyme38*k38*(S21-S19/q37);
​
factor = 1
​
enzyme1 = 1;   enzyme2 = 1;  enzyme3 = 1;  enzyme4 = 1;  enzyme5 = 1;  enzyme6 = 1;
enzyme7 = 1;   enzyme8 = 1;  enzyme9 = 1;  enzyme10 = 1; enzyme11 = 1; enzyme12 = 1;
enzyme13 = 1;  enzyme14 = 1; enzyme15 = 1; enzyme16 = 1; enzyme17 = 1; enzyme18 = 1;
enzyme19 = 1;  enzyme20 = 1; enzyme21 = 1; enzyme22 = 1; enzyme23 = 1; enzyme24 = 1;
enzyme25 = 1;  enzyme26 = 1; enzyme27 = 1; enzyme28 = 1; enzyme29 = 1; enzyme30 = 1;
enzyme31 = 1;  enzyme32 = 1; enzyme33 = 1; enzyme34 = 1; enzyme35 = 1; enzyme36 = 1;
enzyme37 = 1;  enzyme38 = 1; 
​
q0 = 1.5;  q1 = 1.5; q2 = 1.5;   q3 = 1.5;  q4 = 1.5;  q5 = 1.5;  q6 = 1.5;  q7 = 1.5;
q8 = 1.5;  q9 = 1.5; q10 = 1.5;  q11 = 1.5; q12 = 1.5; q13 = 1.5; q14 = 1.5; q15 = 1.5;
q16 = 1.5; q17 = 1.5; q18 = 1.5; q19 = 1.5; q20 = 1.5; q21 = 1.5; q22 = 1.5; q22 = 1.5;
q23 = 1.5; q24 = 1.5; q25 = 1.5; q26 = 1.5; q27 = 1.5; q28 = 1.5; q29 = 1.5; q30 = 1.5;
q31 = 1.5; q32 = 1.5; q33 = 1.5; q34 = 1.5; q35 = 1.5; q36 = 1.5; q37 = 1.5; q38 = 1.5;
q39 = 4.5
​
// Species initializations:
S0 = 0;
S1 = 0;
S4 = 0;
S6 = 0;
S7 = 0;
S8 = 0;
S14 = 0;
S15 = 0;
S19 = 0;
S20 = 0;
S24 = 0;
S26 = 0;
S27 = 0;
S28 = 0;
S30 = 0;
S31 = 0;
S2 = 2;
S3 = 1;
S9 = 2;
S11 = 2;
S12 = 3;
S16 = 6;
S17 = 1;
S21 = 4;
S22 = 6;
S23 = 6;
S25 = 3;
S29 = 4;
S32 = 1;
S33 = 3;
S34 = 5;
​
// Compartment initializations:
compartment_ = 1;
​
// Variable initializations:
k1 = 7
k2 = 0.8
k3 = 0.6
k4 = 3
k5 = 2.4
k6 = 0.4
k7 = 0.9
k8 = 0.2
k9 = 1.1
k10 = 1.7
k11 = 7
k12 = 4.6
k13 = 5.7
k14 = 2.7
k15 = 7.6
k16 = 6.2
k17 = 1.7
k18 = 0.25
k19 = 9.8
k20 = 3
k21 = 5.5
k22 = 8.2
k23 = 7.8
k24 = 3.7
k25 = 2.8
k26 = 1.8
k27 = 6.7
k28 = 4.5
k29 = 6.5
k30 = 9.2
k31 = 3.8
k32 = 4.7
k33 = 5.8
k34 = 1.5
k35 = 0.5
k36 = 3.5
k37 = 4.1
k38 = 5.4
​
Ki20 = 0.01
Ki29 = 0.01
"""
​
modelWithOutRegulation = """
​
species S0, $S1, S4, S6;
species S7, S8, S14, S15;
species S19, S20, S24, S26;
species S27, S28, S30, S31;
species $S2, S3, $S9, $S11;
species $S12, $S16, S17, $S21;
species S22, S23, $S25, S29;
species $S32, S33, S34;
​
// Reactions:
J1: S28 -> $S9;   enzyme1*k1*(S28-S9/q0); 
J2: S24 -> S34;   enzyme2*k2*(S24-S34/q1);
J3: S24 -> S30;   enzyme3*k3*(S24-S30/q2);
J4: S14 -> S0;    enzyme4*k4*(S14-S0/q3);    // Inhibited by S20
J5: S20 -> S23;   enzyme5*k5*(S20-S23/q4);
J6: $S32 -> S24;  enzyme6*k6*(S32-S24/q5);
J7: S29 -> S6;    enzyme7*k7*(S29-S6/q6);
J8: S19 -> S29;   enzyme8*k8*(S19-S29/q7);
J9: S15 -> S27;   enzyme9*k9*(S15-S27/q8);
J10: S34 -> S31;  enzyme10*k10*(S34-S31/q9);
J11: S27 -> S17;  enzyme11*k11*(S27-S17/q10);
J12: S33 -> S20;  enzyme12*k12*(S33-S20/q11);
J13: S26 -> S31;  enzyme13*k13*(S26-S31/q12);  // Inhibited by S29
J14: S26 -> S8;   enzyme14*k14*(S26-S8/q13);
J15: S19 -> S17;  enzyme15*k15*(S19-S17/q14);
J16: S34 -> S15;  enzyme16*k16*(S34-S15/q15);
J17: S4 -> $S12;  enzyme17*k17*(S4-S12/q16);
J18: S27 -> S4;   enzyme18*k18*(S27-S4/q17);
J19: S34 -> S14;  enzyme19*k19*(S34-S14/q18);
J20: S31 -> S17;  enzyme20*k20*(S31-S17/q19);
J21: S29 -> S24;  enzyme21*k21*(S29-S24/q20);
J22: S3 -> $S11;  enzyme22*k22*(S3-S11/q21);
J23: $S32 -> S14; enzyme23*k23*(S32-S14/q22);
J24: S17 -> S28;  enzyme24*k24*(S17-S28/q23);
J25: S23 -> $S2;  enzyme25*k25*(S23-S2/q24);
J26: S8 -> S20;   enzyme26*k26*(S8-S20/q25);
J27: S14 -> S27;  enzyme27*k27*(S14-S27/q26);
J28: S34 -> S31;  enzyme28*k28*(S34-S31/q27);
J29: S22 -> S33;  enzyme29*k29*(S22-S33/q28);
J30: S27 -> S7;   enzyme30*k30*(S27-S7/q29);
J31: S7 -> S17;   enzyme31*k31*(S7-S17/q30);
J32: S30 -> S3;   enzyme32*k32*(S30-S3/q31);
J33: S14 -> $S25; enzyme33*k33*(S14-S25/q32);
J34: S0 -> S22;   enzyme34*k34*(S0-S22/q33);
J35: S6 -> S26;   enzyme35*k35*(S6-S26/q34);
J36: S20 -> $S16; enzyme36*k36*(S20-S16/q35);
J37: S31 -> $S1;  enzyme37*k37*(S31-S1/q36);
J38: $S21 -> S19; enzyme38*k38*(S21-S19/q37);
​
factor = 1
​
enzyme1 = 1;   enzyme2 = 1;  enzyme3 = 1;  enzyme4 = 1;  enzyme5 = 1;  enzyme6 = 1;
enzyme7 = 1;   enzyme8 = 1;  enzyme9 = 1;  enzyme10 = 1; enzyme11 = 1; enzyme12 = 1;
enzyme13 = 1;  enzyme14 = 1; enzyme15 = 1; enzyme16 = 1; enzyme17 = 1; enzyme18 = 1;
enzyme19 = 1;  enzyme20 = 1; enzyme21 = 1; enzyme22 = 1; enzyme23 = 1; enzyme24 = 1;
enzyme25 = 1;  enzyme26 = 1; enzyme27 = 1; enzyme28 = 1; enzyme29 = 1; enzyme30 = 1;
enzyme31 = 1;  enzyme32 = 1; enzyme33 = 1; enzyme34 = 1; enzyme35 = 1; enzyme36 = 1;
enzyme37 = 1;  enzyme38 = 1; 
​
q0 = 1.5;  q1 = 1.5; q2 = 1.5;   q3 = 1.5;  q4 = 1.5;  q5 = 1.5;  q6 = 1.5;  q7 = 1.5;
q8 = 1.5;  q9 = 1.5; q10 = 1.5;  q11 = 1.5; q12 = 1.5; q13 = 1.5; q14 = 1.5; q15 = 1.5;
q16 = 1.5; q17 = 1.5; q18 = 1.5; q19 = 1.5; q20 = 1.5; q21 = 1.5; q22 = 1.5; q22 = 1.5;
q23 = 1.5; q24 = 1.5; q25 = 1.5; q26 = 1.5; q27 = 1.5; q28 = 1.5; q29 = 1.5; q30 = 1.5;
q31 = 1.5; q32 = 1.5; q33 = 1.5; q34 = 1.5; q35 = 1.5; q36 = 1.5; q37 = 1.5; q38 = 1.5;
q39 = 4.5
​
// Species initializations:
S0 = 0;
S1 = 0;
S4 = 0;
S6 = 0;
S7 = 0;
S8 = 0;
S14 = 0;
S15 = 0;
S19 = 0;
S20 = 0;
S24 = 0;
S26 = 0;
S27 = 0;
S28 = 0;
S30 = 0;
S31 = 0;
S2 = 2;
S3 = 1;
S9 = 2;
S11 = 2;
S12 = 3;
S16 = 6;
S17 = 1;
S21 = 4;
S22 = 6;
S23 = 6;
S25 = 3;
S29 = 4;
S32 = 1;
S33 = 3;
S34 = 5;
​
// Compartment initializations:
compartment_ = 1;
​
// Variable initializations:
k1 = 7
k2 = 0.8
k3 = 0.6
k4 = 3
k5 = 2.4
k6 = 0.4
k7 = 0.9
k8 = 0.2
k9 = 1.1
k10 = 1.7
k11 = 7
k12 = 4.6
k13 = 5.7
k14 = 2.7
k15 = 7.6
k16 = 6.2
k17 = 1.7
k18 = 0.25
k19 = 9.8
k20 = 3
k21 = 5.5
k22 = 8.2
k23 = 7.8
k24 = 3.7
k25 = 2.8
k26 = 1.8
k27 = 6.7
k28 = 4.5
k29 = 6.5
k30 = 9.2
k31 = 3.8
k32 = 4.7
k33 = 5.8
k34 = 1.5
k35 = 0.5
k36 = 3.5
k37 = 4.1
k38 = 5.4
​
Ki20 = 0.01
Ki29 = 0.01
"""
​
​
# Ground truth model
r = te.loada(modelWithRegulation)
​
truth_K = [7, 0.8, 0.6, 3, 2.4, 0.4, 0.9, 0.2, 1.1, 1.7,
     7, 4.6, 5.7, 2.7, 7.6, 6.2, 1.7, 0.25, 9.8, 3, 
     5.5, 8.2, 7.8, 3.7, 2.8, 1.8, 6.7, 4.5, 6.5, 9.2,
     3.8, 4.7, 5.8, 1.5, 0.5, 3.5, 4.1, 5.4, ]
    
truth_S = [0.91564855, 2.09909339, 1.07680066, 2.76795806, 2.05978482,
       1.56575563, 1.41386331, 3.38155604, 3.06306073, 1.07036073,
       1.38501943, 2.07323065, 2.28516598, 0.96431429, 1.11917387,
       1.36461269, 3.32420847, 1.36414672, 2.51924525, 0.80311954,
       2.0440679 , 0.94820142]
​
truth_J = [ 6.66282856e+00,  3.50581162e-01,  2.56490726e-01,  9.32610851e-03,
        3.32055336e+00,  1.14570471e-01,  7.67271858e-02,  5.69228603e-01,
        3.48804944e-02,  3.43545366e-01, -1.00035830e+00,  9.32610851e-03,
        4.47873886e-02,  3.19397972e-02,  8.85716966e+00,  3.48804944e-02,
        1.68458764e-01,  1.68458764e-01, -9.37229492e-01, -3.29089532e+00,
        4.92501417e-01,  2.56490726e-01, -3.41929268e-01,  6.66282856e+00,
        3.32055336e+00,  3.19397972e-02,  1.23013249e+00,  9.09384794e-01,
        9.32610851e-03,  2.09691252e+00,  2.09691252e+00,  2.56490726e-01,
       -2.51861735e+00,  9.32610851e-03,  7.67271858e-02, -3.27928745e+00,
        4.58861287e+00,  9.42639826e+00]
    
fitUsingSensitivityMatrix = True
fitUsingConcentrations = True
fitUsingFlux_J1 = False
fitUsingFlux_All = False
​
r.steadyState()
truth_CS = r.getScaledConcentrationControlCoefficientMatrix()
​
​
r = te.loada (modelWithOutRegulation)
​
parameters = ['k1', 'k2', 'k3', 'k4', 'k5', 'k6', 'k7', 'k8', 'k9', 'k10', 'k11', 'k12', 'k13', 'k14', 'k15', 'k16',
 'k17', 'k18', 'k19', 'k20', 'k21', 'k22', 'k23', 'k24', 'k25', 'k26', 'k27', 'k28', 'k29', 'k30', 'k31', 'k32', 'k33',
 'k34', 'k35', 'k36', 'k37', 'k38']
​
nReactions = r.getNumReactions()
reactionIds = r.getReactionIds()
​
def func (x):
    for index, k in enumerate (parameters):
        r.setValue (k, x[index])
    
    diffSqr_sum = 0
    r.steadyState()
    if fitUsingSensitivityMatrix:
       diff = truth_CS - r.getScaledConcentrationControlCoefficientMatrix()
       diffSqr = np.multiply (diff, diff)
       diffSqr_sum = np.sum (diffSqr)
    
    if fitUsingConcentrations:
       diff_S = truth_S - r.getFloatingSpeciesConcentrations()
       diffSqr_S = np.multiply (diff_S, diff_S)
       diffSqr_sum = diffSqr_sum + np.sum (diffSqr_S)
​
#    if fitUsingFlux_J1:
#       diff_J = r.J1 - truth_J1
#       diffSqr_J = diff_J*diff_J
#       diffSqr_sum = diffSqr_sum + diffSqr_J
       
    if fitUsingFlux_All:
       for index, flux in enumerate (reactionIds):
           diff_J = r.getValue (flux) - truth_J[index]
           diffSqr_J = diff_J*diff_J
           diffSqr_sum = diffSqr_sum + diffSqr_J          
       
    return math.sqrt (diffSqr_sum)
​
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
           
        if fitUsingFlux_All:
           for index, flux in enumerate (reactionIds):
               diff_J = rp.getValue (flux) - truth_J[index]
               diffSqr_J = diff_J*diff_J
               diffSqr_sum = diffSqr_sum + diffSqr_J          
           
        return math.sqrt (diffSqr_sum)
    except:
        return 10000
    
​
​
nSpecies = r.getNumFloatingSpecies()
speciesNames = r.getFloatingSpeciesIds()
​
currentFitness = 10000
​
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
    if p < 0.5:
        # Swap in a regulated step
       ratelaw = ratelaw.replace ('$S$', si)
       ratelaw = ratelaw.replace ('Ki$', '0.01')        
    else:
       ratelaw = ratelaw.replace ('/(1 + $S$/Ki$)', '')
       
    #print (ratelaw)
    rp.setKineticLaw (ri, ratelaw, True)      
    
​
def byFitness (elem):
    return elem[1]
​
def copyModel (rx):
   rx.saveState ('c:\\tmp\\r1.txt')
   original = roadrunner.RoadRunner()
   original.loadState ('c:\\tmp\\r1.txt')
   return original
        
        
seed = 12379 #12341
random.seed (seed)
np.random.seed (seed)
nGenerations = 30
popSize = 15   
useTemperature = False
numberOfTrials = 1 
derivedModels = []
​
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
           
        pop.sort(key=byFitness)
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
​
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

  Testing git
