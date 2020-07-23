# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 09:29:05 2020

@author: user
"""

import re

ratelaw = "enzyme1*k1*(S28-S9/q0)"

#ratelaw.replace('S','$$$$')
si = 'S12'


newLaw = '/(1 + ' + si +'/0.01)'

ratelaw = ratelaw + newLaw