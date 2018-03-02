#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 20:59:04 2017

@author: root
"""
import math
import numpy

'''
--- PRELIMINARY CONFIGURATION ---

currently using AR and S optimized for max Wpay with d/bmax = 0.1

'''

print('~ Preliminary tail configuration ~')
print('assuming horizontal tail moment arm of 65 cm')


rho_foam1 = 25.5 #DOW Blue
t_foam = 0.05

c_PV = 0.15 # plane vanilla values
b_PV = 1.52
S_PV = 0.228
l_h_PV = 0.65 # 25.6 inches, taken from plane vanilla 3-view
nose_fraction = 1.1 # distance from tip to leading MAC edge, normalized with c
dh = 10

c = 0.225 # preliminary
b = 2.49
AR = 11.042 
S = 0.5637
CL = 0.93

SM = 0.1

V_h = 0.45
V_v = 0.035
x_cg_norm = 0.365
l_h = 0.975 # scaled up from plane vanilla
l_v = 0.975

m_wing = rho_foam * S * t_foam
m_fuse_0 = .145
m_fuse_l = 0.06
m_fuse_S = 0.045

x_nose = -c*nose_fraction
x_cg_wing = c/2

AR_h_range = numpy.arange(1, 20, 0.001)
for i in AR_h_range:
    AR_h = i
    x_np_norm = 0.25 + (1+2/AR)/(1+2/AR_h) * (1-4/(AR+2)) * V_h
    if abs(x_np_norm - x_cg_norm - SM) < 0.001:
        print("Static margin:", x_np_norm - x_cg_norm)
        print("Horizontal tail AR:", AR_h)
        break
    





S_h_range = numpy.arange(0.001, 0.2, 0.0001)
for i in S_h_range:
    S_h = i
    V_h = (S_h*l_h)/(S*c)
    x_np_norm = 0.25 + (1+2/AR)/(1+2/AR_h) * (1-4/(AR+2)) * V_h
                       
    if abs(x_np_norm - x_cg_norm - SM) < 0.001:
        print("Static margin:", x_np_norm - x_cg_norm)
        print("Horizontal tail S:", S_h)
        break
    
    



S_v = V_v*S*b/l_v
print("Vertical tail S:", round(S_v,4))

B = l_v * dh / (b * CL)
print("B:", round(B, 3))


brat = b/b_PV
Srat = S/S_PV

m_tot = m_fuse_0 + m_fuse_l * brat + m_fuse_S * Srat
