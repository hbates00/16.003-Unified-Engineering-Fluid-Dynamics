# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 18:12:30 2017

@author: Haley and Matthew
"""

#something has changed 

import numpy as np
import math

''' ----------- Global Variables ----------- '''

# Plane Vanilla Parameters
S_pv = 0.228
b_pv = 1.52
e_0 = 0.9136 #straight and level span efficiency
taper = 0.1/0.2
tau = 0.11

# foam properties for weight and deflection estimates
g = 9.81    #m/s^2
rho_air = 1.225     #kg/m^2
mu = 1.8e-5
pfoam = 33.0        #kg/m^2
upsilon = .1823
Efoam = 19e6

T_max = .7      #N
R = 12.5    #m

''' ----------- Functions ----------- '''
# Calculate wing weight -- check
def calc_W_wing(AR, S, taper):
    return ((4/5.0) * pfoam * S * tau * np.sqrt(S / AR) * ((taper**2 + taper + 1) / (taper + 1)**2) * g)
    
# Calculate fuselage weight (kilograms) -- check
def calc_W_fuse(AR, S, taper):
    return (0.145 + ((0.06 * (np.sqrt(S * AR)) / b_pv)) + (0.045 * (S / S_pv))) * g
    
# Calculate profile drag -- check
def calc_prof_drag(AR, S, Cl, W_pay, taper):
    cd_0 = 0.020 * (1 + tau**2)
    cd_1 = -0.005 / (1 + 6 * tau)
    cd_2 = 0.160 / (1 + 60 * tau)
    cd_8 = 1.0
    cl_0 = 1.25 - 3 * tau
    Re_ref = 100000
    a = -0.75
    return (cd_0 + cd_1*(Cl - cl_0) + cd_2*(Cl - cl_0)**2 + cd_8*(Cl - cl_0)**8) * (calc_Re(AR, S, Cl, W_pay, taper) / Re_ref)**a
    
# Calculate CDA_0 -- check
def calc_CDA_0(S):
    return 0.002 + 0.002 * (S / S_pv)

# Calculating Thrust -- check
def calc_T(AR, S, Cl, W_pay, taper, upsilon):   
    return (.5 * rho_air * calc_V(AR, S, Cl, W_pay, taper)**2 * S * calc_CD(AR, S, Cl, W_pay, taper, upsilon))

'''
# Calculate W_Pay
def calc_W_pay(AR, S, Cl):
    return ((T_max / (((calc_CDA_0(S) / S)/Cl) + (calc_prof_drag(AR, S, Cl) / Cl) + (Cl / (np.pi * AR * calc_e(AR, S, Cl))))) - calc_W_fuse(AR, S) - calc_W_wing(AR, S))
'''

# Calculating Drag Coefficient -- check
def calc_CD(AR, S, Cl, W_pay, taper, upsilon): 
    return (calc_CDA_0(S) / S) + calc_prof_drag(AR, S, Cl, W_pay, taper) + ((Cl**2) / (np.pi * AR * calc_e(AR, S, Cl, W_pay, taper, upsilon)))

# Calculate N -- check
def calc_N(AR, S, Cl, W_pay, taper):
    return (1 - ((calc_W_fuse(AR, S, taper) + calc_W_wing(AR, S, taper) + W_pay) / (.5 * rho_air * g * R * S * Cl))**2)**(-1/2.0)

# Calculating Velocity from CL -- check
def calc_V(AR, S, Cl, W_pay, taper):
    return (np.sqrt((2 * (W_pay + calc_W_wing(AR, S, taper) + calc_W_fuse(AR, S, taper))/(Cl * rho_air * S))))

#Calculating time of revolution -- check
def calc_t_rev(AR, S, Cl, W_pay, taper):
    return ((2 * math.pi * R) / calc_V(AR, S, Cl, W_pay, taper))
    
# Caclulate max thrust
def calc_Tmax(AR, S, Cl, W_pay, taper):
    return 1 - 0.08 * (calc_V(AR, S, Cl, W_pay, taper)) - 0.0028 * (calc_V(AR, S, Cl, W_pay, taper))**2

# Calculate delta/b -- check
def calc_deltab(AR, S, Cl, W_pay, taper, Efoam ):
    return (0.018) * calc_N(AR, S, Cl, W_pay, taper) * (1 + taper)**3 * (1 + 2 * taper) * (AR**3 / S) * ((calc_W_fuse(AR, S, taper) + W_pay) / (Efoam * tau * (tau**2 + calc_epsilon()**2)))

# Calculate Reynolds Number -- check
def calc_Re(AR, S, Cl, W_pay, taper):
    return (rho_air * calc_V(AR, S, Cl, W_pay, taper) * np.sqrt(S / AR)) / mu

# Calculate camber/chord ratio
def calc_epsilon():
    return 0.10 - 0.5 * tau
    
# Calculate r_bar, nondimensionalized yaw rate -- check
def calc_r_bar(AR, S, Cl, W_pay, taper):
    return np.sqrt(S * AR) / (2 * R * calc_N(AR, S, Cl, W_pay, taper))

# Calculate Beta, slideslip angle (radians) -- check
def calc_Beta(AR, S, Cl, W_pay, taper, upsilon):
    t1 = Cl / upsilon
    t2 = 1 + (4.0 / AR)
    t3 = calc_r_bar(AR, S, Cl, W_pay, taper) / (2 * np.pi)
    return t1 * t2 * t3

# Calculate span efficiency -- check
def calc_e(AR, S, Cl, W_pay, taper, upsilon):
    return e_0 * (1 - 0.5 * (calc_r_bar(AR, S, Cl, W_pay, taper))**2) * (np.cos(calc_Beta(AR, S, Cl, W_pay, taper, upsilon))**2)

# Calculate Productivity, the objective function -- check
def calc_productivity(W_pay, t_revpay, t_revempty):
    return W_pay / (t_revpay + t_revempty)

AR = 11.2
S = 0.33
W_pay = 0
Cl = 1.22

print calc_e(AR, S, Cl, W_pay, taper, upsilon)
