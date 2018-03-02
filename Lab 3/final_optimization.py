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
upsilon = .175
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



# Optimizing Objective Function
def optimize(deltamax, Efoam):
    
    AR_range = np.arange(8, 18, .2)
    S_range = np.arange(.01, 5, .02)
    Cl_range = np.arange(.8, 1.5, .02)
    W_pay_range = np.arange(1, 5, .2)
    
    payload_vals = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0] # V_pay, N, t_revpay, Cl, CD, e, t_reqpay, deltab, W_pay, productivity
    empty_vals = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0] # V_empty, N, t_revempty, Cl, calc_CD, e, T_reqempty, deltab, 0, productivity
    interest_vals = [0, 0, 0]
    max_productivity = 0
    
    for i in AR_range:
        print i
        for j in S_range:
            for k in Cl_range:
                for l in W_pay_range:
                    
                    AR, S, Cl, W_pay = i, j, k, l
                    
                    V_pay = calc_V(AR, S, Cl, W_pay, taper)
                    t_revpay = calc_t_rev(AR, S, Cl, W_pay, taper) # time of revolution with payload
                    T_maxpay = calc_Tmax(AR, S, Cl, W_pay, taper) # maximum thrust available 
                    CD = calc_CD(AR, S, Cl, W_pay, taper, upsilon)
                    T_reqpay = calc_T(AR, S, Cl, W_pay, taper, upsilon) # Thrust actually used
                    
                    if T_reqpay <= T_maxpay:
                        
                        V_empty = calc_V(AR, S, Cl, 0, taper) # Velocity                      
                        t_revempty = calc_t_rev(AR, S, Cl, 0, taper) # time of revolution for etpty plane
                        T_maxempty = calc_Tmax(AR, S, Cl, 0, taper) # maximum thrust of empty plane   
                        CD = calc_CD(AR, S, Cl, W_pay, taper, upsilon)
                        T_reqempty = calc_T(AR, S, Cl, 0, taper, upsilon) # required thust oempty plane                  
                        
                        if T_reqempty <= T_maxempty: 
                            
                            productivity = calc_productivity(W_pay, t_revpay, t_revempty)
                            
                            if productivity > max_productivity:    
                                
                                if (calc_deltab(AR, S, Cl, W_pay, taper, Efoam) <= deltamax) and (calc_deltab(AR, S, Cl, 0, taper, Efoam) <= deltamax):
                                    
                                    max_productivity = productivity
                                    
                                    payload_vals[0] = V_pay
                                    payload_vals[1] = calc_N(AR, S, Cl, W_pay, taper)
                                    payload_vals[2] = t_revpay
                                    payload_vals[3] = Cl
                                    payload_vals[4] = CD
                                    payload_vals[5] = calc_e(AR, S, Cl, W_pay, taper, upsilon)
                                    payload_vals[6] = T_reqpay
                                    payload_vals[7] = calc_deltab(AR, S, Cl, W_pay, Efoam, taper)
                                    payload_vals[8] = W_pay
                                    payload_vals[9] = productivity
                                    
                                    empty_vals[0] = V_empty
                                    empty_vals[1] = calc_N(AR, S, Cl, 0, taper)
                                    empty_vals[2] = t_revempty
                                    empty_vals[3] = Cl
                                    empty_vals[4] = CD
                                    empty_vals[5] = calc_e(AR, S, Cl, 0, taper, upsilon)
                                    empty_vals[6] = T_reqempty
                                    empty_vals[7] = calc_deltab(AR, S, Cl, 0, Efoam, taper)
                                    empty_vals[8] = 0
                                    empty_vals[9] = productivity
                                    
                                    interest_vals[0] = AR
                                    interest_vals[1] = S
                                    interest_vals[2] = Efoam
                                
                                else:
                                    pass
                            else:
                                pass
                        else:
                            pass
                    else:
                        pass                        
    
    return (payload_vals, empty_vals, interest_vals)


# Optimizing Objective Function, this time for wing parameters
def optimize2(deltamax, Efoam, AR, S, Cl, productivity):
    taper_range = np.arange(.3, 1.01, .01)
    upsilon_range = np.arange(.05, .5, .01)
    W_pay_range = np.arange(1, 6, .05)
    
    payload_vals = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0] # V_pay, N, t_revpay, Cl, CD, e, t_reqpay, deltab, W_pay, productivity
    empty_vals = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0] # V_empty, N, t_revempty, Cl, calc_CD, e, T_reqempty, deltab, 0, productivity
    interest_vals = [0, 0, 0, 0, 0]
    max_productivity = productivity
    
    for j in taper_range:
        print j
        for k in upsilon_range:
            for l in W_pay_range:
                
                taper, upsilon, W_pay = j, k, l
                V_pay = calc_V(AR, S, Cl, W_pay, taper)
                t_revpay = calc_t_rev(AR, S, Cl, W_pay, taper) # time of revolution with payload
                T_maxpay = calc_Tmax(AR, S, Cl, W_pay, taper) # maximum thrust available 
                CD = calc_CD(AR, S, Cl, W_pay, taper, upsilon)
                T_reqpay = calc_T(AR, S, Cl, W_pay, taper, upsilon) # Thrust actually used
                
                if T_reqpay <= T_maxpay:
                    
                    V_empty = calc_V(AR, S, Cl, 0, taper) # Velocity                      
                    t_revempty = calc_t_rev(AR, S, Cl, 0, taper) # time of revolution for etpty plane
                    T_maxempty = calc_Tmax(AR, S, Cl, 0, taper) # maximum thrust of empty plane   
                    CD = calc_CD(AR, S, Cl, W_pay, taper, upsilon)
                    T_reqempty = calc_T(AR, S, Cl, 0, taper, upsilon) # required thust oempty plane                  
                    
                    if T_reqempty <= T_maxempty: 
                        
                        productivity = calc_productivity(W_pay, t_revpay, t_revempty)
                        
                        if productivity > max_productivity:    
                            
                            if (calc_deltab(AR, S, Cl, W_pay, taper, Efoam) <= deltamax) and (calc_deltab(AR, S, Cl, 0, taper, Efoam) <= deltamax):
                                
                                max_productivity = productivity
                                
                                payload_vals[0] = V_pay
                                payload_vals[1] = calc_N(AR, S, Cl, W_pay, taper)
                                payload_vals[2] = t_revpay
                                payload_vals[3] = Cl
                                payload_vals[4] = CD
                                payload_vals[5] = calc_e(AR, S, Cl, W_pay, taper, upsilon)
                                payload_vals[6] = T_reqpay
                                payload_vals[7] = calc_deltab(AR, S, Cl, W_pay, taper, Efoam)
                                payload_vals[8] = W_pay
                                payload_vals[9] = productivity
                                
                                empty_vals[0] = V_empty
                                empty_vals[1] = calc_N(AR, S, Cl, 0, taper)
                                empty_vals[2] = t_revempty
                                empty_vals[3] = Cl
                                empty_vals[4] = CD
                                empty_vals[5] = calc_e(AR, S, Cl, 0, taper, upsilon)
                                empty_vals[6] = T_reqempty
                                empty_vals[7] = calc_deltab(AR, S, Cl, 0, taper, Efoam)
                                empty_vals[8] = 0
                                empty_vals[9] = productivity
                                
                                interest_vals[0] = AR
                                interest_vals[1] = S
                                interest_vals[2] = Efoam
                                interest_vals[3] = taper
                                interest_vals[4] = upsilon
                            
                            else:
                                pass
                        else:
                            pass
                    else:
                        pass
                else:
                    pass 
    
    return (payload_vals, empty_vals, interest_vals)



op_results = optimize(.1, Efoam)
payload_outputs = op_results[0]
empty_outputs = op_results[1]
interest_outputs = op_results[2]

print "With Payload Values:"
print "V = " + str(payload_outputs[0])
print "N = " + str(payload_outputs[1])
print "t_rev = " + str(payload_outputs[2])
print "Cl = " + str(payload_outputs[3])
print "Cd = " + str(payload_outputs[4])
print "e = " + str(payload_outputs[5])
print "T = " + str(payload_outputs[6])
print "delta/b = " + str(payload_outputs[7])
print "W_pay = " + str(payload_outputs[8])
print "Productivity = " + str(payload_outputs[9])
print 
print "No Payload Values:"
print "V = " + str(empty_outputs[0])
print "N = " + str(empty_outputs[1])
print "t_rev = " + str(empty_outputs[2])
print "Cl = " + str(empty_outputs[3])
print "Cd = " + str(empty_outputs[4])
print "e = " + str(empty_outputs[5])
print "T = " + str(empty_outputs[6])
print "delta/b = " + str(empty_outputs[7])
print "W_pay = " + str(empty_outputs[8])
print "Productivity = " + str(empty_outputs[9])
print 
print "Additional Values of Interest:"
print "AR = " + str(interest_outputs[0])
print "S = " + str(interest_outputs[1])
print "E foam = " + str(interest_outputs[2])

AR = interest_outputs[0]
S = interest_outputs[1]
Cl = payload_outputs[3]
productivity = payload_outputs[9]



print 
print 'RUN OPTIMIZATION AGAIN TO OPTIMIZE WING PARAMETERS'
print          

op_results = optimize2(.1, Efoam, AR, S, Cl, productivity)
payload_outputs = op_results[0]
empty_outputs = op_results[1]
interest_outputs = op_results[2]

print "With Payload Values:"
print "V = " + str(payload_outputs[0])
print "N = " + str(payload_outputs[1])
print "t_rev = " + str(payload_outputs[2])
print "Cl = " + str(payload_outputs[3])
print "Cd = " + str(payload_outputs[4])
print "e = " + str(payload_outputs[5])
print "T = " + str(payload_outputs[6])
print "delta/b = " + str(payload_outputs[7])
print "W_pay = " + str(payload_outputs[8])
print "Productivity = " + str(payload_outputs[9])
print 
print "No Payload Values:"
print "V = " + str(empty_outputs[0])
print "N = " + str(empty_outputs[1])
print "t_rev = " + str(empty_outputs[2])
print "Cl = " + str(empty_outputs[3])
print "Cd = " + str(empty_outputs[4])
print "e = " + str(empty_outputs[5])
print "T = " + str(empty_outputs[6])
print "delta/b = " + str(empty_outputs[7])
print "W_pay = " + str(empty_outputs[8])
print "Productivity = " + str(empty_outputs[9])
print 
print "Additional Values of Interest:"
print "AR = " + str(interest_outputs[0])
print "S = " + str(interest_outputs[1])
print "E foam = " + str(interest_outputs[2])
print "taper = " + str(interest_outputs[3])
print "dihedral (in radians) = " + str(interest_outputs[4])




