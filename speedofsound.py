# a python module to calculate gasdynamic speed of sound functions and their inverses
# Karishma Bharti, Srisha Rao M V, IISc, July 2020
# Global conventions and variables
# g=Ratio of Specific heats
# a=Speed of sound
# T=temperature
# rho=density
# M=Mach number
# V=Speed of object

import math
#functions
def soundspeedT(T,R,g):
    """This functions gives the value of speed of sound when temp is given"""
    a_T=math.sqrt(g*R*T)                     # Calculates the speed of sound given g,R,T
    return a_T

def soundspeedPrho(g,P,rho):
    """This functions gives the value of speed of sound when press, density are given"""
    a_prho=math.sqrt(g*P/rho)               # Calculates the speed of sound given g,P,rho
    return a_prho

def Tgivena(a,R,g):
    """This functions gives the temperature when speed of sound is given"""
    T=a*a/(g*R)                             # calculates temperature given a,g,R
    return T

def machnumber(V,a):
    """This functions gives the mach no when speed of object(V) and speed of sound(a) are given"""
    M=V/a                                   # calculates Mach number given V and a
    return M

def machangle(M):
    """This functions gives the mach angle(mue) when mach no is given using trignometric sin relation"""
    mue=math.asin(1.0/M)                    # calculates the Mach angle(mue)
    return mue

def Mgivenmu(mu):
    """This functions gives the mach no when mach angle(mue) is given using trignometric sin relation"""
    M_mue=1/math.sin(mu)                   # calculates the Mach number given the Mach angle(mue)
    return M_mue