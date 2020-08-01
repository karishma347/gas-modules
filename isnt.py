# a python module to calculate gasdynamic isentropic functions and their inverses
# Karishma Bharti, Srisha Rao M V, IISc, July 2020
# Global conventions and variables
# g=Ratio of Specific heats
# TR=temperature ratio
# PR=pressure ratio
# ro=density ratio
# AR=area ratio
# M=mach number

# modules that are imported
# import bisection                                   #when calling in the same folder
from Gasdynamics_module import bisection             #when calling from other folder
#Functions for isentropic relations
def isnt_TR(g, M):
    """This function is giving temp ratio"""
    T_R= 1/(1+((g-1)/2)*(M**2))
    return T_R   #T_R is temp ratio

def isnt_PR(g,M): 
    """This function is giving press ratio"""
    P_R = 1/((1+((g-1)/2)*(M**2))**(g/(g-1)))
    return P_R    #P_R is pressure ratio

def isnt_ro(g,M):
    """This function is giving density ratio"""
    ro_r = 1/((1+((g-1)/2)*(M**2))**(1/(g-1)))
    return ro_r   #ro_r is density ratio

def isnt_AR(g,M):
    """This function is giving area ratio"""
    A_R = (1/M)*(( (2/(g+1))*(1+((g-1)/2)*M**2))**((g+1)/(2*(g-1))))
    return A_R   #A_R is area ratio

def isnt_invTR(g, T_21):
    """This function is giving mach no for given temp ratio"""
    M_TR = ((2*((1/T_21)-1))/(g-1))**0.5
    return M_TR  #M_TR is Mach no for given temp ratio

def isnt_invPR(g,P_21):
    """This function is giving mach no for given press ratio"""
    M_PR = ((2*(((1/P_21)**((g-1)/g))-1))/(g-1))**0.5
    return M_PR  #M_PR is Mach no for given press ratio

def isnt_inv_ro(g,ro_21):
    """This function is giving mach no for given density ratio"""
    M_ro = ((2*(((1/ro_21)**(g-1))-1))/(g-1))**0.5
    return M_ro    #M_ro is Mach no for given density ratio

def isnt_invAR(x,g=1.4):
  """This function is giving mach no for given area ratio using another module bisection""" 
  M_Ar =lambda M: ((1/M)*(( (2/(g+1))*(1+((g-1)/2)*M**2))**((g+1)/(2*(g-1)))))-x     
  M_sub= bisection.bis(0.01,1,M_Ar)     #M_sub gives subsonic mach no for given area ratio
  M_sup= bisection.bis(1,20,M_Ar)       #M_sup gives supersonic mach no for given area ratio
  return M_sub,M_sup
