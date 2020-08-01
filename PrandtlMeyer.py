# a python module to calculate gasdynamic oblique shock waves functions and their inverses
# Karishma Bharti, Srisha Rao M V, IISc, July 2020
# Global conventions and variables
# g=Ratio of Specific heats
# M=Mach no
# mue=prandtl meyer function
# theta= deflection angle

# modules that are imported
from Gasdynamics_module import bisection           #when calling from other folder
# import bisection                                 #when calling in the same folder
import math
def Prandtl_Meyer(g,M):
    g_r = (g+1)/(g-1)
    Mr =M**2-1
    mue_PM  = g_r**0.5*(math.atan((Mr/g_r)**0.5))*(180/math.pi)-(math.atan(Mr**0.5))*(180/math.pi)
    return mue_PM                                  #mue_PM is giving Prandtl Meyer Function

def inv_Prandtl_Meyer(x,g=1.4):
    g_r = (g+1)/(g-1)
    M_PM  =lambda M : (g_r**0.5*(math.atan(((M**2-1)/g_r)**0.5))*(180/math.pi)-(math.atan((M**2-1)**0.5))*(180/math.pi))-x
    M_mu = bisection.bis(1,20,M_PM)
    return M_mu                                    #M_mu is giving mach no for given prandtl meyer function

def theta_PM(g,M1,M2):
    """theta=mue(M2)-mue(M1)"""
    mue_M1 = Prandtl_Meyer(g,M1)
    mue_M2 = Prandtl_Meyer(g,M2)
    theta = mue_M2-mue_M1
    return theta                                   #theta is deflection angle
