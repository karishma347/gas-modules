# a python module to calculate gasdynamic oblique shock waves functions and their inverses
# Karishma Bharti, Srisha Rao M V, IISc, July 2020
# Global conventions and variables
# g=Ratio of Specific heats
# M=Mach no
# mue=prandtl meyer function
# theta= deflection angle

# modules that are imported
import math
def Prandtl_Meyer(g,M):
    g_r = (g+1)/(g-1)
    Mr =M**2-1
    mue_PM  = g_r**0.5*(math.atan((Mr/g_r)**0.5))*(180/math.pi)-(math.atan(Mr**0.5))*(180/math.pi)
    return mue_PM    #mue_PM is giving Prandtl Meyer Function

def inv_Prandtl_Meyer(M,x,g):
    g_r = (g+1)/(g-1)
    Mr =M**2-1
    M_PM  = (g_r**0.5*(math.atan((Mr/g_r)**0.5))*(180/math.pi)-(math.atan(Mr**0.5))*(180/math.pi))-x
    return M_PM    #M_PM is giving mach no for givem prandtl meyer function

def theta_PM(g,M1,M2):
    """theta=mue(M2)-mue(M1)"""
    mue_M1 = Prandtl_Meyer(g,M1)
    mue_M2 = Prandtl_Meyer(g,M2)
    theta = mue_M2-mue_M1
    return theta       #theta is deflection angle

