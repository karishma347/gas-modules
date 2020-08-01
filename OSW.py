# a python module to calculate gasdynamic oblique shock waves functions and their inverses
# Karishma Bharti, Srisha Rao M V, IISc, July 2020
# Global conventions and variables
# g=Ratio of Specific heats
# TR=temperature ratio
# PR=pressure ratio
# dS=change in density
# ro=density ratio
# PO=stagnation pressure ratio

# modules that are imported
import math
from Gasdynamics_module import bisection       #when calling from other folder
# import bisection                             #when calling in the same folder
#Functions for Oblique Shock waves Relations
def OSW_Mn1(M1, g, beta):
    """This function is giving normal component of upstream mach no (Mn1=M1*sin(beta))"""
    b = beta*math.pi/180
    Mn_1_OSW = M1*math.sin(b)
    return Mn_1_OSW    #Mn_1_OSW is Mn1

def OSW_Mn2(M1,g,beta):
    """This function is giving normal component of downstream mach no"""
    Mn_1 = OSW_Mn1(M1, g, beta)
    Mn_2_OSW = ((1+(((g-1)/2)*(Mn_1**2))) / ((g*Mn_1**2)-((g-1)/2)))**0.5
    return Mn_2_OSW    #Mn_2_OSW is Mn2

def OSW_theta(M1,g,beta):
    """This function is giving theta using theta-beta-M relation"""
    b = beta*math.pi/180
    f11=((M1**2)*(math.sin(b)**2))-1
    f21=(M1**2*(g+math.cos(2*b)))+2
    theta = (math.atan(2*(1/math.tan(b))*(f11/f21)))*(180/math.pi)
    return theta

def OSW_M2(M1,g,beta):
    """This function is giving downstream mach no"""
    Mn2 = OSW_Mn2(M1,g,beta)
    b = beta*math.pi/180
    th=OSW_theta(M1,g,beta)
    t = th*math.pi/180
    M2_OSW = Mn2/math.sin(b-t)
    return M2_OSW    #M2_OSW is M2 

def OSW_PR(M1, g,beta):
    """This function is giving press ratio"""
    Mn_1 = OSW_Mn1(M1,g,beta)
    PR_OSW = (1+((2*g)/(g+1))*(Mn_1**2-1))
    return PR_OSW    #PR_OSW is pressure ratio

def OSW_ro(M1, g,beta):
    """This function is giving density ratio"""
    Mn_1 = OSW_Mn1(M1,g,beta)
    ro_OSW = ((g+1)*Mn_1**2) /(2+(g-1)*Mn_1**2)
    return ro_OSW    #ro_OSW is density ratio

def OSW_TR(M1,g,beta):
    """This function is giving temp ratio"""
    Mn_1 = OSW_Mn1(M1,g,beta)
    TR_OSW = (1 + ((2*g/(g+1))*(Mn_1**2-1)))*((2+((g-1)*Mn_1**2))/((g+1)*Mn_1**2))
    return TR_OSW    #TR_OSW is temp ratio

def OSW_M1(M1n,g,beta):
    """This function is giving M1 from Mn1"""
    b = beta*math.pi/180
    M1_OSW= M1n/math.sin(b)   #M1_OSW is upstream mach no
    return M1_OSW

def OSW_invPR(Pr_OSW,g,beta):
    """This function is giving inverse of press ratio as normal component of upstream mach no and upstream mach no"""
    M1n_PR_OSW = (((Pr_OSW-1)*((g+1)/(2*g)))+1)**0.5   #M1n_PR_OSW is normal component of upstream Mach no for given press ratio
    b = beta*math.pi/180
    M1_PR_OSW= M1n_PR_OSW/math.sin(b)  #M1_PR_OSW upstream Mach no for given press ratio
    return M1n_PR_OSW, M1_PR_OSW    

def OSW_inv_ro(rhor_OSW,g,beta):
    """This function is giving inverse of density ratio as normal component of upstream mach no and upstream mach no"""
    M1n_ro_OSW = ((2*rhor_OSW)/(g+1-rhor_OSW*(g-1)))**0.5  #M1n_ro_OSW is normal component of upstream Mach no for given density ratio
    b = beta*math.pi/180
    M1_ro_OSW= M1n_ro_OSW/math.sin(b)    #M1_ro_OSW is upstream Mach no for given density ratio
    return M1n_ro_OSW, M1_ro_OSW     

def OSW_M2n(M1n,g):
    """This function is giving normal component downstream mach no for normal component of upstream mach no"""
    M2n_OSW = ((1+(((g-1)/2)*(M1n**2)))/((g*(M1n**2))-((g-1)/2)))**0.5
    return M2n_OSW     #M2n_OSW is M2n

def OSW_invTR(x,beta,g=1.4):
    """This function is giving inverse of temp ratio as normal component of upstream mach no and upstream mach no, using the bisection module"""
    M1n_TR_OSW =lambda Mn_1: ((1 + (((2*g)/(g+1))*(Mn_1**2-1)))*((2+((g-1)*Mn_1**2))/((g+1)*Mn_1**2))-x)
    Mn_1_TR_OSW=bisection.bis(0.5,20, M1n_TR_OSW)   #Mn_1_TR_OSW is  normal component of upstream Mach no for given temp ratio
    b = beta*math.pi/180
    M1_TR_OSW= Mn_1_TR_OSW/math.sin(b)   #M1_TR_OSW is upstream Mach no for given temp ratio
    return Mn_1_TR_OSW, M1_TR_OSW
 

