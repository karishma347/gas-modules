# a python module to calculate gasdynamic normal shock waves functions and their inverses
# Karishma Bharti, Srisha Rao M V, IISc, July 2020
# Global conventions and variables
# g=Ratio of Specific heats
# TR=temperature ratio
# PR=pressure ratio
# dS=change in density
# ro=density ratio
#PO=stagnation pressure ratio

# modules that are imported
import numpy as np
from Gasdynamics_module import bisection #when calling from other folder
# import bisection #when calling in the same folder
#Functions for Normal Shock Waves relations
def NSW_TR(g, M):
    """This function is giving temp ratio"""
    TR_NSW = (1 + ((2*g/(g+1))*(M**2-1)))*((2+((g-1)*M**2))/((g+1)*M**2))
    return TR_NSW    #TR_NSW is temp ratio

def NSW_PR(g, M):
    """This function is giving press ratio"""
    PR_NSW = (1+((2*g)/(g+1))*(M**2-1))
    return PR_NSW    #PR_NSW is pressure ratio

def NSW_ro(g, M):
    """This function is giving density ratio"""
    ro_NSW = ((g+1)*M**2) /(2+(g-1)*M**2)
    return ro_NSW    #ro_NSW is density ratio

def NSW_M2(g,M):
    """This function is giving downstream mach no"""
    M2_NSW = ((1+(((g-1)/2)*(M**2))) / ((g*M**2)-((g-1)/2)))**0.5
    return M2_NSW    #M_2_NSW is M2

def NSW_PO(g,M):
    """This function is giving stag press ratio using change in entropy relation"""
    t = (1 + ((2*g/(g+1))*(M**2-1)))*((2+((g-1)*M**2))/((g+1)*M**2))
    p = (1+((2*g)/(g+1))*(M**2-1))
    dS  = (1005*np.log(t)) - (287*np.log(p))
    DS = (-1*dS)/287
    POr_NSW = np.exp(DS)
    return POr_NSW     #POr_NSW is stagnation pressure ratio

def NSW_DS(g,M):
    """This function is giving change in entropy"""
    t = (1 + ((2*g/(g+1))*(M**2-1)))*((2+((g-1)*M**2))/((g+1)*M**2))   #t is temp ratio used in formula
    p = (1+((2*g)/(g+1))*(M**2-1))    #p is press ratio used in formula
    dS_NSW  = (1005*np.log(t)) - (287*np.log(p))
    return dS_NSW     #dS_NSW is change in entropy
    
def NSW_invPR(Pr_NSW,g):
    """This function is giving inverse of press ratio"""
    M_PR_NSW = (((Pr_NSW-1)*((g+1)/(2*g)))+1)**0.5
    return M_PR_NSW     #M_PR_NSW is Mach no for given press ratio

def NSW_inv_ro(rhor_NSW,g):
    """This function is giving inverse of density ratio"""
    Mro_NSW = ((2*rhor_NSW)/(g+1-rhor_NSW*(g-1)))**0.5
    return Mro_NSW     #Mro_NSW is Mach no for given density ratio

def NSW_invM2(x,g=1.4):    #x is M2
    """This function is giving inverse of downstream mach no(M2) as upstream mach no(M) using bisection module"""
    M2_NSW = lambda M: (((1+(((g-1)/2)*(M**2)))/((g*(M**2))-((g-1)/2)))**0.5)-x
    return bisection.bis(0.5,20,M2_NSW) #M2_NSW is M2 found by calling bisection module

def NSW_invTR(x,g=1.4):   #x is temp ratio
    """This function is giving inverse of temp ratio as upstream mach no(M) using bisection module"""
    M_TR_NSW = lambda M: ((1 + (((2*g)/(g+1))*(M**2-1)))*((2+((g-1)*M**2))/((g+1)*M**2))-x)
    return bisection.bis(0.5,20,M_TR_NSW)     #M_TR_NSW is Mach no for given temp ratio found by calling bisection module

def NSW_inv_entropy(x,g=1.4):   #x is change in entropy
    """This function is giving inverse of change in entropy as upstream mach no(M) using bisection module"""
    M_DS_NSW = lambda M: ((1005*np.log((1 + (((2*g)/(g+1))*(M**2-1)))*((2+((g-1)*M**2))/((g+1)*M**2)))) - (287*np.log(1+((2*g)/(g+1))*(M**2-1))))-x
    return bisection.bis(0.5,20,M_DS_NSW)       #M_DS_NSW is Mach no for given entropy found by calling bisection module

def NSW_invPO(x,g=1.4):     #x is stag press ratio
    """This function is giving inverse of stag press ratio as upstream mach no(M) using bisection module"""
    M_POr_NSW = lambda M:(np.exp(((-1*((1005*np.log((1 + ((2*g/(g+1))*(M**2-1)))*((2+((g-1)*M**2))/((g+1)*M**2)))) - (287*np.log((1+((2*g)/(g+1))*(M**2-1))))))/287)))-x        
    return bisection.bis(0.5,20,M_POr_NSW)      #Mach no for given press ratio found by calling bisection

