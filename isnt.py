# a python module to calculate gasdynamic isentropic functions and their inverses
# Karishma Bharti, Srisha Rao M V, IISc, July 2020
# Global conventions and variables
# g=Ratio of Specific heats
# TR=temperature ratio
# PR=pressure ratio
# ro=density ratio
#AR=area ratio
#M=mach no

# modules that are imported
import bisection

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

"""This function is giving mach no for given area ratio using another module bisection"""
"""method 1"""
#def isnt_invAR(x,g):
#    def isnt__AR(M,x,g):
#        M_Ar = ((1/M)*(( (2/(g+1))*(1+((g-1)/2)*M**2))**((g+1)/(2*(g-1)))))-x
#        return M_Ar   #M_Ar is mach no for given area ratio 
#    return bisection.bisection_NSW(0.01,1,1.5,1.4,isnt__AR)
#w=isnt_invAR(1.6,1.4)
#print(w)

"""method 2"""
def isnt_invAR(M,x,g):
    M_Ar = ((1/M)*(( (2/(g+1))*(1+((g-1)/2)*M**2))**((g+1)/(2*(g-1)))))-x
    return M_Ar   #M_Ar is mach no for given area ratio

""""method 3"""
#def isnt__AR(M,x,g):
#    M_Ar = ((1/M)*(( (2/(g+1))*(1+((g-1)/2)*M**2))**((g+1)/(2*(g-1)))))-x
#    return M_Ar   #M_Ar is mach no for given area ratio   
#
#def isnt_invAR(x,g):
#    return bisection.bisection_NSW(0.01,1,1.5,1.4,isnt__AR)
##w=isnt_invAR(1.6,1.4)
##print(w)