#import Gasdynamics_module.isnt as isnt
#import bisection as b
from Gasdynamics_module import OSW, isnt,NSW, PrandtlMeyer as PM
t_T = isnt.isnt_TR(1.4, 4)
print('Isen: Temp ratio is', t_T,isnt.isnt_TR.__doc__)
print(isnt.isnt_TR.__doc__)
p_P = isnt.isnt_PR(1.4, 4)
print('Isen:Press ratio is', p_P)
ro_RO = isnt.isnt_ro(1.4, 4)
print('Isen:Density ratio is',ro_RO)
Ar = isnt.isnt_AR(1.4, 4)
print('Isen:Area ratio is',Ar)
M_t = isnt.isnt_invTR(1.4,.834)
print('Isen:Mach no for given temp ratio is',M_t)
M_p = isnt.isnt_invPR(1.4,.528)
print('Isen:Mach no for given press ratio is',M_p)
M_ro = isnt.isnt_inv_ro(1.4,.6339)
print('Isen:Mach no for given desity ratio is',M_ro)
M_sub,M_sup=isnt.isnt_invAR(1.6)
print('Isen:Mach no for given subsonic area ratio is',M_sub )
print('Isen:Mach no for given subsonic area ratio is',M_sup)
t_21_NSW = NSW.NSW_TR(1.4, 1.66)
print('NSW: Temp ratio is',t_21_NSW)
p_21_NSW = NSW.NSW_PR(1.4, 1.66)
print('NSW:Press ratio is',p_21_NSW)
ro_21_NSW = NSW.NSW_ro(1.4,1.66)
print('NSW:Density ratio is',ro_21_NSW)
M2_NSW = NSW.NSW_M2(1.4,1.66)
print('NSW:M2 is',M2_NSW)
PO_21_NSW = NSW.NSW_PO(1.4,1.66)
print('NSW:Stag press ratio is',PO_21_NSW)
DS_NSW = NSW.NSW_DS(1.4,12.8)
print('NSW:Change in entropy is',DS_NSW)
M_p_NSW = NSW.NSW_invPR(7.8,1.4)
print('NSW:Mach no for given press ratio is',M_p_NSW)
M_ro_NSW = NSW.NSW_inv_ro(3.5,1.4)
print('NSW:Mach no for given desity ratio is',M_ro_NSW)
M2_NSW=NSW.NSW_invM2(0.65)
print('NSW:Mach no for given M2 is',M2_NSW)
M_T_NSW = NSW.NSW_invTR(2.2179)
print('NSW:Mach no for given temp ratio is',M_T_NSW)
M_DS_NSW = NSW.NSW_inv_entropy(2000)
print('NSW:Mach no for given change in entropy is',M_DS_NSW)
M_st_Pr_NSW = NSW.NSW_invPO(0.0211)
print('NSW:Mach no for given stag press ratio is',M_st_Pr_NSW)
Mn1 = OSW.OSW_Mn1(5,1.4,19)
print('OSW:Mn1 for given M1 is',Mn1)
Mn2 = OSW.OSW_Mn2(5,1.4,19)
print('OSW:Mn2 for given M1 is',Mn2)
theta=OSW.OSW_theta(4,1.4,50)
print('OSW:theta from theta-beta-M relation is', theta)
M_2 = OSW.OSW_M2(10,1.4,29.8)
print('OSW:M2 is',M_2)
pr = OSW.OSW_PR(5,1.4,19.37)
print('OSW:Press ratio is',pr)
ro_21_OSW = OSW.OSW_ro(5,1.4,19.37)
print('OSW:Density ratio is',ro_21_OSW)
t_21_OSW = OSW.OSW_TR(5,1.4,19.37)
print('OSW: Temp ratio is',t_21_OSW)
M_1_OSW=OSW.OSW_M1(1.6278,1.4,19)
print('OSW:M1 for given Mn1 is',M_1_OSW)
M1n_p_OSW,M1_p_OSW = OSW.OSW_invPR(7.04,1.4,30)
print('OSW:M1n for given press ratio is',M1n_p_OSW)
print('OSW:M1 for given press ratio is',M1_p_OSW)
M1n_ro_OSW, M1_ro_OSW = OSW.OSW_inv_ro(3.31,1.4,30)
print('OSW:M1n for given density ratio is',M1n_ro_OSW)
print('OSW:M1 for given density ratio is',M1_ro_OSW)
M2n_NSW = OSW.OSW_M2n(1.65,1.4)
print('OSW:M2n for given M1n is ',M2n_NSW)
M1n_T_OSW, M1_T_OSW= OSW.OSW_invTR(2.94,30,g=1.4)
print('OSW:Mach no for given temp ratio is',M1n_T_OSW)
print('OSW:Mach no for given temp ratio is',M1_T_OSW)
mue = PM.Prandtl_Meyer(1.67,2.8)
print('PM: mue is', mue)
M_mue =PM.inv_Prandtl_Meyer(47.7898)
print('PM: Mach no for given mue is',M_mue)
theta_PMF = PM.theta_PM(1.4, 10,6.4)
print('PM:theta for PMF is', theta_PMF)