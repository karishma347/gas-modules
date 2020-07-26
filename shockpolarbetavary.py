#Constant mach no, beta varying from 0 to 90
import numpy as np
import math
import matplotlib.pyplot as plt
def solver_OSW(g,M1):
    beta = list(range(0,90))
    beta = np.array(beta)
    b = beta*math.pi/180
    Mn_1 = M1*np.sin(b)
    pr_OSW = (1+((2*g)/(g+1))*(Mn_1**2-1))
    return pr_OSW    #pr_OSW is pressure ratio for OSW relation
p_r_OSW = solver_OSW(1.4,4)
print(p_r_OSW)
print('next')

def shockpolar(g,M1):
    b = list(range(0,90))
    b = np.array(b)
    beta = b*math.pi/180
    f1=((M1**2)*(np.sin(beta)**2))-1
    f2=(M1**2*(g+np.cos(2*beta)))+2
    theta_SP = (np.arctan(2*(1/np.tan(beta))*np.divide(f1,f2)))*(180/math.pi)
    return theta_SP
theta_OSW = shockpolar(1.4,4)
theta_OSW_list=[]
p_r_OSW_list=[]
for i in range(len(theta_OSW)):
    if theta_OSW[i]>0:
        print(theta_OSW)
        theta_OSW_list.append(theta_OSW[i])
        p_r_OSW_list.append(p_r_OSW[i]) 
theta_OSW_list=np.array(theta_OSW_list)
p_r_OSW_list=np.array(p_r_OSW_list)
plt.plot(theta_OSW_list,p_r_OSW_list)



