#Constant mach no, beta varying from 0 to 90 and -90 to 0
import numpy as np
import math
import matplotlib.pyplot as plt
def solver_OSW(g,M1):
    beta1 = list(range(0,90))
    beta1 = np.array(beta1)
    b1 = beta1*math.pi/180
    beta2 = list(range(-90,0))
    beta2 = np.array(beta2)
    b2 = beta2*math.pi/180
    Mn_1_1 = M1*np.sin(b1)
    Mn_1_2 = M1*np.sin(b2)
    pr_OSW_1 = (1+((2*g)/(g+1))*(Mn_1_1**2-1))
    pr_OSW_2 = (1+((2*g)/(g+1))*(Mn_1_2**2-1))
    return pr_OSW_1, pr_OSW_2   #pr_OSW is pressure ratio for OSW relation
p_r_OSW_1,pr_OSW_2 = solver_OSW(1.4,4)
print(p_r_OSW_1,pr_OSW_2)
print('next')

def shockpolar(g,M1):
    b1 = list(range(0,90))
    b1 = np.array(b1)
    beta1 = b1*math.pi/180
    f11=((M1**2)*(np.sin(beta1)**2))-1
    f21=(M1**2*(g+np.cos(2*beta1)))+2
    theta_SP1 = (np.arctan(2*(1/np.tan(beta1))*np.divide(f11,f21)))*(180/math.pi)
    b2 = list(range(-90,0))
    b2 = np.array(b2)
    beta2 = b2*math.pi/180
    f12=((M1**2)*(np.sin(beta2)**2))-1
    f22=(M1**2*(g+np.cos(2*beta2)))+2
    theta_SP2 = (np.arctan(2*(1/np.tan(beta2))*np.divide(f12,f22)))*(180/math.pi)
    return theta_SP1,theta_SP2
theta_OSW1, theta_OSW2 = shockpolar(1.4,4)
theta_OSW_list1=[]
p_r_OSW_list1=[]
for i in range(len(theta_OSW1)):
    if theta_OSW1[i]>=0:
        print(theta_OSW1)
        theta_OSW_list1.append(theta_OSW1[i])
        p_r_OSW_list1.append(p_r_OSW_1[i]) 
theta_OSW_list1=np.array(theta_OSW_list1)
p_r_OSW_list1=np.array(p_r_OSW_list1)
theta_OSW_list2=[]
p_r_OSW_list2=[]
for i in range(len(theta_OSW2)):
    if theta_OSW2[i]<=0:
        print(theta_OSW2,'o')
        theta_OSW_list2.append(theta_OSW2[i])
        p_r_OSW_list2.append(pr_OSW_2[i]) 
theta_OSW_list2=np.array(theta_OSW_list2)
p_r_OSW_list2=np.array(p_r_OSW_list2)
plt.plot(theta_OSW_list1,p_r_OSW_list1)
plt.plot(theta_OSW_list2,p_r_OSW_list2)



