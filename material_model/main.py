#Author: Waleed Ahmad Mirza
#Reference: doi:/10.1016/S0045-7825(96)01232-7
#Date: 29/06/2022

import numpy as np
import matplotlib.pyplot as plt
########################################
#Notes:
#The units of stress is MPA and Temperature is Celcius
#######################################

## Material parameters
#Elasticity
E=1000
#Temperature dependent super-elastic parameters.
C_AS= 1
C_SA= 1
T_ASf =10
T_ASs =70
T_SAs=90
T_SAf=130
#The parameters below set rates
beta_AS=10
beta_SA=10
#Material parameter that goes with the inelastic strain
epsi_L=10
#temperature of the material
T=160
##Numerical parameters
#tolerance
tol=1E-6
#time step
delta_t = 1E-8
#quantity on x-axis
applied_strain = np.linspace(0, 0.3, 100)
#quantity on y-axis
stress =  np.zeros(len(applied_strain))

##Calculate the material paramters for the evolution equation
# Parameters for transformation between austenite and martensite phase
R_ASs=C_AS*T_ASs
R_ASf=C_AS*T_ASf
R_SAs=C_SA*T_SAs
R_SAf=C_SA*T_SAf
##Initial guess for strain=0
#initial stress
sigma = 0
#internal variable at the first time step.
xi_s_n0 = 0
#strain at the first time step
epsilon_s_n0=0
#Phase transformation variables at n time step
F_SAn = 0
F_ASn = 0
 #Other variables
stress_idx=0
debug= True
## Strain loop
for strain_n in applied_strain:

    #print(applied_strain[stress_idx],len(applied_strain))
    #Phase transformation to martensite
    F_AS=sigma - C_AS*T
    F_ASs=F_AS-R_ASs
    F_ASf=F_AS-R_ASf
    #Phase transformation to austini
    F_SA=sigma-C_SA*T
    F_SAs=F_SA-R_SAs
    F_SAf=F_SA-R_SAf

     #Calculate change in the transformation parameter
    dot_F_SA      = F_SA - F_SAn
    dot_F_AS      = F_AS - F_ASn


    #Define the transformation parameter H^{AS}
    if  (F_ASs >0  and  dot_F_AS>0): # and F_ASf<0
        H_AS=1
    else:
        H_AS=0
    #Define the transformation parameter H^{SA}
    if  (F_SAs<0 and F_SAf>0 and  dot_F_SA<0):
        H_SA=1
    else:
        H_SA=0


    #An aux parameter used in the time-intergration
    X=1-H_AS*beta_AS*dot_F_AS/(F_ASf)**2 - H_SA*beta_SA*dot_F_SA/(F_SAf)**2

    # Internal variable at time step n after Backward Euler
    epsilon_s    = (H_AS*beta_AS*dot_F_AS/(F_ASf)**2 + epsilon_s_n0)/X
    #Calculate stress
    sigma =E*(strain_n - epsi_L*epsilon_s)

    #Stress from the free energ funtional
    stress[stress_idx] =sigma

    if (debug):
        print(F_ASs,",  ", F_ASf, ",  ", dot_F_AS, ",   ", H_AS,",     ",H_SA,",   ", sigma,",   ",F_AS)

            #, H_AS, H_SA, sigma,  epsilon_s)

    #Updating the internal variables
    epsilon_s_n0= epsilon_s
    F_SAn = F_SA
    F_ASn = F_AS
    stress_idx=stress_idx+1

plt.plot(applied_strain, stress)
plt.show()
#print(stress)






