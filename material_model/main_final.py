import numpy as np
import matplotlib.pyplot as plt

# Material properties
E = 1000
eps_l = 0.1
CAS = 1
CSA = 1
TAS_f = 10.000575443
TAS_s = 70.000575443
TSA_s = 90.000575443
TSA_f = 130.000575443
betaSA = 10
betaAS = 10

toll = 1e-6
toll_ = 1e-3

# Boundary Conditions
T = 160
eps_1 = np.arange(0, 0.15+toll, 0.00001)
eps_2 = np.arange(0.15,-toll, -0.00001)
eps_3 = np.arange(0.,0.175-toll, 0.00001)
eps_4 = np.arange(0.175,-toll, -0.00001)
eps_5 = np.arange(0.,0.2-toll, 0.00001)
eps_6 = np.arange(0.2,-toll, -0.00001)
eps_7 = np.arange(0.,0.225-toll, 0.00001)
eps_8 = np.arange(0.225,-toll, -0.00001)
eps_9 = np.arange(0.,0.28-toll, 0.00001)
eps_10 = np.arange(0.28,-toll, -0.00001)
eps = np.concatenate((eps_1,eps_2,eps_3,eps_4,eps_5,eps_6,eps_7,eps_8,eps_9,eps_10))
print(eps)
# internal variable
xi_s = np.zeros(eps.size)

# Initialize stress storage
sig = np.zeros(eps.size)


FSA_n = CSA*T
FAS_n = CAS*T


for i in range(0, eps.size - 1):

    # Compute trial state
    eps_e_tr = eps[i + 1] - eps_l * xi_s[i]
    sig_tr = E * eps_e_tr

    # Check for Phase transformation at trail state
    FAS_tr = sig_tr - CAS * T
    FAS_s_tr = sig_tr - CAS * (T-TAS_s)
    FAS_f_tr = sig_tr - CAS * (T-TAS_f)

    if (FAS_s_tr > 0   and xi_s[i]<1 and  FAS_tr-FAS_n>0 )  :
        HAS = 1
    else:
        HAS = 0

    FSA_tr = sig_tr - CSA * T
    FSA_s_tr = sig_tr - CSA*(T-   TSA_s)
    FSA_f_tr = sig_tr - CSA*(T-   TSA_f)

    if (FSA_s_tr  < 0    and  xi_s[i]> 0  and FSA_tr-FSA_n<0):
        HSA = 1
    else:
        HSA = 0
#    print("The values of tranformation coefficients HAS  are ",HAS , " because FAS_s_tr*FSA_f_tr= ", FAS_s_tr*FAS_f_tr , " and HSA is ", HSA, "because FSA_s_tr*FSA_f_tr ", FSA_s_tr*FSA_f_tr," for load-step ", i+1)
    print("HAS is ", HAS, "because FAS_f_tr", FAS_f_tr,"FAS_s_tr ", FAS_s_tr, "  xi_s ", xi_s[i], ", ",FAS_s_tr > 0 )
    print("HSA is ", HSA)
  # In case no phase transformation (Elastic step)
    if (HAS == 0) and (HSA == 0):
        # Update internal variable and the next trial value
        xi_s[i + 1] = xi_s[i]
        sig[i + 1] = sig_tr
        FAS_n = FAS_tr
        FSA_n = FSA_tr
    else:

        xi_s_i = xi_s[i]

        eps_e = eps[i + 1] - eps_l * xi_s_i
        sig[i + 1] = E * eps_e
        FSA = sig[i + 1] - CSA * T

        FSA_n = sig[i] - CSA * T
        FSA_f = FSA + CSA * TSA_f

        FAS = sig[i + 1] - CAS * T
        FAS_n = sig[i] + CAS * T
        FAS_f = FAS + CAS * TAS_f

        #Residual










        #R = (xi_s_i - xi_s[i])*(FAS_f**2)  - HSA*(betaSA * xi_s_i) * (FSA - FSA_n) \
        #    - HAS * (betaAS * (1 - xi_s_i))*(FAS - FAS_n)

        #Hessian
        #G = (FAS_f)**2 -   2*(xi_s_i - xi_s[i])*FAS_f*E*eps_l - HSA*betaSA  * (FSA - FSA_n) \
        #    + HSA * E * eps_l * betaSA * xi_s_i  \
        #    + HAS * betaAS * (FAS - FAS_n) + HAS * betaAS * eps_l * E * (1 - xi_s_i)
        # For calculating relative tolerance
        iter=0


        RHS_AS = (xi_s_i - xi_s[i]) - HAS*betaAS*(FAS-FAS_n)*(1.0 - xi_s_i)/(FAS_f**2)
        RHS_AS = RHS_AS -  HSA*betaSA*xi_s_i*(FSA-FSA_n)/(FSA_f**2)


        RHS_SA = 0 # (FSA - CSA*TSA_f)**2*dh_SA - HSA*betaSA*(FSA-FSA_n)*xi_s_i
        R=np.sqrt( RHS_AS**2  + RHS_SA**2)
        NR_RHS = np.zeros(2)
        NR_MAT = np.zeros((2,2))
        NR_SOL = np.zeros(2)
        while abs(R) >= toll:


            eps_e = eps[i + 1] - eps_l * xi_s_i
            sig[i + 1] = E * eps_e
            #print(sig[i+1], R)
            FSA = sig[i + 1] - CSA * T
            FSA_n = sig[i] - CSA * T
            FSA_f = FSA + CSA * TSA_f

            FAS = sig[i + 1] - CAS * T
            FAS_n = sig[i] - CAS * T
            FAS_f = FAS + CAS * TAS_f

#            R = (xi_s_i - xi_s[i])*(FAS_f**2)  - HSA*(betaSA * xi_s_i) * (FSA - FSA_n) \
#             - HAS * (betaAS * (1 - xi_s_i))*(FAS - FAS_n)

#            G = (FAS_f)**2 -   2*(xi_s_i - xi_s[i])*FAS_f*E*eps_l - HSA*betaSA  * (FSA - FSA_n) \
#                 + HSA * E * eps_l * betaSA * xi_s_i  \
#                 + HAS * betaAS * (FAS - FAS_n) + HAS * betaAS * eps_l * E * (1 - xi_s_i)

            RHS_AS = (xi_s_i - xi_s[i]) - HAS*betaAS*(FAS-FAS_n)*(1.0 - xi_s_i)/(FAS_f**2)
            RHS_AS = RHS_AS -  HSA*betaSA*xi_s_i*(FSA-FSA_n)/(FSA_f**2)

            RHS_SA = 0 # xi_s_i-xi_s[i]   ; #  FSA_f**2*dh_SA - HSA*betaSA*(FSA-FSA_n)*xi_s_i;

            #Assemble the RHS
            NR_RHS[0] = RHS_AS;
            NR_RHS[1] = RHS_SA;

            NR_MAT[0,0] = 1 + HAS*betaAS*E*eps_l*(1.0 - xi_s_i)/(FAS_f**2)  + HAS*betaAS*(FAS-FAS_n)/(FAS_f**2) +2*E*eps_l*HAS*betaAS*(FAS-FAS_n)*(1.0 - xi_s_i)/(FAS_f**3)
            NR_MAT[0,0] = NR_MAT[0,0] - HSA*betaSA*(FSA-FSA_n)/(FSA_f**2) + HSA*betaSA*xi_s_i*E*eps_l/(FSA_f**2) - 2*E*eps_l*HSA*betaSA*xi_s_i*(FSA-FSA_n)/(FSA_f**3)

            ##      HAS*betaAS*(FAS-FAS_n) +E*eps_l*(betaAS*HAS*(1-xi_s_i)-2*(xi_s_i-xi_s[i])*FAS_f) + FAS_f**2;
            NR_MAT[0,1] =  0 # HAS*betaAS*(FAS-FAS_n) +E*eps_l*(betaAS*HAS*(1-xi_s_i)-2*(xi_s_i-xi_s[i])*FAS_f) ;
            NR_MAT[1,0] =  0 #-HSA*betaSA*(FSA-FSA_n) +E*eps_l*(betaSA*HSA*xi_s_i-2*dh_SA*FSA_f) ;
            NR_MAT[1,1] =  0 #-HSA*betaSA*(FSA-FSA_n) +E*eps_l*(betaSA*HAS*xi_s_i-2*dh_SA*FSA_f) + FSA_f**2;

            #NR_SOL =   -np.linalg.inv(NR_MAT).dot(NR_RHS)

#            dh_AS = dh_AS-NR_SOL[0];
#            dh_SA = dh_SA-NR_SOL[1];


            xi_s_i = xi_s_i -  RHS_AS/NR_MAT[0,0]

            R=np.sqrt( RHS_AS**2)

            iter=iter+1
            print("     NR iteration   ", iter," has residual ", R )

        xi_s[i+1] = xi_s_i
    print("     The converged value of stress is ", sig[i+1], "and the value of the internal variable is  ", xi_s[i+1])
    print("              ")

fig = plt.plot(eps, sig,'k')
plt.xlabel('Strain',fontsize=18)
plt.ylabel('Stress [MPa]',fontsize=18)
plt.title('Uniaxis loading',fontsize=20)
plt.grid()
plt.savefig('cyclic_loading.svg')

load_steps= np.arange(0, len(eps)/10000,1e-4 )
fig = plt.plot(load_steps, eps,'k')
plt.ticklabel_format(style='sci')
plt.title('Strain load function',fontsize=20)
plt.xlabel('Load steps',fontsize=18)
plt.ylabel('Strain',fontsize=18)
plt.grid()
plt.savefig('strain_load.svg')






















#fig, (ax1, ax2) = plt.subplots(1, 2)
#fig.suptitle('Horizontally stacked subplots')
#ax1.plot(eps, sig)
#ax2.plot(eps, sig)
#plt.show()

















