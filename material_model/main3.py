import numpy as np
import matplotlib.pyplot as plt
import sys
# Material properties
E = 1000
eps_l = 0.1
CAS = 1
CSA = 1
TAS_s = 70
TAS_f = 10
TSA_s = 90
TSA_f = 130
betaSA = 10
betaAS = 10

toll = 1e-8

# Boundary Conditions
T = 160
eps = np.arange(0, 0.30000001, 0.01)
#eps_2 = np.arange(0.30000001,0, -0.01)
#eps   = np.concatenate((eps_1, eps_2))

# internal variable
xi_s = np.zeros(eps.size)

# Initialize stress storage
sig = np.zeros(eps.size)

#Initially sigma=0
FAS = -CAS*T
FSA = -CSA*T

for i in range(0, eps.size - 1):

    # Compute trial state
    eps_e_tr = eps[i + 1] - eps_l * xi_s[i]
    sig_tr = E * eps_e_tr

    # Check for Phase transformation at trail state
    FAS_tr = sig_tr - CAS * T
    FAS_s_tr = FAS_tr - CAS * TAS_s
    FAS_f_tr = FAS_tr - CAS * TAS_f

#    print("*** ", FAS_s_tr)
#    if FAS_s_tr < 0 and FAS_tr<FAS:
#        HAS = 0
#    else:
#        HAS = 1

    if (FAS_s_tr*FAS_f_tr<0 ):
        HAS = 1
    else:
        HAS = 0

    FSA_tr = sig_tr - CSA * T
    FSA_s_tr = FSA_tr - CSA * TSA_s
    FSA_f_tr = FSA_tr - CSA * TSA_f


#    if FSA_s_tr > 0 and FSA_tr>FSA:
#        HSA = 0
#    else:
#        HSA = 1

    if (FSA_s_tr*FSA_f_tr<0):
        HSA = 1
    else:
        HSA = 0



    print("The values of tranformation coefficients HAS AND HSA are ",HAS, HSA,  " for load-step ", i+1)

    # In case no phase transformation (Elastic step)
    if (HAS == 0) and (HSA == 0):
        # Update internal variable and the next trial value
        xi_s[i + 1] = xi_s[i]
        sig[i + 1] = sig_tr
        FSA = FSA_tr
        FAS = FAS_tr

    else:
        xi_s_i = xi_s[i]

        eps_e = eps[i + 1] - eps_l * xi_s_i
        sig[i + 1] = E * eps_e
        FSA = sig[i + 1] - CSA * T

        FSA_n = sig[i] - CSA * T
        FSA_f = FSA - CSA * TSA_f

        FAS = sig[i + 1] - CAS * T
        FAS_n = sig[i] - CAS * T
        FAS_f = FAS - CAS * TAS_f

        #Residual
        R = (xi_s_i - xi_s[i])*(FAS_f**2)  - HSA*(betaSA * xi_s_i) * (FSA - FSA_n) \
            - HAS * (betaAS * (1 - xi_s_i))*(FAS - FAS_n)

        #Hessian
        G = (FAS_f)**2 -   2*(xi_s_i - xi_s[i])*FAS_f*E*eps_l - HSA*betaSA  * (FSA - FSA_n) \
            + HSA * E * eps_l * betaSA * xi_s_i  \
            + HAS * betaAS * (FAS - FAS_n) + HAS * betaAS * eps_l * E * (1 - xi_s_i)
        # For calculating relative tolerance
        iter=0
        while (abs(R) >= toll and iter<20):

            xi_s_i = xi_s_i - R/G
            eps_e = eps[i + 1] - eps_l * xi_s_i
            sig[i + 1] = E * eps_e
            #print(sig[i+1], R)
            FSA = sig[i + 1] - CSA * T
            FSA_n = sig[i] - CSA * T
            FSA_f = FSA - CSA * TSA_f

            FAS = sig[i + 1] - CAS * T
            FAS_n = sig[i] - CAS * T
            FAS_f = FAS - CAS * TAS_f

            R = (xi_s_i - xi_s[i])*(FAS_f**2)  - HSA*(betaSA * xi_s_i) * (FSA - FSA_n) \
             - HAS * (betaAS * (1 - xi_s_i))*(FAS - FAS_n)

            G = (FAS_f)**2 -   2*(xi_s_i - xi_s[i])*FAS_f*E*eps_l - HSA*betaSA  * (FSA - FSA_n) \
                 + HSA * E * eps_l * betaSA * xi_s_i  \
                 + HAS * betaAS * (FAS - FAS_n) + HAS * betaAS * eps_l * E * (1 - xi_s_i)

            iter=iter+1
            print("     NR iteration   ", iter," has residual ", R )

       # xi_s_i = min(1,xi_s_i)
        xi_s[i+1] = xi_s_i
        print("     The converged value of stress is ", sig[i+1], "and the value of the internal variable is  ", xi_s_i)
        print("              ")

plt.plot(eps, sig)
plt.show()
