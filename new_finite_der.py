import subprocess
from numpy import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scipy as spy
from scipy import spatial
from itertools import product
from numpy import save

numz = 301
#fiducial parameter values as given in table 3 of Planck 2015 (arxiv 1502.01589)
 #updated accroding to table 1 of Planck 2018 (arXiv 1807.06209) -- Yang
fv = [0.02233,0.1198,0.0001,10**-27,67.37,0.9652]
#parameter confidence intervals as given in table 3 of Planck 2015 (arxiv 1502.01589)
 #updated accroding to table 1 of Planck 2018 (arXiv 1807.06209) -- Yang
ci = [0.00015,0.0012,'n/a','n/a',0.54,0.0042]
#number of axioncamb runs for each parameter
# runs = 25

#testrun with fiducial values
subprocess.call('../axionCAMB-master/camb params_backup.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0]*(1+0.02),fv[1],fv[2],fv[3],fv[4],fv[5]), shell = True)
subprocess.call('../axionCAMB-master/camb params.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0]*(1+0.02),fv[1],fv[2],fv[3],fv[4],fv[5]), shell = True)




omegab = np.round(np.linspace(fv[0]-10*ci[0],fv[0]+10*ci[0],runs),5)
omegac = np.round(np.linspace(fv[1]-10*ci[1],fv[1]+10*ci[1],runs),4)
fa = np.linspace(0.0001,0.025,runs)
omegaax = fv[1]*fa
omegac2 = fv[1]-omegaax
fa_zip = zip(omegac2,omegaax)
h = np.round(np.linspace(fv[4]-ci[4],fv[4]+ci[4],runs),2)
ns = np.round(np.linspace(fv[5]-ci[5],fv[5]+ci[5],runs),4)

epsilon = np.linspace(0.002,0.008,8)

# print(len(omegab))
# print(len(omegac))


# omegab ##
for i in range(runs):
for e in range(len(epsilon)):
    subprocess.call('../axionCAMB-master/camb params_backup.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0]*(1+epsilon[e]),fv[1],fv[2],fv[3],fv[4],fv[5]), shell = True)
    subprocess.call('../axionCAMB-master/camb params.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0]*(1+epsilon[e]),fv[1],fv[2],fv[3],fv[4],fv[5]), shell = True)
    for k in range(numz):
        subprocess.call('mv test_matterpower%s.dat ../Axion-main/new_derivative/omegab_%s_upstep%s.dat' %(k,k,e), shell=True)

print('partly done')

# for i in range(runs):
for e in range(len(epsilon)):
    subprocess.call('../axionCAMB-master/camb params_backup.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0]*(1-epsilon[e]),fv[1],fv[2],fv[3],fv[4],fv[5]), shell = True)
    subprocess.call('../axionCAMB-master/camb params.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0]*(1-epsilon[e]),fv[1],fv[2],fv[3],fv[4],fv[5]), shell = True)
    for k in range(numz):
        subprocess.call('mv test_matterpower%s.dat ../Axion-main/new_derivative/omegab_%s_downstep%s.dat' %(k,k,e), shell=True)

print('done generating files')

# der = []
# for i in range(runs):
#     l1 = []
#     for e in range(len(epsilon)):
#         l2 = []
#         for k in range(numz):
#             diffder1 = loadtxt('../Axion-main/new_derivative/omegab_%s_%s_upstep%s.dat' %(i,k,e),unpack=True,usecols=[1])
#             diffder2 = loadtxt('../Axion-main/new_derivative/omegab_%s_%s_downstep%s.dat' %(i,k,e),unpack=True,usecols=[1])
#             diffder = (diffder1-diffder2)/(2*epsilon[e]*omegab[i])
#             l2.append(diffder)
#         l1.append(l2)
#     der.append(l1)
#
# print(np.shape(der))
# save('../Axion-main/new_derivative/der_omegab.npy',der)
#
# print('done calculating derivatives')


# omegac ##
epsilon = np.linspace(0.002,0.008,8)

# for i in range(runs):
for e in range(len(epsilon)):
    subprocess.call('../axionCAMB-master/camb params_backup.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1]*(1+epsilon[e]),fv[2],fv[3],fv[4],fv[5]), shell = True)
    subprocess.call('../axionCAMB-master/camb params.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1]*(1+epsilon[e]),fv[2],fv[3],fv[4],fv[5]), shell = True)
    for k in range(numz):
        subprocess.call('mv test_matterpower%s.dat ../Axion-main/new_derivative/omegac_%s_upstep%s.dat' %(k,k,e), shell=True)

print('partly done')

# for i in range(runs):
for e in range(len(epsilon)):
    subprocess.call('../axionCAMB-master/camb params_backup.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1]*(1+epsilon[e]),fv[2],fv[3],fv[4],fv[5]), shell = True)
    subprocess.call('../axionCAMB-master/camb params.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1]*(1-epsilon[e]),fv[2],fv[3],fv[4],fv[5]), shell = True)
    for k in range(numz):
        subprocess.call('mv test_matterpower%s.dat ../Axion-main/new_derivative/omegac_%s_downstep%s.dat' %(k,k,e), shell=True)

print('done generating files')

# der = []
# for i in range(runs):
#     l1 = []
#     for e in range(len(epsilon)):
#         l2 = []
#         for k in range(numz):
#             diffder1 = loadtxt('../Axion-main/new_derivative/omegac_%s_%s_upstep%s.dat'%(i,k,e),unpack=True,usecols=[1])
#             diffder2 = loadtxt('../Axion-main/new_derivative/omegac_%s_%s_downstep%s.dat'%(i,k,e),unpack=True,usecols=[1])
#             diffder = (diffder1-diffder2)/(2*epsilon[e]*omegac[i])
#             l2.append(diffder)
#         l1.append(l2)
#     der.append(l1)
#
# print(np.shape(der))
# save('../Axion-main/new_derivative/der_omegac.npy',der)

# print('done calculating derivatives')



# fa ##
# accuracy boost, sample boost, lsample boost in ini file
# 1-side der, use smaller step size but with higher percision
# 2-side der, use fraction step size
epsilon = np.linspace(10,30,8)

# for i in range(runs):
for e in range(len(epsilon)):
    subprocess.call('../axionCAMB-master/camb params_backup.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2]*(1+epsilon[e]),fv[3],fv[4],fv[5]), shell = True)
    subprocess.call('../axionCAMB-master/camb params.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2]*(1+epsilon[e]),fv[3],fv[4],fv[5]), shell = True)
    for k in range(numz):
        subprocess.call('mv test_matterpower%s.dat ../Axion-main/new_derivative/fa_%s_upstep%s.dat' %(k,k,e), shell=True)

print('partly done')

# for i in range(runs):
for e in range(len(epsilon)):
    subprocess.call('../axionCAMB-master/camb params_backup.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2],fv[3],fv[4],fv[5]), shell = True)
    subprocess.call('../axionCAMB-master/camb params.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2],fv[3],fv[4],fv[5]), shell = True)
    for k in range(numz):
        subprocess.call('mv test_matterpower%s.dat ../Axion-main/new_derivative/fa_%s_downstep%s.dat' %(k,k,e), shell=True)

print('done generating files')

# der = []
# for i in range(runs):
#     l1 = []
#     for e in range(len(epsilon)):
#         l2 = []
#         for k in range(numz):
#             diffder1 = loadtxt('../Axion-main/new_derivative/fa_%s_%s_upstep%s.dat'%(i,k,e),unpack=True,usecols=[1])
#             diffder2 = loadtxt('../Axion-main/new_derivative/fa_%s_%s_downstep%s.dat'%(i,k,e),unpack=True,usecols=[1])
#             diffder = (diffder1-diffder2)/(2*epsilon[e]*fa[i])
#             l2.append(diffder)
#         l1.append(l2)
#     der.append(l1)
#
# print(np.shape(der))
# save('../Axion-main/new_derivative/der_fa.npy',der)
#
# print('done calculating derivatives')



# ns ##
epsilon = np.linspace(0.005,0.04,8)

# for i in range(runs):
for e in range(len(epsilon)):
    subprocess.call('../axionCAMB-master/camb params_backup.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2],fv[3],fv[4],fv[5]*(1+epsilon[e])), shell = True)
    subprocess.call('../axionCAMB-master/camb params.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2],fv[3],fv[4],fv[5]*(1+epsilon[e])), shell = True)
    for k in range(numz):
        subprocess.call('mv test_matterpower%s.dat ../Axion-main/new_derivative/ns_%s_upstep%s.dat' %(k,k,e), shell=True)

print('partly done')

# for i in range(runs):
for e in range(len(epsilon)):
    subprocess.call('../axionCAMB-master/camb params_backup.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2],fv[3],fv[4],fv[5]*(1-epsilon[e])), shell = True)
    subprocess.call('../axionCAMB-master/camb params.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2],fv[3],fv[4],fv[5]*(1-epsilon[e])), shell = True)
    for k in range(numz):
        subprocess.call('mv test_matterpower%s.dat ../Axion-main/new_derivative/ns_%s_downstep%s.dat' %(k,k,e), shell=True)

print('done generating files')

# der = []
# for i in range(runs):
#     l1 = []
#     for e in range(len(epsilon)):
#         l2 = []
#         for k in range(numz):
#             diffder1 = loadtxt('../Axion-main/new_derivative/ns_%s_%s_upstep%s.dat'%(i,k,e),unpack=True,usecols=[1])
#             diffder2 = loadtxt('../Axion-main/new_derivative/ns_%s_%s_downstep%s.dat'%(i,k,e),unpack=True,usecols=[1])
#             diffder = (diffder1-diffder2)/(2*epsilon[e]*ns[i])
#             l2.append(diffder)
#         l1.append(l2)
#     der.append(l1)
#
# print(np.shape(der))
# save('../Axion-main/new_derivative/der_ns.npy',der)
#
# print('done calculating derivatives')



## h ##
epsilon = np.linspace(0.01,0.08,8)

# for i in range(runs):
for e in range(len(epsilon)):
    subprocess.call('../axionCAMB-master/camb params_backup.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2],fv[3],fv[4]*(1+epsilon[e]),fv[5]), shell = True)
    subprocess.call('../axionCAMB-master/camb params.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2],fv[3],fv[4]*(1+epsilon[e]),fv[5]), shell = True)
    for k in range(numz):
        subprocess.call('mv test_matterpower%s.dat ../Axion-main/new_derivative/h_%s_upstep%s.dat' %(k,k,e), shell=True)

print('partly done')

# for i in range(runs):
for e in range(len(epsilon)):
    subprocess.call('../axionCAMB-master/camb params_backup.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2],fv[3],fv[4]*(1-epsilon[e]),fv[5]), shell = True)
    subprocess.call('../axionCAMB-master/camb params.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2],fv[3],fv[4]*(1-epsilon[e]),fv[5]), shell = True)
    for k in range(numz):
        subprocess.call('mv test_matterpower%s.dat ../Axion-main/new_derivative/h_%s_downstep%s.dat' %(k,k,e), shell=True)

print('done generating files')

# der = []
# for i in range(runs):
#     l1 = []
#     for e in range(len(epsilon)):
#         l2 = []
#         for k in range(numz):
#             diffder1 = loadtxt('../Axion-main/new_derivative/h_%s_%s_upstep%s.dat'%(i,k,e),unpack=True,usecols=[1])
#             diffder2 = loadtxt('../Axion-main/new_derivative/h_%s_%s_downstep%s.dat'%(i,k,e),unpack=True,usecols=[1])
#             diffder = (diffder1-diffder2)/(2*epsilon[e]*h[i])
#             l2.append(diffder)
#         l1.append(l2)
#     der.append(l1)
#
# print(np.shape(der))
# save('../Axion-main/new_derivative/der_h.npy',der)

print('done calculating derivatives')
