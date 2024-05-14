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

# omegab = np.round(np.linspace(fv[0]-10*ci[0],fv[0]+10*ci[0],runs),5)
# omegac = np.round(np.linspace(fv[1]-10*ci[1],fv[1]+10*ci[1],runs),4)
# fa = np.linspace(0.0001,0.025,runs)
# omegaax = fv[1]*fa
# omegac2 = fv[1]-omegaax
# fa_zip = zip(omegac2,omegaax)
# h = np.round(np.linspace(fv[4]-ci[4],fv[4]+ci[4],runs),2)
# ns = np.round(np.linspace(fv[5]-ci[5],fv[5]+ci[5],runs),4)


## omegab ##

epsilon = np.linspace(0.002,0.008,8)

mergekd = []
mergepkd = []

mergeku = []
mergepku = []

# for i in range(runs):
k1 = []
pk1 = []
for e in range(len(epsilon)):
    k2 = []
    pk2 = []
    for y in range(numz):
        kd = loadtxt('../Axion-main/new_derivative/omegab_%s_downstep%s.dat'%(y,7-e),unpack=True,usecols=[0])
        pkd = file = loadtxt('../Axion-main/new_derivative/omegab_%s_downstep%s.dat'%(y,7-e),unpack=True,usecols=[1])
        k2.append(kd)
        pk2.append(pkd)
    k1.append(k2)
    pk1.append(pk2)
mergekd.append(k1)
mergepkd.append(pk1)


# for i in range(runs):
ku1 = []
pku1 = []
for e in range(len(epsilon)):
    ku2 = []
    pku2 = []
    for y in range(numz):
        ku = loadtxt('../Axion-main/new_derivative/omegab_%s_upstep%s.dat'%(y,e),unpack=True,usecols=[0])
        pku = loadtxt('../Axion-main/new_derivative/omegab_%s_upstep%s.dat'%(y,e),unpack=True,usecols=[1])
        ku2.append(ku)
        pku2.append(pku)
    ku1.append(ku2)
    pku1.append(pku2)
mergeku.append(ku1)
mergepku.append(pku1)


print(np.shape(mergekd))
print(np.shape(mergeku))
# save('../Axion-main/ns_mergek2.npy',mergek)
# save('../Axion-main/ns_mergepk2.npy',mergepk)
fullk = np.concatenate((np.array(mergekd),np.array(mergeku)),axis=1)
fullpk = np.concatenate((np.array(mergepkd),np.array(mergepku)),axis=1)
print(np.shape(fullk))
print(np.shape(fullpk))
save('../Axion-main/new_omegab_k.npy',fullk)
save('../Axion-main/new_omegab_pk.npy',fullpk)
print('done creating k and pk omegab')



der = []
# for i in range(runs):
l1 = []
for e in range(len(epsilon)):
    l2 = []
    for k in range(numz):
        diffder1 = loadtxt('../Axion-main/new_derivative/omegab_%s_upstep%s.dat'%(k,e),unpack=True,usecols=[1])
        diffder2 = loadtxt('../Axion-main/new_derivative/omegab_%s_downstep%s.dat'%(k,e),unpack=True,usecols=[1])
        diffder = (diffder1-diffder2)/(2*epsilon[e]*fv[0])
        l2.append(diffder)
    l1.append(l2)
der.append(l1)

print(np.shape(der))
save('../Axion-main/new_omegab_derivative.npy',der)

print('done calculating derivatives of omegab')






## omegac ##

epsilon = np.linspace(0.002,0.008,8)


mergekd = []
mergepkd = []

mergeku = []
mergepku = []

# for i in range(runs):
k1 = []
pk1 = []
for e in range(len(epsilon)):
    k2 = []
    pk2 = []
    for y in range(numz):
        kd = loadtxt('../Axion-main/new_derivative/omegac_%s_downstep%s.dat'%(y,7-e),unpack=True,usecols=[0])
        pkd = file = loadtxt('../Axion-main/new_derivative/omegac_%s_downstep%s.dat'%(y,7-e),unpack=True,usecols=[1])
        k2.append(kd)
        pk2.append(pkd)
    k1.append(k2)
    pk1.append(pk2)
mergekd.append(k1)
mergepkd.append(pk1)


# for i in range(runs):
ku1 = []
pku1 = []
for e in range(len(epsilon)):
    ku2 = []
    pku2 = []
    for y in range(numz):
        ku = loadtxt('../Axion-main/new_derivative/omegac_%s_upstep%s.dat'%(y,e),unpack=True,usecols=[0])
        pku = loadtxt('../Axion-main/new_derivative/omegac_%s_upstep%s.dat'%(y,e),unpack=True,usecols=[1])
        ku2.append(ku)
        pku2.append(pku)
    ku1.append(ku2)
    pku1.append(pku2)
mergeku.append(ku1)
mergepku.append(pku1)


print(np.shape(mergekd))
print(np.shape(mergeku))
# save('../Axion-main/ns_mergek2.npy',mergek)
# save('../Axion-main/ns_mergepk2.npy',mergepk)
fullk = np.concatenate((np.array(mergekd),np.array(mergeku)),axis=1)
fullpk = np.concatenate((np.array(mergepkd),np.array(mergepku)),axis=1)
print(np.shape(fullk))
print(np.shape(fullpk))
save('../Axion-main/new_omegac_k.npy',fullk)
save('../Axion-main/new_omegac_pk.npy',fullpk)
print('done creating k and pk omegac')



der = []
# for i in range(runs):
l1 = []
for e in range(len(epsilon)):
    l2 = []
    for k in range(numz):
        diffder1 = loadtxt('../Axion-main/new_derivative/omegac_%s_upstep%s.dat'%(k,e),unpack=True,usecols=[1])
        diffder2 = loadtxt('../Axion-main/new_derivative/omegac_%s_downstep%s.dat'%(k,e),unpack=True,usecols=[1])
        diffder = (diffder1-diffder2)/(2*epsilon[e]*fv[1])
        l2.append(diffder)
    l1.append(l2)
der.append(l1)

print(np.shape(der))
save('../Axion-main/new_omegac_derivative.npy',der)

print('done calculating derivatives of omegac')




## fa ##

epsilon = np.linspace(10,30,8)


mergekd = []
mergepkd = []

mergeku = []
mergepku = []

# for i in range(runs):
k1 = []
pk1 = []
for e in range(len(epsilon)):
    k2 = []
    pk2 = []
    for y in range(numz):
        kd = loadtxt('../Axion-main/new_derivative/fa_%s_downstep%s.dat'%(y,7-e),unpack=True,usecols=[0])
        pkd = file = loadtxt('../Axion-main/new_derivative/fa_%s_downstep%s.dat'%(y,7-e),unpack=True,usecols=[1])
        k2.append(kd)
        pk2.append(pkd)
    k1.append(k2)
    pk1.append(pk2)
mergekd.append(k1)
mergepkd.append(pk1)


# for i in range(runs):
ku1 = []
pku1 = []
for e in range(len(epsilon)):
    ku2 = []
    pku2 = []
    for y in range(numz):
        ku = loadtxt('../Axion-main/new_derivative/fa_%s_upstep%s.dat'%(y,e),unpack=True,usecols=[0])
        pku = loadtxt('../Axion-main/new_derivative/fa_%s_upstep%s.dat'%(y,e),unpack=True,usecols=[1])
        ku2.append(ku)
        pku2.append(pku)
    ku1.append(ku2)
    pku1.append(pku2)
mergeku.append(ku1)
mergepku.append(pku1)


print(np.shape(mergekd))
print(np.shape(mergeku))
# save('../Axion-main/ns_mergek2.npy',mergek)
# save('../Axion-main/ns_mergepk2.npy',mergepk)
fullk = np.concatenate((np.array(mergekd),np.array(mergeku)),axis=1)
fullpk = np.concatenate((np.array(mergepkd),np.array(mergepku)),axis=1)
print(np.shape(fullk))
print(np.shape(fullpk))
save('../Axion-main/new_fa_k.npy',fullk)
save('../Axion-main/new_fa_pk.npy',fullpk)
print('done creating k and pk fa')



der = []
# for i in range(runs):
l1 = []
for e in range(len(epsilon)):
    l2 = []
    for k in range(numz):
        diffder1 = loadtxt('../Axion-main/new_derivative/fa_%s_upstep%s.dat'%(k,e),unpack=True,usecols=[1])
        diffder2 = loadtxt('../Axion-main/new_derivative/fa_%s_downstep%s.dat'%(k,e),unpack=True,usecols=[1])
        diffder = (diffder1-diffder2)/(2*epsilon[e]*fv[2])
        l2.append(diffder)
    l1.append(l2)
der.append(l1)

print(np.shape(der))
save('../Axion-main/new_fa_derivative.npy',der)

print('done calculating derivatives of fa')




## ns ##

epsilon = np.linspace(0.005,0.04,8)


mergekd = []
mergepkd = []

mergeku = []
mergepku = []

# for i in range(runs):
k1 = []
pk1 = []
for e in range(len(epsilon)):
    k2 = []
    pk2 = []
    for y in range(numz):
        kd = loadtxt('../Axion-main/new_derivative/ns_%s_downstep%s.dat'%(y,7-e),unpack=True,usecols=[0])
        pkd = file = loadtxt('../Axion-main/new_derivative/ns_%s_downstep%s.dat'%(y,7-e),unpack=True,usecols=[1])
        k2.append(kd)
        pk2.append(pkd)
    k1.append(k2)
    pk1.append(pk2)
mergekd.append(k1)
mergepkd.append(pk1)


# for i in range(runs):
ku1 = []
pku1 = []
for e in range(len(epsilon)):
    ku2 = []
    pku2 = []
    for y in range(numz):
        ku = loadtxt('../Axion-main/new_derivative/ns_%s_upstep%s.dat'%(y,e),unpack=True,usecols=[0])
        pku = loadtxt('../Axion-main/new_derivative/ns_%s_upstep%s.dat'%(y,e),unpack=True,usecols=[1])
        ku2.append(ku)
        pku2.append(pku)
    ku1.append(ku2)
    pku1.append(pku2)
mergeku.append(ku1)
mergepku.append(pku1)


print(np.shape(mergekd))
print(np.shape(mergeku))
# save('../Axion-main/ns_mergek2.npy',mergek)
# save('../Axion-main/ns_mergepk2.npy',mergepk)
fullk = np.concatenate((np.array(mergekd),np.array(mergeku)),axis=1)
fullpk = np.concatenate((np.array(mergepkd),np.array(mergepku)),axis=1)
print(np.shape(fullk))
print(np.shape(fullpk))
save('../Axion-main/new_ns_k.npy',fullk)
save('../Axion-main/new_ns_pk.npy',fullpk)
print('done creating k and pk ns')



der = []
# for i in range(runs):
l1 = []
for e in range(len(epsilon)):
    l2 = []
    for k in range(numz):
        diffder1 = loadtxt('../Axion-main/new_derivative/ns_%s_upstep%s.dat'%(k,e),unpack=True,usecols=[1])
        diffder2 = loadtxt('../Axion-main/new_derivative/ns_%s_downstep%s.dat'%(k,e),unpack=True,usecols=[1])
        diffder = (diffder1-diffder2)/(2*epsilon[e]*fv[5])
        l2.append(diffder)
    l1.append(l2)
der.append(l1)

print(np.shape(der))
save('../Axion-main/new_ns_derivative.npy',der)

print('done calculating derivatives of ns')




## h ##


mergekd = []
mergepkd = []

mergeku = []
mergepku = []

# for i in range(runs):
k1 = []
pk1 = []
for e in range(len(epsilon)):
    k2 = []
    pk2 = []
    for y in range(numz):
        kd = loadtxt('../Axion-main/new_derivative/h_%s_downstep%s.dat'%(y,7-e),unpack=True,usecols=[0])
        pkd = file = loadtxt('../Axion-main/new_derivative/h_%s_downstep%s.dat'%(y,7-e),unpack=True,usecols=[1])
        k2.append(kd)
        pk2.append(pkd)
    k1.append(k2)
    pk1.append(pk2)
mergekd.append(k1)
mergepkd.append(pk1)


# for i in range(runs):
ku1 = []
pku1 = []
for e in range(len(epsilon)):
    ku2 = []
    pku2 = []
    for y in range(numz):
        ku = loadtxt('../Axion-main/new_derivative/h_%s_upstep%s.dat'%(y,e),unpack=True,usecols=[0])
        pku = loadtxt('../Axion-main/new_derivative/h_%s_upstep%s.dat'%(y,e),unpack=True,usecols=[1])
        ku2.append(ku)
        pku2.append(pku)
    ku1.append(ku2)
    pku1.append(pku2)
mergeku.append(ku1)
mergepku.append(pku1)


print(np.shape(mergekd))
print(np.shape(mergepku))
print(np.shape(mergeku))
# save('../Axion-main/ns_mergek2.npy',mergek)
# save('../Axion-main/ns_mergepk2.npy',mergepk)
fullk = np.concatenate((np.array(mergekd),np.array(mergeku)),axis=1)
fullpk = np.concatenate((np.array(mergepkd),np.array(mergepku)),axis=1)
print(np.shape(fullk))
print(np.shape(fullpk))
save('../Axion-main/new_h_k.npy',fullk)
save('../Axion-main/new_h_pk.npy',fullpk)
print('done creating k and pk h')



epsilon = np.linspace(0.01,0.08,8)

der = []
# for i in range(runs):
l1 = []
for e in range(len(epsilon)):
    l2 = []
    for k in range(numz):
        diffder1 = loadtxt('../Axion-main/new_derivative/h_%s_upstep%s.dat'%(k,e),unpack=True,usecols=[1])
        diffder2 = loadtxt('../Axion-main/new_derivative/h_%s_downstep%s.dat'%(k,e),unpack=True,usecols=[1])
        diffder = (diffder1-diffder2)/(2*epsilon[e]*fv[4])
        l2.append(diffder)
    l1.append(l2)
der.append(l1)

print(np.shape(der))
save('../Axion-main/new_h_derivative.npy',der)

print('done calculating derivatives of h')
