import subprocess
from numpy import *
import numpy as np
import scipy as spy
from scipy import spatial
from scipy import interpolate
from numpy import save

import matplotlib
import matplotlib.pyplot as plt

numz = 301

fv = [0.02233,0.1198,0.0001,10**-27,67.37,0.9652]
ci = [0.00015,0.0012,'n/a','n/a',0.54,0.0042]

runs = 1



## omegab ##
epsilon = np.linspace(0.002,0.008,8)
for i in range(runs):
    subprocess.call('../axionCAMB-master/camb params_backup.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2],fv[3],fv[4],fv[5]), shell = True)
    # subprocess.call('../axionCAMB-master/camb params.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2],fv[3],fv[4],fv[5]), shell = True)
    # subprocess.call('../axionCAMB-master/camb params_change.ini 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2],fv[3],fv[4],fv[5]), shell = True)
    # subprocess.call('../axionCAMB-master/axion_background.o 1 %s 2 %s 3 %s 4 %s 5 %s 6 %s' %(fv[0],fv[1],fv[2],fv[3],fv[4],fv[5]), shell = True)



aosc_params = loadtxt('aosc_param.dat', unpack=True)
# for i in range(10):
#     param = loadtxt('aosc_param.dat', unpack=True, usecols=[i])
#     aosc_params.append(param)

#print(aosc_params)
# print(np.shape(aosc_params[6]))

# v_1 = aosc_params[0]
# v_2 = aosc_params[1]

# omegah2_regm, Params%omegah2_rad, omegah2_lambda, &
#      maxion_twiddle, omk*hsq, a_arr(i), v_vec(1:2,i), lhsqcont_massive, Params%Nu_mass_eigenstates, mass_correctors
omegah2_regm = aosc_params[2]
print('matter',omegah2_regm)
omegah2_rad = aosc_params[3]
print('radiation',omegah2_rad)
omegah2_lambda = aosc_params[4]
print('dark energy',omegah2_lambda)
maxion_twiddle = aosc_params[5]
omkhsq = aosc_params[6]
# a = np.exp(aosc_params[7])
a = aosc_params[7]
print('a', a)
lhsqcont_massive = aosc_params[8]
Nu_mass_eigenstates = aosc_params[9]
mass_correctors = aosc_params[10]
v_1 = loadtxt('littleh.dat', unpack=True, usecols=[2])
v_1 = v_1[-1]
print('v1', v_1)
v_2 = loadtxt('littleh.dat', unpack=True, usecols=[3])
v_2 = v_2[-1]
print('v2', v_2)


littlehfunc=omegah2_rad/(a**4.0)+(omegah2_regm/(a**3.0))+(lhsqcont_massive*mass_correctors+Nu_mass_eigenstates)/(a**4.0)
littlehfunc=littlehfunc+omegah2_lambda+((v_2/a)**2.0)
littlehfunc=littlehfunc+(maxion_twiddle*v_1)**2.0
littlehfunc=littlehfunc*(a**2.0)+omkhsq
littlehfunc=sqrt(littlehfunc)
print(((v_2/a)**2.0)+(maxion_twiddle*v_1)**2.0)
print('littlehfunc',littlehfunc)




#
# aosc_params = []
# with open('aosc_param.dat', 'r') as f:
#     d = f.readlines()
#     for i in d:
#         k = i.rstrip().split(",")
#         aosc_params.append([float(i) if is_float(i) else i for i in k])
#
# print(np.shape(aosc_params))


littleh = loadtxt('littleh.dat', unpack=True, usecols=[1])
scalefactor = loadtxt('littleh.dat', unpack=True, usecols=[0])
# print(scalefactor[0],scalefactor[4998])
redshift = 1/scalefactor - 1
# print('h min:',min(littleh))
# print('h max:',max(littleh))
# print('a min:',min(scalefactor))
# print('a max:',max(scalefactor))
# print('z min:',min(redshift))
# print('z max:',max(redshift))
# print("shape of littleh array",np.shape(littleh))
print(littleh[0],littleh[4998])
# littleh = np.array(littleh, dtype='O')
# littleh = np.genfromtxt('littleh.dat',
#                      dtype=None,
#                      delimiter=',')
# print(littleh)
# print(scalefactor)
# print(type(littleh))
# print(np.shape(littleh))
'''



zlist = np.linspace(0,4.5,num=numz)
# zlist = zlist[::-1]
# print(zlist)

hubble = littleh*100

print(np.shape(hubble))
print(np.shape(redshift))


f = interpolate.interp1d(redshift, hubble, fill_value='extrapolate')
H = f(zlist)
print('H min:',min(H))
print('H max:',max(H))
# print(f(5000))
# print(f(zlist))
'''



'''
# plt.plot(redshift,hubble)
plt.plot(zlist,f(zlist))

# plt.yscale('log')
plt.ylabel('H')
plt.xlabel('z')
plt.title('H interpolated through the same z values for other parameters')
# plt.gca().invert_xaxis()
plt.savefig('H_check.pdf')
'''
