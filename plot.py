import numpy as np
import time
import os
import matplotlib
from matplotlib import pyplot as plt


#set fonts
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'serif','serif':['Palatino']})


k = np.loadtxt('nonlin_PS_with_axion_full.txt',unpack=True, usecols=[0])
Pk = np.loadtxt('nonlin_PS_with_axion_full.txt',unpack=True, usecols=[2])
plt.plot(k,Pk, label=r'$\Lambda$CDM cosmology')


k = np.loadtxt('nonlin_PS_with_axion_full.txt',unpack=True, usecols=[0])
Pk = np.loadtxt('nonlin_PS_with_axion_full.txt',unpack=True, usecols=[1])
plt.plot(k,Pk, label=r'$f_a=100\%$ $m_a=10^{-21}$ eV', linestyle='dashed')


k = np.loadtxt('nonlin_PS_with_axion_full_3.txt',unpack=True, usecols=[0])
Pk = np.loadtxt('nonlin_PS_with_axion_full_3.txt',unpack=True, usecols=[1])
plt.plot(k,Pk, label=r'$f_a=100\%$ $m_a=10^{-23}$ eV', linestyle='dashed')


k = np.loadtxt('nonlin_PS_with_axion_1.txt',unpack=True, usecols=[0])
Pk = np.loadtxt('nonlin_PS_with_axion_1.txt',unpack=True, usecols=[1])
plt.plot(k,Pk, label=r'$f_a=10\%$ $m_a=10^{-21}$ eV',linestyle='dashdot')


k = np.loadtxt('nonlin_PS_with_axion_3.txt',unpack=True, usecols=[0])
Pk = np.loadtxt('nonlin_PS_with_axion_3.txt',unpack=True, usecols=[1])
plt.plot(k,Pk, label=r'$f_a=10\%$ $m_a=10^{-23}$ eV',linestyle='dashdot')


k = np.loadtxt('nonlin_PS_with_axion_8.txt',unpack=True, usecols=[0])
Pk = np.loadtxt('nonlin_PS_with_axion_8.txt',unpack=True, usecols=[1])
plt.plot(k,Pk, label=r'$f_a=10\%$ $m_a=10^{-28}$ eV',linestyle='dashdot')
# for i in range(2):
#     j = i+2
#     ma = j+20
#
#     k = np.loadtxt('nonlin_PS_with_axion_full_%s.txt' %(j), unpack=True, usecols=[0])
#     Pk = np.loadtxt('nonlin_PS_with_axion_full_%s.txt' %(j), unpack=True, usecols=[1])
#
#     plt.plot(k,Pk, label=r'$m_a=10^{-%s}$ eV' %(ma))
#     print(j, 'done')

# print(np.shape(k),np.shape(Pk))



# plt.plot(power_spec_dic_ax['k'],PS_LCDM_matter_nonlin[0], label=r'non-linear power spectrum in $\Lambda$CDM cosmology')
plt.xscale('log')
plt.yscale('log')
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel(r'k',fontsize=24)
plt.ylabel(r'P(k)',fontsize=24)
plt.title(r'Non-linear Power Spectra',fontsize=30)
plt.legend(fontsize='medium')
plt.savefig('Pk_poster.pdf',bbox_inches='tight')
plt.show()
