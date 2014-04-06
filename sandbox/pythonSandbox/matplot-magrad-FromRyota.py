import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pylab
import math
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from scipy import genfromtxt
import matplotlib as mpl
mpl.rc('figure.subplot',bottom=0.12)

params = {'axes.labelsize': 20,
    'text.fontsize': 10,
    'legend.fontsize': 10,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20
#,'figure.figsize': [8, 6.0]
}
#    'text.usetex': True

pylab.rcParams.update(params)

fig = plt.figure()
ax = fig.add_subplot(111)
#print   fig.subplotpars.top
#print   fig.subplotpars.wspace

####################################################################

# Size - Luminosity _z8

data7 = np.loadtxt('/Users/KawamataRyota/Dropbox/Reports@Dropbox/Gravitational_Lensing@Dropbox/glafic/Abell2744/Plot_Abell/catalogue-dropout_z7_wo-average.cat', comments='#')
mag7=data7[:,74]
rad7=data7[:,73]
sigmasfr7=data7[:,84]
beta7= -2.0 + 4.39*(data7[:,16]-data7[:,30])
a1600_7 = 4.43 + beta7*1.99
sfr7=data7[:,81]
ellip7=data7[:,87-1]
sigmasfr7_corr=data7[:,84] * 10**(a1600_7/2.5)


xerr_u7 = data7[:,79]
xerr_l7 = data7[:,78]
yerr_u7 = data7[:,77]
yerr_l7 = data7[:,76]


data7_b1 = np.loadtxt('/Users/KawamataRyota/Dropbox/Reports@Dropbox/Gravitational_Lensing@Dropbox/glafic/Abell2744/Plot_Abell/catalogue_z7-b1.cat', comments='#')
mag7_b1=data7_b1[:,74]
rad7_b1=data7_b1[:,73]
sigmasfr7_b1=data7_b1[:,84]
beta7_b1= -2.0 + 4.39*(data7_b1[:,16]-data7_b1[:,30])
sfr7_b1=data7_b1[:,81]
ellip7_b1=data7_b1[:,87-1]

xerr_u7_b1 = data7_b1[:,79]
xerr_l7_b1 = data7_b1[:,78]
yerr_u7_b1 = data7_b1[:,77]
yerr_l7_b1 = data7_b1[:,76]



data7_b2 = np.loadtxt('/Users/KawamataRyota/Dropbox/Reports@Dropbox/Gravitational_Lensing@Dropbox/glafic/Abell2744/Plot_Abell/catalogue_z7-b2.cat', comments='#')
mag7_b2=data7_b2[:,74]
rad7_b2=data7_b2[:,73]
sigmasfr7_b2=data7_b2[:,84]
beta7_b2= -2.0 + 4.39*(data7_b2[:,16]-data7_b2[:,30])
sfr7_b2=data7_b2[:,81]
ellip7_b2=data7_b2[:,87-1]

xerr_u7_b2 = data7_b2[:,79]
xerr_l7_b2 = data7_b2[:,78]
yerr_u7_b2 = data7_b2[:,77]
yerr_l7_b2 = data7_b2[:,76]

data8_b1 = np.loadtxt('/Users/KawamataRyota/Dropbox/Reports@Dropbox/Gravitational_Lensing@Dropbox/glafic/Abell2744/Plot_Abell/catalogue_z8-b1.cat', comments='#')
mag8_b1=data8_b1[:,74]
rad8_b1=data8_b1[:,73]
sigmasfr8_b1=data8_b1[:,84]
beta8_b1= -2.0 + 4.39*(data8_b1[:,16]-data8_b1[:,30])
sfr8_b1=data8_b1[:,81]
ellip8_b1=data8_b1[:,87-1]

xerr_u8_b1 = data8_b1[:,79]
xerr_l8_b1 = data8_b1[:,78]
yerr_u8_b1 = data8_b1[:,77]
yerr_l8_b1 = data8_b1[:,76]



#data8_b2 = np.loadtxt('/Users/KawamataRyota/Dropbox/Reports@Dropbox/Gravitational_Lensing@Dropbox/glafic/Abell2744/Plot_Abell/catalogue_z8-b2.cat', comments='#')
#mag8_b2=data8_b2[:,74]
#rad8_b2=data8_b2[:,73]
#sigmasfr8_b2=data8_b2[:,84]
#beta8_b2= -2.0 + 4.39*(data8_b2[:,16]-data8_b2[:,30])
#sfr8_b2=data8_b2[:,81]
#ellip8_b2=data8_b2[:,87-1]

#xerr_u8_b2 = data8_b2[:,79]
#xerr_l8_b2 = data8_b2[:,78]
#yerr_u8_b2 = data8_b2[:,77]
#yerr_l8_b2 = data8_b2[:,76]


data8 = np.loadtxt('/Users/KawamataRyota/Dropbox/Reports@Dropbox/Gravitational_Lensing@Dropbox/glafic/Abell2744/Plot_Abell/catalogue-dropout_z8_wo-average.cat', comments='#')
mag8=data8[:,74]
rad8=data8[:,73]
sigmasfr8=data8[:,84]
sfr8=data8[:,81]
beta8= -2.0 + 8.98*(data8[:,24-1]-data8[:,31-1])
a1600_8 = 4.43 + beta8*1.99
ellip8=data8[:,87-1]
sigmasfr8_corr=data8[:,84] * 10**(a1600_8/2.5)


xerr_u8 = data8[:,79]
xerr_l8 = data8[:,78]
yerr_u8 = data8[:,77]
yerr_l8 = data8[:,76]


#data_multiple7 = np.loadtxt('multiple_catalogue_z7.cat', comments='#')
#size_multiple7 = data_multiple7[:,74-1]
#mag_multiple7 = data_multiple7[:,75-1]
#beta_multiple7=data_multiple7[:,17-1]-data_multiple7[:,31-1]


data_multiple8 = np.loadtxt('multiple_catalogue_z8.cat', comments='#')
size_multiple8 = data_multiple8[:,74-1]
mag_multiple8 = data_multiple8[:,75-1]
beta_multiple8= -2.0 + 8.98*(data_multiple8[:,24-1]-data_multiple8[:,31-1])

xerr_u_mag_multiple8 = data_multiple8[:,79]
xerr_l_mag_multiple8 = data_multiple8[:,78]
yerr_u_size_multiple8 = data_multiple8[:,77]
yerr_l_size_multiple8 = data_multiple8[:,76]

data_not_multiple8 = np.loadtxt('not-multiple_catalogue_z8.cat', comments='#')
size_not_multiple8 = data_not_multiple8[:,74-1]
mag_not_multiple8 = data_not_multiple8[:,75-1]
beta_not_multiple8= -2.0 + 8.98*(data_not_multiple8[:,24-1]-data_not_multiple8[:,31-1])

xerr_u_mag_not_multiple8 = data_not_multiple8[:,79]
xerr_l_mag_not_multiple8 = data_not_multiple8[:,78]
yerr_u_size_not_multiple8 = data_not_multiple8[:,77]
yerr_l_size_not_multiple8 = data_not_multiple8[:,76]


data_not_multiple7 = np.loadtxt('not-multiple_catalogue_z7.cat', comments='#')
size_not_multiple7 = data_not_multiple7[:,74-1]
mag_not_multiple7 = data_not_multiple7[:,75-1]
beta_not_multiple7= -2.0 + 4.39*(data_not_multiple7[:,17-1]-data_not_multiple7[:,31-1])

xerr_u_mag_not_multiple7 = data_not_multiple7[:,79]
xerr_l_mag_not_multiple7 = data_not_multiple7[:,78]
yerr_u_size_not_multiple7 = data_not_multiple7[:,77]
yerr_l_size_not_multiple7 = data_not_multiple7[:,76]




data_UDF = np.loadtxt('Muv_re/UDF12_ydropouts.data', comments='#')
x_UDF=data_UDF[:,0]
y_UDF=data_UDF[:,2]
raderr_UDF = data_UDF[:,1]
magerr_UDF = data_UDF[:,3]


data_UDF_b1 = np.loadtxt('Muv_re/UDF12_ydropouts_b1.data', comments='#')
x_UDF_b1=data_UDF_b1[0]
y_UDF_b1=data_UDF_b1[2]
raderr_UDF_b1 = data_UDF_b1[1]
magerr_UDF_b1 = data_UDF_b1[3]

data_UDF_b1_ave = np.loadtxt('Muv_re/UDF12_ydropouts_b1_ave.data', comments='#')
x_UDF_b1_ave=data_UDF_b1_ave[0]
y_UDF_b1_ave=data_UDF_b1_ave[2]
raderr_UDF_b1_ave = data_UDF_b1_ave[1]
magerr_UDF_b1_ave = data_UDF_b1_ave[3]


data_UDF_stack = np.loadtxt('Muv_re/UDF12_ydropouts_stack.data', comments='#')
x_UDF_stack=data_UDF_stack[:,0]
y_UDF_stack=data_UDF_stack[:,2]
raderr_UDF_stack = data_UDF_stack[:,1]
magerr_UDF_stack = data_UDF_stack[:,3]

data_multiple8_UDF = np.loadtxt('Muv_re/UDF12_ydropouts_multiple.data', comments='#')
size_multiple8_UDF = data_multiple8_UDF[2]
mag_multiple8_UDF = data_multiple8_UDF[0]

data_multiple7_UDF = np.loadtxt('Muv_re/UDF12_zdropouts_multiple.data', comments='#')
size_multiple7_UDF = data_multiple7_UDF[:,2]
mag_multiple7_UDF = data_multiple7_UDF[:,0]



data_UDF_z7 = np.loadtxt('Muv_re/UDF12_zdropouts.data', comments='#')
x_UDF_z7=data_UDF_z7[:,0]
y_UDF_z7=data_UDF_z7[:,2]
raderr_UDF_z7 = data_UDF_z7[:,1]
magerr_UDF_z7 = data_UDF_z7[:,3]



data_UDF_b1_z7 = np.loadtxt('Muv_re/UDF12_zdropouts_b1.data', comments='#')
x_UDF_b1_z7=data_UDF_b1_z7[:,0]
y_UDF_b1_z7=data_UDF_b1_z7[:,2]
raderr_UDF_b1_z7 = data_UDF_b1_z7[:,1]
magerr_UDF_b1_z7 = data_UDF_b1_z7[:,3]



data_UDF_b1_ave_z7 = np.loadtxt('Muv_re/UDF12_zdropouts_b1_ave.data', comments='#')
x_UDF_b1_ave_z7=data_UDF_b1_ave_z7[0]
y_UDF_b1_ave_z7=data_UDF_b1_ave_z7[2]
raderr_UDF_b1_ave_z7 = data_UDF_b1_ave_z7[1]
magerr_UDF_b1_ave_z7 = data_UDF_b1_ave_z7[3]



data_UDF_stack_z7 = np.loadtxt('Muv_re/UDF12_zdropouts_stack.data', comments='#')
x_UDF_stack_z7=data_UDF_stack_z7[:,0]
y_UDF_stack_z7=data_UDF_stack_z7[:,2]
raderr_UDF_stack_z7 = data_UDF_stack_z7[:,1]
magerr_UDF_stack_z7 = data_UDF_stack_z7[:,3]




data_UDF_bright_z7 = np.loadtxt('UDF12_zdropouts_bright.data', comments='#')
x_UDF_bright_z7=data_UDF_bright_z7[:,0]
y_UDF_bright_z7=data_UDF_bright_z7[:,2]
raderr_UDF_bright_z7 = data_UDF_bright_z7[:,1]
magerr_UDF_bright_z7 = data_UDF_bright_z7[:,3]



data_UDF_bright_z8 = np.loadtxt('UDF12_ydropouts_bright.data', comments='#')
x_UDF_bright_z8=data_UDF_bright_z8[:,0]
y_UDF_bright_z8=data_UDF_bright_z8[:,2]
raderr_UDF_bright_z8 = data_UDF_bright_z8[:,1]
magerr_UDF_bright_z8 = data_UDF_bright_z8[:,3]

data_UDF_faint_z7 = np.loadtxt('UDF12_zdropouts_faint.data', comments='#')
x_UDF_faint_z7=data_UDF_faint_z7[:,0]
y_UDF_faint_z7=data_UDF_faint_z7[:,2]
raderr_UDF_faint_z7 = data_UDF_faint_z7[:,1]
magerr_UDF_faint_z7 = data_UDF_faint_z7[:,3]



data_UDF_faint_z8 = np.loadtxt('UDF12_ydropouts_faint.data', comments='#')
x_UDF_faint_z8=data_UDF_faint_z8[:,0]
y_UDF_faint_z8=data_UDF_faint_z8[:,2]
raderr_UDF_faint_z8 = data_UDF_faint_z8[:,1]
magerr_UDF_faint_z8 = data_UDF_faint_z8[:,3]


data_UDF_bin1_z7 = np.loadtxt('UDF12_zdropouts_bin1.data', comments='#')
x_UDF_bin1_z7=data_UDF_bin1_z7[:,0]
y_UDF_bin1_z7=data_UDF_bin1_z7[:,2]
raderr_UDF_bin1_z7 = data_UDF_bin1_z7[:,1]
magerr_UDF_bin1_z7 = data_UDF_bin1_z7[:,3]

#data_UDF_bin2_z7 = np.loadtxt('UDF12_zdropouts_bin2.data', comments='#')
#x_UDF_bin2_z7=data_UDF_bin2_z7[:,0]
#y_UDF_bin2_z7=data_UDF_bin2_z7[:,2]
#raderr_UDF_bin2_z7 = data_UDF_bin2_z7[:,1]
#magerr_UDF_bin2_z7 = data_UDF_bin2_z7[:,3]

data_UDF_bin3_z7 = np.loadtxt('UDF12_zdropouts_bin3.data', comments='#')
x_UDF_bin3_z7=data_UDF_bin3_z7[:,0]
y_UDF_bin3_z7=data_UDF_bin3_z7[:,2]
raderr_UDF_bin3_z7 = data_UDF_bin3_z7[:,1]
magerr_UDF_bin3_z7 = data_UDF_bin3_z7[:,3]

#data_UDF_bin4_z7 = np.loadtxt('UDF12_zdropouts_bin4.data', comments='#')
#x_UDF_bin4_z7=data_UDF_bin4_z7[:,0]
#y_UDF_bin4_z7=data_UDF_bin4_z7[:,2]
#raderr_UDF_bin4_z7 = data_UDF_bin4_z7[:,1]
#magerr_UDF_bin4_z7 = data_UDF_bin4_z7[:,3]



data7_bin1 = np.loadtxt('catalogue_z7-bin1.cat', comments='#')
mag7_bin1=data7_bin1[:,74]
rad7_bin1=data7_bin1[:,73]

data7_bin2 = np.loadtxt('catalogue_z7-bin2.cat', comments='#')
mag7_bin2=data7_bin2[:,74]
rad7_bin2=data7_bin2[:,73]

data7_bin3 = np.loadtxt('catalogue_z7-bin3.cat', comments='#')
mag7_bin3=data7_bin3[:,74]
rad7_bin3=data7_bin3[:,73]

data7_bin4 = np.loadtxt('catalogue_z7-bin4.cat', comments='#')
mag7_bin4=data7_bin4[:,74]
rad7_bin4=data7_bin4[:,73]





ax.errorbar(mag8,rad8, yerr=[yerr_l8,yerr_u8],xerr=[xerr_l8,xerr_u8],fmt='.',color="r",label="This Work",elinewidth=1,alpha = .6)
plt.scatter(mag_not_multiple8,size_not_multiple8,s=400,color="r",marker=',',alpha = .6)

plt.scatter(mag_multiple8,size_multiple8,s=1200,color="r",marker='*',alpha = .6)
plt.scatter(-mag_multiple8,size_multiple8,s=100,color="k",marker='*',label="Multiple Component",alpha = .6)

ax.errorbar(x_UDF,y_UDF, yerr=raderr_UDF,xerr=magerr_UDF,fmt='.',color="k",label="Ono+13",elinewidth=1,alpha = .6)
plt.scatter(x_UDF,y_UDF,s=100,color="k",marker=',',alpha = .6)

ax.errorbar(x_UDF_b1,y_UDF_b1, yerr=raderr_UDF_b1,xerr=magerr_UDF_b1,fmt='.',color="k",elinewidth=1,alpha = .6)
#plt.scatter(x_UDF_b1,y_UDF_b1,s=400,color="k",marker='.',alpha = .6)

#ax.errorbar(x_UDF_b1_ave,y_UDF_b1_ave, yerr=raderr_UDF_b1_ave,xerr=magerr_UDF_b1_ave,fmt='o',color="g",label="Ono+13 average",elinewidth=2,alpha = .6)
#plt.scatter(x_UDF_b1_ave,y_UDF_b1_ave,s=400,color="g",marker=',',alpha = .6)

ax.errorbar(x_UDF_stack,y_UDF_stack, yerr=raderr_UDF_stack,xerr=magerr_UDF_stack,fmt='.',color="k",label="Ono+13 stacked",elinewidth=1,alpha = .6)
plt.scatter(x_UDF_stack,y_UDF_stack,s=400,color="k",marker='.',alpha = .6)
plt.scatter(x_UDF_stack,y_UDF_stack,s=300,color="w",marker='.')


plt.scatter(mag_multiple8_UDF,size_multiple8_UDF,s=300,color="k",marker='*',alpha = .6)



data_Halo_z8 = np.loadtxt('mag-r_halo8.dat', comments='#')
x_halo_z8=data_Halo_z8[:,0]
y_halo_z8=data_Halo_z8[:,3]
#plt.plot(x_halo_z8,y_halo_z8/20,color="#A9A9A9",linewidth=0.5)
#plt.plot(x_halo_z8,y_halo_z8/25,color="#A9A9A9",linewidth=0.5)
#plt.plot(x_halo_z8,y_halo_z8/30,color="#A9A9A9",linewidth=0.5)
#plt.plot(x_halo_z8,y_halo_z8/35,color="#A9A9A9",linewidth=0.5)
#plt.plot(x_halo_z8,y_halo_z8/40,color="#A9A9A9",linewidth=0.5)
#plt.plot(x_halo_z8,y_halo_z8/50,color="#A9A9A9",linewidth=0.5)
#plt.plot(x_halo_z8,y_halo_z8/60,color="#A9A9A9",linewidth=0.5)


rad7_bright_combine = np.r_[rad7_b1,y_UDF_bright_z7]
rad7_faint_combine = np.r_[rad7_b2,y_UDF_faint_z7]

rad8_bright_combine = np.r_[rad8_b1,y_UDF_bright_z8]
rad8_faint_combine = np.r_[y_UDF_faint_z8]

mag7_bright_combine = np.r_[mag7_b1,x_UDF_bright_z7]
mag7_faint_combine = np.r_[mag7_b2,x_UDF_faint_z7]

mag8_bright_combine = np.r_[mag8_b1,x_UDF_bright_z8]
mag8_faint_combine = np.r_[x_UDF_faint_z8]


mag7_bin1_combine = np.r_[mag7_bin1,x_UDF_bin1_z7]
rad7_bin1_combine = np.r_[rad7_bin1,y_UDF_bin1_z7]

mag7_bin2_combine = np.r_[mag7_bin2]
rad7_bin2_combine = np.r_[rad7_bin2]

mag7_bin3_combine = np.r_[mag7_bin3,x_UDF_bin3_z7]
rad7_bin3_combine = np.r_[rad7_bin3,y_UDF_bin3_z7]

mag7_bin4_combine = np.r_[mag7_bin4]
rad7_bin4_combine = np.r_[rad7_bin4]







ax.errorbar(np.average(mag8_bright_combine),np.average(rad8_bright_combine), yerr=math.sqrt(np.var(rad8_bright_combine)),xerr=math.sqrt(np.var(mag8_bright_combine)),fmt='.',color="b",label="This Work + Ono+13 (average)",elinewidth=2,alpha = 1.0)
#plt.scatter(np.average(mag8_bright_combine),np.average(rad8_bright_combine),s=800,color="b",marker='o',alpha = .3)

ax.errorbar(np.average(mag8_faint_combine),np.average(rad8_faint_combine), yerr=math.sqrt(np.var(rad8_faint_combine)),xerr=math.sqrt(np.var(mag8_faint_combine)),fmt='.',color="b",elinewidth=2,alpha = 1.0)
#plt.scatter(np.average(mag8_faint_combine),np.average(rad8_faint_combine),s=800,color="b",marker='o',alpha = .3)


plt.scatter(-18.3,1.01,s=200,color="k",marker=',',alpha = .6)
plt.scatter(-18.3,0.94,s=200,color="k",marker='*',alpha = .6)
plt.scatter(-18.3,0.87,s=200,color="k",marker='o',alpha = .6)
plt.scatter(-18.3,0.87,s=150,color="w",marker='o',alpha = 1.0)
plt.scatter(-18.3,0.80,s=300,color="k",marker='o',alpha = .3)

plt.text(-18,0.99, 'Bright', ha = 'center', va = 'bottom', size=10)
plt.text(-18,0.92, 'Multiple', ha = 'center', va = 'bottom', size=10)
plt.text(-18,0.85, 'Stack', ha = 'center', va = 'bottom', size=10)
plt.text(-18,0.78, 'Average', ha = 'center', va = 'bottom', size=10)

plt.text(-18,1.13, 'This Work', ha = 'center', va = 'bottom', size=12,color='r')
plt.text(-18,1.06, 'Ono_13', ha = 'center', va = 'bottom', size=12,color='c')






plt.xticks([-21,-20,-19,-18])
plt.yticks([0.0, 0.4, 0.8,1.2])
plt.xlim(-21.5,-17.5)
plt.ylim(0,1.25)
plt.xlabel("$M_{UV}$",fontname='serif',fontsize=20)
plt.ylabel("$R_e$ [kpc]",fontname='serif',fontsize=20)
#plt.legend(loc='upper right')

filename = "magrad8.jpg"
plt.savefig(filename, dpi=300)

#############################################################################

# Size - Luminostity _ z7


pylab.clf()
fig = plt.figure()
ax = fig.add_subplot(111)


ax.errorbar(mag7,rad7, yerr=[yerr_l7,yerr_u7],xerr=[xerr_l7,xerr_u7],fmt='.',color="r",label="This Work",elinewidth=1,alpha = .6)
plt.scatter(mag7,rad7,s=200,color="r",marker=',',alpha = .6)




ax.errorbar(x_UDF_z7,y_UDF_z7, yerr=raderr_UDF_z7,xerr=magerr_UDF_z7,fmt=',',color="c",label="Ono+13",elinewidth=1,alpha = .6)
plt.scatter(x_UDF_z7,y_UDF_z7,s=200,color="c",marker=',',alpha = .6)


ax.errorbar(x_UDF_b1_z7,y_UDF_b1_z7, yerr=raderr_UDF_b1_z7,xerr=magerr_UDF_b1_z7,fmt=',',color="c",elinewidth=1,alpha = .6)
plt.scatter(x_UDF_b1_z7,y_UDF_b1_z7,s=10,color="c",marker=',',alpha = .6)

#ax.errorbar(x_UDF_b1_ave_z7,y_UDF_b1_ave_z7, yerr=raderr_UDF_b1_ave_z7,xerr=magerr_UDF_b1_ave_z7,fmt='o',color="g",label="Ono+13 average",elinewidth=2)
#plt.scatter(x_UDF_b1_ave_z7,y_UDF_b1_ave_z7,s=400,color="g",marker=',')


ax.errorbar(x_UDF_stack_z7,y_UDF_stack_z7, yerr=raderr_UDF_stack_z7,xerr=magerr_UDF_stack_z7,fmt='.',color="c",label="Ono+13 stacked",elinewidth=2,alpha = .6)
plt.scatter(x_UDF_stack_z7,y_UDF_stack_z7,s=200,color="c",marker='o',alpha = .6)
plt.scatter(x_UDF_stack_z7,y_UDF_stack_z7,s=150,color="w",marker='o')


plt.scatter(mag_multiple7_UDF,size_multiple7_UDF,s=400,color="c",marker='*',alpha = .6)
plt.scatter(-mag_multiple8,size_multiple8,s=100,color="k",marker='*',label="Multiple Component",alpha = .6)



#ax.errorbar(np.average(mag7_bright_combine),np.average(rad7_bright_combine), yerr=math.sqrt(np.var(rad7_bright_combine)),xerr=math.sqrt(np.var(mag7_bright_combine)),fmt='.',color="k",label="This Work + Ono+13 (average)",elinewidth=2,alpha = .6)
#plt.scatter(np.average(mag7_bright_combine),np.average(rad7_bright_combine),s=600,color="k",marker='o',alpha = .3)
#ax.errorbar(np.average(mag7_faint_combine),np.average(rad7_faint_combine), yerr=math.sqrt(np.var(rad7_faint_combine)),xerr=math.sqrt(np.var(mag7_faint_combine)),fmt='.',color="k",elinewidth=2,alpha = .6)
#plt.scatter(np.average(mag7_faint_combine),np.average(rad7_faint_combine),s=600,color="k",marker='o',alpha = .3)

ax.errorbar(np.average(mag7_bin1_combine),np.average(rad7_bin1_combine), yerr=math.sqrt(np.var(rad7_bin1_combine)),xerr=math.sqrt(np.var(mag7_bin1_combine)),fmt='.',color="k",label="This Work + Ono+13 (average)",elinewidth=2,alpha = .6)
plt.scatter(np.average(mag7_bin1_combine),np.average(rad7_bin1_combine),s=800,color="k",marker='o',alpha = .3)

ax.errorbar(np.average(mag7_bin2_combine),np.average(rad7_bin2_combine), yerr=math.sqrt(np.var(rad7_bin2_combine)),xerr=math.sqrt(np.var(mag7_bin2_combine)),fmt='.',color="k",elinewidth=2,alpha = .6)
plt.scatter(np.average(mag7_bin2_combine),np.average(rad7_bin2_combine),s=800,color="k",marker='o',alpha = .3)

ax.errorbar(np.average(mag7_bin3_combine),np.average(rad7_bin3_combine), yerr=math.sqrt(np.var(rad7_bin3_combine)),xerr=math.sqrt(np.var(mag7_bin3_combine)),fmt='.',color="k",elinewidth=2,alpha = .6)
plt.scatter(np.average(mag7_bin3_combine),np.average(rad7_bin3_combine),s=800,color="k",marker='o',alpha = .3)

ax.errorbar(np.average(mag7_bin4_combine),np.average(rad7_bin4_combine), yerr=math.sqrt(np.var(rad7_bin4_combine)),xerr=math.sqrt(np.var(mag7_bin4_combine)),fmt='.',color="k",elinewidth=2,alpha = .6)
plt.scatter(np.average(mag7_bin4_combine),np.average(rad7_bin4_combine),s=800,color="k",marker='o',alpha = .3)




data_Halo_z7 = np.loadtxt('mag-r_halo7.dat', comments='#')
x_halo_z7=data_Halo_z7[:,0]
y_halo_z7=data_Halo_z7[:,3]
#plt.plot(x_halo_z7,y_halo_z7/20,color="#A9A9A9",linewidth=0.5)
#plt.plot(x_halo_z7,y_halo_z7/25,color="#A9A9A9",linewidth=0.5)
#plt.plot(x_halo_z7,y_halo_z7/30,color="#A9A9A9",linewidth=0.5)
#plt.plot(x_halo_z7,y_halo_z7/35,color="#A9A9A9",linewidth=0.5)
#plt.plot(x_halo_z7,y_halo_z7/40,color="#A9A9A9",linewidth=0.5)
#plt.plot(x_halo_z7,y_halo_z7/50,color="#A9A9A9",linewidth=0.5)
#plt.plot(x_halo_z7,y_halo_z7/60,color="#A9A9A9",linewidth=0.5)

plt.scatter(-18.3,1.01,s=200,color="k",marker=',',alpha = .6)
plt.scatter(-18.3,0.94,s=200,color="k",marker='*',alpha = .6)
plt.scatter(-18.3,0.87,s=200,color="k",marker='o',alpha = .6)
plt.scatter(-18.3,0.87,s=150,color="w",marker='o',alpha = 1.0)
plt.scatter(-18.3,0.80,s=300,color="k",marker='o',alpha = .3)

plt.text(-18,0.99, 'Bright', ha = 'center', va = 'bottom', size=10)
plt.text(-18,0.92, 'Multiple', ha = 'center', va = 'bottom', size=10)
plt.text(-18,0.85, 'Stack', ha = 'center', va = 'bottom', size=10)
plt.text(-18,0.78, 'Average', ha = 'center', va = 'bottom', size=10)

plt.text(-18,1.13, 'This Work', ha = 'center', va = 'bottom', size=12,color='r')
plt.text(-18,1.06, 'Ono_13', ha = 'center', va = 'bottom', size=12,color='c')




plt.xticks([-21,-20,-19,-18])
plt.yticks([0.0, 0.4, 0.8,1.2])
plt.xlim(-21.5,-17.5)
plt.ylim(0,1.25)
plt.xlabel("$M_{UV}$",fontname='serif',fontsize=20)
plt.ylabel("$R_e$ [kpc]",fontname='serif',fontsize=20)
#plt.legend(loc='upper right')

filename = "magrad7.jpg"
plt.savefig(filename, dpi=300)





#############################################################################

# Halo Radius - Stellar Radius z ~ 7

pylab.clf()
fig = plt.figure()
ax = fig.add_subplot(111)


plt.scatter(mag7,rad7,c='r',label = "This Work",s=200)
plt.scatter(x_UDF_z7,y_UDF_z7,c="c", label = "Ono+13",s=200)
plt.scatter(x_UDF_b1_z7,y_UDF_b1_z7,c="c",s=200)
plt.scatter(mag_multiple7_UDF,size_multiple7_UDF,s=800,color="c",marker='*')
plt.scatter(-mag_multiple7_UDF,size_multiple7_UDF,s=200,color="k",marker='*',label = "Multiple Component")


#plt.plot(x_halo_z7,y_halo_z7/20,color="#A9A9A9",linewidth=0.5,label = "R_halo / 20")
#plt.plot(x_halo_z7,y_halo_z7/25,color="#A9A9A9",linewidth=0.5,label = "R_halo / 25")
#plt.plot(x_halo_z7,y_halo_z7/30,color="#A9A9A9",linewidth=0.5,label = "R_halo / 30")
#plt.plot(x_halo_z7,y_halo_z7/35,color="#A9A9A9",linewidth=0.5,label = "R_halo / 35")
#plt.plot(x_halo_z7,y_halo_z7/40,color="#A9A9A9",linewidth=0.5,label = "R_halo / 40")
#plt.plot(x_halo_z7,y_halo_z7/50,color="#A9A9A9",linewidth=0.5,label = "R_halo / 50")
#plt.plot(x_halo_z7,y_halo_z7/60,color="#A9A9A9",linewidth=0.5,label = "R_halo / 60")
#plt.plot(x_halo_z7,y_halo_z7/100,color="#A9A9A9",linewidth=0.5,label = "R_halo / 100")

plt.plot(x_halo_z7,y_halo_z7/20,color="k",linewidth=0.5,label = "R_halo / 20")
#plt.plot(x_halo_z7,y_halo_z7/25,color="#A9A9A9",linewidth=0.5,label = "R_halo / 25")
#plt.plot(x_halo_z7,y_halo_z7/30,color="#A9A9A9",linewidth=0.5,label = "R_halo / 30")
#plt.plot(x_halo_z7,y_halo_z7/35,color="#A9A9A9",linewidth=0.5,label = "R_halo / 35")
plt.plot(x_halo_z7,y_halo_z7/40,color="k",linewidth=0.5,label = "R_halo / 40")
#plt.plot(x_halo_z7,y_halo_z7/50,color="#A9A9A9",linewidth=0.5,label = "R_halo / 50")
#plt.plot(x_halo_z7,y_halo_z7/60,color="#A9A9A9",linewidth=0.5,label = "R_halo / 60")
plt.plot(x_halo_z7,y_halo_z7/100,color="k",linewidth=0.5,label = "R_halo / 100")

plt.xlim(-21.5,-17.5)
plt.ylim(0,1.2)
plt.xticks([-21,-20,-19,-18])
plt.yticks([0.0, 0.4, 0.8,1.2])
plt.legend(loc='upper right')
plt.xlabel("$M_{UV}$",fontname='serif',fontsize=20)
plt.ylabel("$R_e$ [kpc]",fontname='serif',fontsize=20)


filename = "r_halo-r_stellar_z7.jpg"
plt.savefig(filename, dpi=300)



#############################################################################

# Halo Radius - Stellar Radius z ~ 8

pylab.clf()
fig = plt.figure()
ax = fig.add_subplot(111)


plt.scatter(mag8,rad8,c='r', label = "This Work",s=200)
plt.scatter(x_UDF,y_UDF,c="c", label = "Ono+13",s=200)
plt.scatter(x_UDF_b1,y_UDF_b1,c="c",s=200)
plt.scatter(mag_multiple8,size_multiple8,c='r', s=800,marker='*')
plt.scatter(mag_multiple8_UDF,size_multiple8_UDF,s=800,color="c",marker='*')
plt.scatter(-mag_multiple7_UDF,size_multiple7_UDF,s=200,color="k",marker='*',label = "Multiple Component")


plt.plot(x_halo_z8,y_halo_z8/20,color="k",linewidth=0.5,label = "R_halo / 20")
#plt.plot(x_halo_z8,y_halo_z8/25,color="#A9A9A9",linewidth=0.5,label = "R_halo / 25")
#plt.plot(x_halo_z8,y_halo_z8/30,color="#A9A9A9",linewidth=0.5,label = "R_halo / 30")
#plt.plot(x_halo_z8,y_halo_z8/35,color="#A9A9A9",linewidth=0.5,label = "R_halo / 35")
plt.plot(x_halo_z8,y_halo_z8/40,color="k",linewidth=0.5,label = "R_halo / 40")
#plt.plot(x_halo_z8,y_halo_z8/50,color="#A9A9A9",linewidth=0.5,label = "R_halo / 50")
#plt.plot(x_halo_z8,y_halo_z8/60,color="#A9A9A9",linewidth=0.5,label = "R_halo / 60")
plt.plot(x_halo_z8,y_halo_z8/100,color="k",linewidth=0.5,label = "R_halo / 100")


plt.xlim(-21.5,-17.5)
plt.ylim(0,1.2)
plt.xticks([-21,-20,-19,-18])
plt.yticks([0.0, 0.4, 0.8,1.2])
plt.legend(loc='upper right')
plt.xlabel("$M_{UV}$",fontname='serif',fontsize=20)
plt.ylabel("$R_e$ [kpc]",fontname='serif',fontsize=20)


filename = "r_halo-r_stellar_z8.jpg"
plt.savefig(filename, dpi=300)




#############################################################################

# Radius - Redshift


pylab.clf()
fig = plt.figure()
ax = fig.add_subplot(111)

#data_average_rad6 = np.loadtxt('average_rad_z6.dat', comments='#')
#redshift_average_rad6=data_average_rad6[1-1]
#rad6_average_rad6=data_average_rad6[2-1]
#yerr_l_average_rad6=data_average_rad6[3-1]
#yerr_u_average_rad6=data_average_rad6[4-1]
#ax.errorbar(redshift_average_rad6,rad6_average_rad6, yerr=[yerr_l_average_rad6,yerr_u_average_rad6],fmt='o',color="r",label="This Work")

data_average_rad7 = np.loadtxt('average_rad_z7.dat', comments='#')
redshift_average_rad7=data_average_rad7[1-1]
rad7_average_rad7=data_average_rad7[2-1]
yerr_l_average_rad7=data_average_rad7[3-1]
yerr_u_average_rad7=data_average_rad7[4-1]

#ax.errorbar(redshift_average_rad7,rad7_average_rad7, yerr=[[yerr_l_average_rad7],[yerr_l_average_rad7]],fmt='o',color="r",label="This Work",elinewidth=2)


data_average_rad8 = np.loadtxt('average_rad_z8.dat', comments='#')
redshift_average_rad8=data_average_rad8[1-1]
rad8_average_rad8=data_average_rad8[2-1]
yerr_l_average_rad8=data_average_rad8[3-1]
yerr_u_average_rad8=data_average_rad8[4-1]

#ax.errorbar(redshift_average_rad8,rad8_average_rad8, yerr=[[yerr_l_average_rad8],[yerr_u_average_rad8]],fmt='o',color="r",elinewidth=2)


data_average_rad8_UDFb1ave = np.loadtxt('redshift_re/UDF12_ydropouts_b1_ave.data', comments='#')
redshift_average_rad8_UDFb1ave=data_average_rad8_UDFb1ave[1-1]
rad8_average_rad8_UDFb1ave=data_average_rad8_UDFb1ave[2-1]
yerr_l_average_rad8_UDFb1ave=data_average_rad8_UDFb1ave[3-1]
yerr_u_average_rad8_UDFb1ave=data_average_rad8_UDFb1ave[3-1]


data_average_rad7_UDFb1ave = np.loadtxt('redshift_re/UDF12_zdropouts_b1_ave.data', comments='#')
redshift_average_rad7_UDFb1ave=data_average_rad7_UDFb1ave[1-1]
rad7_average_rad7_UDFb1ave=data_average_rad7_UDFb1ave[2-1]
yerr_l_average_rad7_UDFb1ave=data_average_rad7_UDFb1ave[3-1]
yerr_u_average_rad7_UDFb1ave=data_average_rad7_UDFb1ave[3-1]

data_average_rad_UDF_Oeschgalfit = np.loadtxt('redshift_re/Oesch2010_Fig3_bright_galfit.data', comments='#')
redshift_average_rad_UDF_Oeschgalfit=data_average_rad_UDF_Oeschgalfit[:,1-1]
rad_average_rad_UDF_Oeschgalfit=data_average_rad_UDF_Oeschgalfit[:,2-1]
yerr_l_average_rad_UDF_Oeschgalfit=data_average_rad_UDF_Oeschgalfit[:,3-1]
yerr_u_average_rad_UDF_Oeschgalfit=data_average_rad_UDF_Oeschgalfit[:,3-1]

data_average_rad_UDF_Oeschgalfit_lowz = np.loadtxt('redshift_re/Oesch2010_Fig3_bright_lowz.data', comments='#')
redshift_average_rad_UDF_Oeschgalfit_lowz=data_average_rad_UDF_Oeschgalfit_lowz[:,1-1]
rad_average_rad_UDF_Oeschgalfit_lowz=data_average_rad_UDF_Oeschgalfit_lowz[:,2-1]
yerr_l_average_rad_UDF_Oeschgalfit_lowz=data_average_rad_UDF_Oeschgalfit_lowz[:,3-1]
yerr_u_average_rad_UDF_Oeschgalfit_lowz=data_average_rad_UDF_Oeschgalfit_lowz[:,3-1]

ax.errorbar(redshift_average_rad8_UDFb1ave,rad8_average_rad8_UDFb1ave, yerr=[[yerr_l_average_rad8_UDFb1ave],[yerr_u_average_rad8_UDFb1ave]],fmt='.',color="c",label="Ono+13",elinewidth=4,alpha = .6)
plt.scatter(redshift_average_rad8_UDFb1ave,rad8_average_rad8_UDFb1ave,s=100,color="c",marker='o',alpha = .6)

ax.errorbar(redshift_average_rad_UDF_Oeschgalfit+0.1,rad_average_rad_UDF_Oeschgalfit, yerr=[yerr_l_average_rad_UDF_Oeschgalfit,yerr_u_average_rad_UDF_Oeschgalfit],fmt='.',color="g",label="Oesch+10",elinewidth=4,alpha = .6)
plt.scatter(redshift_average_rad_UDF_Oeschgalfit+0.1,rad_average_rad_UDF_Oeschgalfit,s=100,color="g",marker='o',alpha = .6)


ax.errorbar(redshift_average_rad7_UDFb1ave,rad7_average_rad7_UDFb1ave, yerr=[[yerr_l_average_rad7_UDFb1ave],[yerr_u_average_rad7_UDFb1ave]],fmt='.',color="c",elinewidth=4,alpha = .6)
plt.scatter(redshift_average_rad7_UDFb1ave,rad7_average_rad7_UDFb1ave,s=100,color="c",marker='o',alpha = .6)

ax.errorbar(redshift_average_rad_UDF_Oeschgalfit_lowz,rad_average_rad_UDF_Oeschgalfit_lowz, yerr=[yerr_l_average_rad_UDF_Oeschgalfit_lowz,yerr_u_average_rad_UDF_Oeschgalfit_lowz],fmt='.',color="k",label="Bouwens+04",elinewidth=4,alpha = .6)
plt.scatter(redshift_average_rad_UDF_Oeschgalfit_lowz,rad_average_rad_UDF_Oeschgalfit_lowz,s=100,color="k",marker='o',alpha = .6)


plt.scatter(7.0,np.average(rad7_bright_combine),s=400,color="r",marker='o',alpha = .6)
ax.errorbar(7.0,np.average(rad7_bright_combine), yerr=np.std(rad7_bright_combine),fmt='.',color="r",label="This Work + Ono+13",elinewidth=4,alpha = .6)

ax.errorbar(8.0,np.average(rad8_bright_combine), yerr=np.std(rad8_bright_combine),fmt='.',color="r",elinewidth=4,alpha = .6)
plt.scatter(8.0,np.average(rad8_bright_combine),s=400,color="r",marker='o',alpha = .6)



plt.xlim(1.8,9.0)
plt.ylim(0.1,3)
plt.yscale("log")
plt.xlabel("Redshift",fontname='serif',fontsize=20)
plt.ylabel("$R_e$ [kpc]",fontname='serif',fontsize=20)
plt.legend(loc='upper right')
plt.minorticks_on()

filename = "radius-redshift.jpg"
plt.savefig(filename, dpi=300)



#############################################################################

# Sigma_SFR - Redshift


pylab.clf()
fig = plt.figure()
ax = fig.add_subplot(111)

data_average_sigmasfr7 = np.loadtxt('/Users/KawamataRyota/Dropbox/Reports@Dropbox/Gravitational_Lensing@Dropbox/glafic/Abell2744/SFR_Abell/SFR-7_combine_bright.dat', comments='#')
redshift_average_sigmasfr7=data_average_sigmasfr7[:,1-1]+0.1
sigmasfr_average_sfr7=data_average_sigmasfr7[:,10-1]
sigmasfr_average_sigmasfr7=data_average_sigmasfr7[:,13-1]
ax.errorbar(np.average(redshift_average_sigmasfr7),np.average(sigmasfr_average_sigmasfr7), yerr=np.std(sigmasfr_average_sigmasfr7),fmt='.',color="r",label="This Work + Ono+13",elinewidth=4,alpha = .6)
plt.scatter(np.average(redshift_average_sigmasfr7),np.average(sigmasfr_average_sigmasfr7),s=200,color="r",marker='o',alpha = .6)



data_average_sigmasfr8 = np.loadtxt('/Users/KawamataRyota/Dropbox/Reports@Dropbox/Gravitational_Lensing@Dropbox/glafic/Abell2744/SFR_Abell/SFR-8_combine_bright.dat', comments='#')
redshift_average_sigmasfr8=data_average_sigmasfr8[:,1-1]+0.1
sigmasfr_average_sfr8=data_average_sigmasfr8[:,10-1]
sigmasfr_average_sigmasfr8=data_average_sigmasfr8[:,13-1]
ax.errorbar(np.average(redshift_average_sigmasfr8),np.average(sigmasfr_average_sigmasfr8), yerr=np.std(sigmasfr_average_sigmasfr8),fmt='.',color="r",elinewidth=4,alpha = .6)
plt.scatter(np.average(redshift_average_sigmasfr8),np.average(sigmasfr_average_sigmasfr8),s=200,color="r",marker='o',alpha = .6)


data_average_sigmasfr7_UDFb1ave = np.loadtxt('/Users/KawamataRyota/Dropbox/Reports@Dropbox/Gravitational_Lensing@Dropbox/glafic/Abell2744/Plot_Abell/redshift_SigmaSFR/UDF12_zdropouts_b1_ave.data', comments='#')
redshift_average_sigmasfr7_UDFb1ave=data_average_sigmasfr7_UDFb1ave[1-1]
sigmasfr_average_sigmasfr7_UDFb1ave=data_average_sigmasfr7_UDFb1ave[2-1]
yerr_l_average_sigmasfr7_UDFb1ave=data_average_sigmasfr7_UDFb1ave[3-1]
yerr_u_average_sigmasfr7_UDFb1ave=data_average_sigmasfr7_UDFb1ave[3-1]
ax.errorbar(redshift_average_sigmasfr7_UDFb1ave,sigmasfr_average_sigmasfr7_UDFb1ave, yerr=[[yerr_l_average_sigmasfr7_UDFb1ave],[yerr_u_average_sigmasfr7_UDFb1ave]],fmt='.',color="c",label="Ono+13",elinewidth=4,alpha = .6)
plt.scatter(redshift_average_sigmasfr7_UDFb1ave,sigmasfr_average_sigmasfr7_UDFb1ave,s=200,color="c",marker='o',alpha = .6)


data_average_sigmasfr8_UDFb1ave = np.loadtxt('/Users/KawamataRyota/Dropbox/Reports@Dropbox/Gravitational_Lensing@Dropbox/glafic/Abell2744/Plot_Abell/redshift_SigmaSFR/UDF12_ydropouts_b1_ave.data', comments='#')
redshift_average_sigmasfr8_UDFb1ave=data_average_sigmasfr8_UDFb1ave[1-1]
sigmasfr_average_sigmasfr8_UDFb1ave=data_average_sigmasfr8_UDFb1ave[2-1]
yerr_l_average_sigmasfr8_UDFb1ave=data_average_sigmasfr8_UDFb1ave[3-1]
yerr_u_average_sigmasfr8_UDFb1ave=data_average_sigmasfr8_UDFb1ave[3-1]
ax.errorbar(redshift_average_sigmasfr8_UDFb1ave,sigmasfr_average_sigmasfr8_UDFb1ave, yerr=[[yerr_l_average_sigmasfr8_UDFb1ave],[yerr_u_average_sigmasfr8_UDFb1ave]],fmt='.',color="c",elinewidth=4,alpha = .6)
plt.scatter(redshift_average_sigmasfr8_UDFb1ave,sigmasfr_average_sigmasfr8_UDFb1ave,s=200,color="c",marker='o',alpha = .6)




data_average_sigmasfr_Oesch = np.loadtxt('/Users/KawamataRyota/Dropbox/Reports@Dropbox/Gravitational_Lensing@Dropbox/glafic/Abell2744/Plot_Abell/redshift_SigmaSFR/Oesch2010_Fig5.data', comments='#')
redshift_average_sigmasfr_Oesch=data_average_sigmasfr_Oesch[:,1-1]
sigmasfr_average_sigmasfr_Oesch=data_average_sigmasfr_Oesch[:,2-1]

yerr_l_average_sigmasfr_Oesch=data_average_sigmasfr_Oesch[:,2-1]-data_average_sigmasfr_Oesch[:,4-1]
yerr_u_average_sigmasfr_Oesch=data_average_sigmasfr_Oesch[:,3-1]-data_average_sigmasfr_Oesch[:,2-1]


ax.errorbar(redshift_average_sigmasfr_Oesch,sigmasfr_average_sigmasfr_Oesch, yerr=[yerr_l_average_sigmasfr_Oesch,yerr_u_average_sigmasfr_Oesch],fmt='.',color="g",label="Oesch+10",elinewidth=4,alpha = .6)
plt.scatter(redshift_average_sigmasfr_Oesch,sigmasfr_average_sigmasfr_Oesch,s=200,color="g",marker='o',alpha = .6)



plt.xlim(0.001,10.0)
plt.ylim(0.01,40)
ax.set_yscale('log')
plt.xlabel("Redshift",fontname='serif',fontsize=20)
plt.ylabel("${\Sigma_{SFR}}$",fontname='serif',fontsize=20)
plt.legend(loc='upper left')

filename = "sigmasfr-redshift.jpg"
plt.savefig(filename, dpi=300)



###############################################################################

# Size - Luminosity _ Abell w/beta

pylab.clf()

plt.scatter(mag8,rad8,c=beta8,vmin=-4.0,vmax=0.0, label = "This Work")
plt.scatter(mag7,rad7,c=beta7,vmin=-4.0,vmax=0.0)
plt.scatter(mag_multiple8,size_multiple8,c=beta_multiple8,vmin=-4.0,vmax=0.0, label = "This Work",s=400,marker='*')


plt.xlim(-21.5,-17.5)
plt.ylim(0,1.5)
plt.legend(loc='upper right')
plt.colorbar()


filename = "magsize_beta.jpg"
plt.savefig(filename, dpi=300)

###############################################################################

# Beta - Size _ Abell

pylab.clf()


plt.scatter(beta8,rad8,c=beta8,vmin=-4.0,vmax=0.0, label = "This Work")
plt.scatter(beta7,rad7,c=beta7,vmin=-4.0,vmax=0.0)
plt.scatter(beta_multiple8,size_multiple8,c=beta_multiple8,vmin=-4.0,vmax=0.0, label = "This Work",s=400,marker='*')


#plt.xlim(-21.5,-18.5)
#plt.ylim(0,1.5)
plt.legend(loc='upper left')
plt.colorbar()
plt.xlabel(r"$\beta$",fontname='serif',fontsize=20)
plt.ylabel("$R_e$ [kpc]",fontname='serif',fontsize=20)

filename = "betasize_beta.jpg"
plt.savefig(filename, dpi=300)


###############################################################################

# Beta - Size _ UDF

pylab.clf()

data_UDF_betasize7 = np.loadtxt('beta-size_UDF_z7.dat', comments='#')
data_UDF_betasize8 = np.loadtxt('beta-size_UDF_z8.dat', comments='#')

beta_UDF_betasize7 = -2.0 + 4.39*(data_UDF_betasize7[:,6-1] - data_UDF_betasize7[:,8-1])
beta_UDF_betasize8 = -2.0 + 8.98*(data_UDF_betasize8[:,7-1] - data_UDF_betasize8[:,8-1])
a1600_UDF_betasize7 = 4.43 + beta_UDF_betasize7*1.99
a1600_UDF_betasize8 = 4.43 + beta_UDF_betasize8*1.99
size_UDF_betasize7 = data_UDF_betasize7[:,3-1]
size_UDF_betasize8 = data_UDF_betasize8[:,3-1]
mag_UDF_betasize7 = data_UDF_betasize7[:,2-1]
mag_UDF_betasize8 = data_UDF_betasize8[:,2-1]
sigmasfr_UDF_betasize7 = 2.7*10**11/data_UDF_betasize7[:,3-1]**2*10**(-(data_UDF_betasize7[:,2-1]+48.6)/2.5)
sigmasfr_UDF_betasize8 = 2.7*10**11/data_UDF_betasize8[:,3-1]**2*10**(-(data_UDF_betasize8[:,2-1]+48.6)/2.5)
sigmasfr_corr_UDF_betasize7 = 2.7*10**11/data_UDF_betasize7[:,3-1]**2*10**(-(data_UDF_betasize7[:,2-1]+48.6-a1600_UDF_betasize7)/2.5)
sigmasfr_corr_UDF_betasize8 = 2.7*10**11/data_UDF_betasize8[:,3-1]**2*10**(-(data_UDF_betasize8[:,2-1]+48.6-a1600_UDF_betasize8)/2.5)

plt.scatter(beta_UDF_betasize8,size_UDF_betasize8,c=beta_UDF_betasize8,vmin=-4.0,vmax=0.0, label = "This Work")
plt.scatter(beta_UDF_betasize7,size_UDF_betasize7,c=beta_UDF_betasize7,vmin=-4.0,vmax=0.0, label = "This Work")

#plt.xlim(-21.5,-18.5)
#plt.ylim(0,1.5)
#plt.legend(loc='upper right')
plt.colorbar()


filename = "betasize_beta_UDF.jpg"
plt.savefig(filename, dpi=300)


###############################################################################

# Beta - Size _ Combine

pylab.clf()


#plt.scatter(beta8,rad8,c='r',vmin=-4.0,vmax=0.0, label = "This Work",s=200)
#plt.scatter(beta7,rad7,c='r',vmin=-4.0,vmax=0.0,s=200)
#plt.scatter(beta_UDF_betasize8,size_UDF_betasize8,c='c',vmin=-4.0,vmax=0.0, label = "Ono+13",marker='o',s=200)
#plt.scatter(beta_UDF_betasize7,size_UDF_betasize7,c='c',vmin=-4.0,vmax=0.0,marker='o',s=200)


plt.scatter(beta8,rad8,c=mag8,vmin=-22,vmax=-18, label = "This Work",s=200)
plt.scatter(beta7,rad7,c=mag7,vmin=-22,vmax=-18,s=200)
plt.scatter(beta_UDF_betasize8,size_UDF_betasize8,c=mag_UDF_betasize8,vmin=-22,vmax=-18, label = "Ono+13",marker='o',s=200)
plt.scatter(beta_UDF_betasize7,size_UDF_betasize7,c=mag_UDF_betasize7,vmin=-22,vmax=-18,marker='o',s=200)


plt.colorbar()

#plt.xlim(-21.5,-18.5)
#plt.ylim(0,1.5)
plt.legend(loc='upper left')

plt.xlabel(r"$\beta$",fontname='serif',fontsize=20)
plt.ylabel("$R_e$ [kpc]",fontname='serif',fontsize=20)


filename = "betasize_beta_combine.jpg"
plt.savefig(filename, dpi=300)


corr = np.corrcoef(np.r_[beta8,beta7,beta_UDF_betasize8,beta_UDF_betasize7], np.r_[rad8,rad7,size_UDF_betasize8,size_UDF_betasize7])
print(corr[0,1])



###############################################################################

# Beta - Mag _ Combine

pylab.clf()


plt.scatter(beta8,mag8,c='r',vmin=-4.0,vmax=0.0, label = "This Work",s=200)
plt.scatter(beta7,mag7,c='r',vmin=-4.0,vmax=0.0,s=200)
#plt.scatter(beta_multiple8,mag_multiple8,c='r',vmin=-4.0,vmax=0.0,s=800,marker='*')

plt.scatter(beta_UDF_betasize8,mag_UDF_betasize8,c='c',vmin=-4.0,vmax=0.0, label = "Ono+13",marker='o',s=200)
plt.scatter(beta_UDF_betasize7,mag_UDF_betasize7,c='c',vmin=-4.0,vmax=0.0,marker='o',s=200)



#plt.xlim(-21.5,-18.5)
#plt.ylim(0,1.5)
plt.legend(loc='upper right')

plt.xlabel(r"$\beta$",fontname='serif',fontsize=20)
plt.ylabel("$M_{UV}$",fontname='serif',fontsize=20)


filename = "betamag_beta_combine.jpg"
plt.savefig(filename, dpi=300)
###############################################################################

# Beta - Sigmasfr _ Combine

pylab.clf()


#plt.scatter(beta8,sigmasfr8,c='r',vmin=0,vmax=50.0, label = "This Work",s=200)
#plt.scatter(beta7,sigmasfr7,c='r',vmin=0,vmax=50.0,s=200)
#plt.scatter(beta_multiple8,mag_multiple8,c='r',vmin=-4.0,vmax=0.0,s=800,marker='*')
#plt.scatter(beta_UDF_betasize8,sigmasfr_UDF_betasize8,c='c',vmin=0,vmax=50.0, label = "Ono+13",marker='o',s=200)
#plt.scatter(beta_UDF_betasize7,sigmasfr_UDF_betasize7,c='c',vmin=0,vmax=50.0,marker='o',s=200)

#plt.scatter(beta8,sigmasfr8,c=rad8,vmin=0,vmax=1.0, label = "This Work",s=200)
#plt.scatter(beta7,sigmasfr7,c=rad7,vmin=0,vmax=1.0,s=200)
#plt.scatter(beta_multiple8,mag_multiple8,c='r',vmin=-4.0,vmax=0.0,s=800,marker='*')
#plt.scatter(beta_UDF_betasize8,sigmasfr_UDF_betasize8,c=size_UDF_betasize8,vmin=0,vmax=1.0, label = "Ono+13",marker='o',s=200)
#plt.scatter(beta_UDF_betasize7,sigmasfr_UDF_betasize7,c=size_UDF_betasize7,vmin=0,vmax=1.0,marker='o',s=200)

plt.scatter(beta8,sigmasfr8,c=mag8,vmin=-21,vmax=-18, label = "This Work",s=200,cmap = 'jet_r'
)
plt.scatter(beta7,sigmasfr7,c=mag7,vmin=-21,vmax=-18,s=200,cmap = 'jet_r'
)
#plt.scatter(beta_multiple8,mag_multiple8,c='r',vmin=-4.0,vmax=0.0,s=800,marker='*')
plt.scatter(beta_UDF_betasize8,sigmasfr_UDF_betasize8,c=mag_UDF_betasize8,vmin=-21,vmax=-18, label = "Ono+13",marker='o',s=200,cmap = 'jet_r'
)
plt.scatter(beta_UDF_betasize7,sigmasfr_UDF_betasize7,c=mag_UDF_betasize7,vmin=-21,vmax=-18,marker='o',s=200,cmap = 'jet_r'
)


#ax.set_yscale('log')
plt.yscale('log')

#plt.xlim(-5.0,0.0)
plt.ylim(0.8,4000)
plt.legend(loc='upper right')
plt.colorbar()

plt.xlabel(r"$\beta$",fontname='serif',fontsize=20)
plt.ylabel("${\Sigma_{SFR}}$",fontname='serif',fontsize=20)


filename = "betasigmasfr_beta_combine.jpg"
plt.savefig(filename, dpi=300)


pylab.clf()


plt.scatter(beta8,sigmasfr8_corr,c=mag8,vmin=-21,vmax=-18, label = "This Work",s=200,cmap = 'jet_r')
plt.scatter(beta7,sigmasfr7_corr,c=mag7,vmin=-21,vmax=-18,s=200,cmap = 'jet_r')
plt.scatter(beta_UDF_betasize8,sigmasfr_corr_UDF_betasize8,c=mag_UDF_betasize8,vmin=-21,vmax=-18, label = "Ono+13",marker='o',s=200,cmap = 'jet_r')
plt.scatter(beta_UDF_betasize7,sigmasfr_corr_UDF_betasize7,c=mag_UDF_betasize7,vmin=-21,vmax=-18,marker='o',s=200,cmap = 'jet_r')


plt.yscale('log')

#plt.xlim(-5.0,0.0)
plt.ylim(0.8,4000)
plt.legend(loc='upper right')
plt.colorbar()

plt.xlabel(r"$\beta$",fontname='serif',fontsize=20)
plt.ylabel("${\Sigma_{SFR}}$",fontname='serif',fontsize=20)


filename = "betasigmasfr_beta_combine_corrected.jpg"
plt.savefig(filename, dpi=300)



##################################

pylab.clf()


plt.scatter(beta8,sigmasfr8_corr,c='r',vmin=-21,vmax=-18, label = "This Work",s=200,cmap = 'jet_r')
plt.scatter(beta7,sigmasfr7_corr,c='r',vmin=-21,vmax=-18,s=200,cmap = 'jet_r')
plt.scatter(beta_UDF_betasize8,sigmasfr_corr_UDF_betasize8,c='r',vmin=-21,vmax=-18, label = "Ono+13",marker='o',s=200,cmap = 'jet_r')
plt.scatter(beta_UDF_betasize7,sigmasfr_corr_UDF_betasize7,c='r',vmin=-21,vmax=-18,marker='o',s=200,cmap = 'jet_r')


plt.scatter(beta8,sigmasfr8,c='c',vmin=-21,vmax=-18, label = "This Work",s=200,cmap = 'jet_r')
plt.scatter(beta7,sigmasfr7,c='c',vmin=-21,vmax=-18,s=200,cmap = 'jet_r')
plt.scatter(beta_UDF_betasize8,sigmasfr_UDF_betasize8,c='c',vmin=-21,vmax=-18, label = "Ono+13",marker='o',s=200,cmap = 'jet_r')
plt.scatter(beta_UDF_betasize7,sigmasfr_UDF_betasize7,c='c',vmin=-21,vmax=-18,marker='o',s=200,cmap = 'jet_r')


plt.yscale('log')

#plt.xlim(-5.0,0.0)
plt.ylim(0.8,4000)
plt.legend(loc='upper right')

plt.xlabel(r"$\beta$",fontname='serif',fontsize=20)
plt.ylabel("${\Sigma_{SFR}}$",fontname='serif',fontsize=20)


filename = "check.jpg"
plt.savefig(filename, dpi=300)


###############################################################################

# Mag - Sigmasfr _ Combine

pylab.clf()


#plt.scatter(beta8,sigmasfr8,c='r',vmin=0,vmax=50.0, label = "This Work",s=200)
#plt.scatter(beta7,sigmasfr7,c='r',vmin=0,vmax=50.0,s=200)
#plt.scatter(beta_multiple8,mag_multiple8,c='r',vmin=-4.0,vmax=0.0,s=800,marker='*')
#plt.scatter(beta_UDF_betasize8,sigmasfr_UDF_betasize8,c='c',vmin=0,vmax=50.0, label = "Ono+13",marker='o',s=200)
#plt.scatter(beta_UDF_betasize7,sigmasfr_UDF_betasize7,c='c',vmin=0,vmax=50.0,marker='o',s=200)

#plt.scatter(beta8,sigmasfr8,c=rad8,vmin=0,vmax=1.0, label = "This Work",s=200)
#plt.scatter(beta7,sigmasfr7,c=rad7,vmin=0,vmax=1.0,s=200)
#plt.scatter(beta_multiple8,mag_multiple8,c='r',vmin=-4.0,vmax=0.0,s=800,marker='*')
#plt.scatter(beta_UDF_betasize8,sigmasfr_UDF_betasize8,c=size_UDF_betasize8,vmin=0,vmax=1.0, label = "Ono+13",marker='o',s=200)
#plt.scatter(beta_UDF_betasize7,sigmasfr_UDF_betasize7,c=size_UDF_betasize7,vmin=0,vmax=1.0,marker='o',s=200)

plt.scatter(mag8,sigmasfr8,c="r",vmin=-21,vmax=-18, label = "This Work",s=200,cmap = 'jet_r')
plt.scatter(mag7,sigmasfr7,c="r",vmin=-21,vmax=-18,s=200,cmap = 'jet_r')
#plt.scatter(beta_multiple8,mag_multiple8,c='r',vmin=-4.0,vmax=0.0,s=800,marker='*')
plt.scatter(mag_UDF_betasize8,sigmasfr_UDF_betasize8,c="c",vmin=-21,vmax=-18, label = "Ono+13",marker='o',s=200,cmap = 'jet_r')
plt.scatter(mag_UDF_betasize7,sigmasfr_UDF_betasize7,c="c",vmin=-21,vmax=-18,marker='o',s=200,cmap = 'jet_r')


#ax.set_yscale('log')
plt.yscale('log')

#plt.xlim(-5.0,0.0)
plt.ylim(0.8,4000)
plt.legend(loc='upper right')

plt.xlabel("$M_{UV}$",fontname='serif',fontsize=20)
plt.ylabel("${\Sigma_{SFR}}$",fontname='serif',fontsize=20)


filename = "magsigmasfr_beta_combine.jpg"
plt.savefig(filename, dpi=300)


###############################################################################

# Beta - Mag _ Abell

pylab.clf()


plt.scatter(beta8,mag8,c=beta8,vmin=-4.0,vmax=0.0, label = "This Work")
plt.scatter(beta7,mag7,c=beta7,vmin=-4.0,vmax=0.0)
plt.scatter(beta_multiple8,mag_multiple8,c=beta_multiple8,vmin=-4.0,vmax=0.0, label = "This Work",s=400,marker='*')


#plt.xlim(-21.5,-18.5)
#plt.ylim(0,1.5)
#plt.legend(loc='upper right')
plt.colorbar()

plt.xlabel(r"$\beta$",fontname='serif',fontsize=20)
plt.ylabel("$M_{UV}$",fontname='serif',fontsize=20)

filename = "betamag_beta.jpg"
plt.savefig(filename, dpi=300)

corr = np.corrcoef(np.r_[beta8,beta7,beta_UDF_betasize8,beta_UDF_betasize7], np.r_[mag8,mag7,mag_UDF_betasize8,mag_UDF_betasize7])
print(corr[0,1])


###############################################################################
pylab.clf()

data_UDF_betasize7_previous = np.loadtxt('beta-size_UDF_z7_previous.dat', comments='#')
data_UDF_betasize8_previous = np.loadtxt('beta-size_UDF_z8_previous.dat', comments='#')

beta_UDF_betasize7_previous = -2.0 + 4.39*(data_UDF_betasize7_previous[:,3-1])
beta_UDF_betasize8_previous = -2.0 + 4.39*(data_UDF_betasize8_previous[:,3-1])
size_UDF_betasize7_previous = data_UDF_betasize7_previous[:,2-1]
size_UDF_betasize8_previous = data_UDF_betasize8_previous[:,2-1]

plt.scatter(beta_UDF_betasize8,size_UDF_betasize8,c='r',vmin=-4.0,vmax=0.0, label = "This Work")
plt.scatter(beta_UDF_betasize7,size_UDF_betasize7,c='r',vmin=-4.0,vmax=0.0)
plt.scatter(beta_UDF_betasize8_previous,size_UDF_betasize8_previous,c='c',vmin=-4.0,vmax=0.0, label = "This Work")
plt.scatter(beta_UDF_betasize7_previous,size_UDF_betasize7_previous,c='c',vmin=-4.0,vmax=0.0)

#plt.xlim(-21.5,-18.5)
#plt.ylim(0,1.5)
#plt.legend(loc='upper right')


filename = "betasize_compare.jpg"
plt.savefig(filename, dpi=300)












