import numpy as np
import matplotlib.pyplot as plt
from pylab import figure, show, legend, xlabel
#from scipy.optimize import curve_fit
import os.path
import sys

stuff = sys.argv[1].split(',')
for item in stuff:
    exec(item)

plt.rcParams.update({'font.size': 16})
plt.rcParams["figure.figsize"] = [20, 9]

fig1 = figure()
ax = fig1.add_subplot(111)

thirtyx = np.linspace(0,500,500)
thirtyy = -(thirtyx-500.)*np.tan(np.deg2rad(30.))
ax.plot(thirtyx[:],thirtyy[:],c='black')

sfile = "../run/cone.dat"
print("Reading "+sfile)
xcor = np.loadtxt(sfile, usecols=(0))
zcor = np.loadtxt(sfile, usecols=(1))

ax.plot(xcor,zcor,label='cone')
ax.plot(xcor,zcor,'x')

sfile = "../run/baselayer.dat"
zcor = np.loadtxt(sfile, usecols=(1))
ax.plot(xcor,zcor,'1',label='base l.')
sfile = "../run/toplayer.dat"
zcor = np.loadtxt(sfile, usecols=(1))
ax.plot(xcor,zcor,'2',label='top l.')


ax.set_xlabel('distance [m]')
ax.set_ylabel('height [m]')
ax.set_ylim(bottom=0., top=1000.)
ax.set_xlim(right=17000., left=0.)

plt.legend()
plt.savefig("../figs/Cone.png", bbox_inches='tight')

