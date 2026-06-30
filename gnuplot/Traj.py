import numpy as np
import matplotlib.pyplot as plt
from pylab import figure, show, legend, xlabel
#from scipy.optimize import curve_fit
import os.path
import sys

Nd=3
Nv=3
Na=3

stuff = sys.argv[1].split(',')
for item in stuff:
    exec(item)

print("Nd, Nv, Na ",Nd,Nv,Na)
plt.rcParams.update({'font.size': 16})
plt.rcParams["figure.figsize"] = [20, 9]
fig1 = figure()
ax = fig1.add_subplot(111)
#ax2 = fig1.add_subplot(122)

for id in range(1,Nd+1):
    for iv in range(1,Nv+1):
        for ia in range(1,Na+1):

            sfile = "../run/traj"+str(id)+"_"+str(iv)+"_"+str(ia)+".dat"
            print("Reading "+sfile)
            time = np.loadtxt(sfile, usecols=(0))
            xcor = np.loadtxt(sfile, usecols=(1))
            zcor = np.loadtxt(sfile, usecols=(2))
            vxcor = np.loadtxt(sfile, usecols=(3))
            vzcor = np.loadtxt(sfile, usecols=(4))
            dpar = np.loadtxt(sfile, usecols=(5), max_rows=1)
            vpar = np.loadtxt(sfile, usecols=(6), max_rows=1)
            apar = np.loadtxt(sfile, usecols=(7), max_rows=1)

            #ax.plot(xcor,zcor,label=("d= %.1f [cm], v= %.1f [m/s], a=%.1f [°], dur= %.1f [s]" % (dpar*100.,vpar,apar,time[-1])))
            ax.plot(xcor,zcor,label=("d= %.1f [cm], v= %.1f [m/s], a=%.1f [°]" % (dpar*100.,vpar,apar)))

ax.set_xlabel('distance [m]')
ax.set_ylabel('height [m]')
ax.set_ylim(bottom=0.)
ax.legend() #bbox_to_anchor=(0.6, 0.34), framealpha=1.0)

#plt.show()
plt.savefig("../figs/Traj.png", bbox_inches='tight')

