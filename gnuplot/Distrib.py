import numpy as np
import matplotlib.pyplot as plt
from pylab import figure, show, legend, xlabel
from scipy.optimize import curve_fit
import os.path
import sys

dname='dcoll'
fitdat=True

stuff = sys.argv[1].split(',')
for item in stuff:
    exec(item)

plt.rcParams.update({'font.size': 16})
plt.rcParams["figure.figsize"] = [20, 9]

fig1 = figure()
ax = fig1.add_subplot(111)

sfile = "../run/"+dname+".dat"
print("Reading "+sfile)
xcor = np.loadtxt(sfile, usecols=(0))
zcor = np.loadtxt(sfile, usecols=(1))

ax.plot(xcor,zcor,'o',fillstyle='none',label='data')

if fitdat:
    def flognor(x, A, B, C):
        return  C*np.exp(-0.5*np.power(((np.log10(x)-np.log10(A))/B),2))
    def flognormal(x, A, B, C):
        return  C*np.exp(-0.5*np.power(((np.log(x)-np.log(A))/B),2))/x
    def fnormal(x, A, B):
        return  B*np.exp(-0.5*np.power(x/A,2))*np.sin(np.deg2rad(x))
        
if(dname=='dcoll'): 
    ax.set_xlabel('diameter [m]')
    if fitdat:
        guess = np.array([0.01, 0.5, max(zcor[:])])
        pars, covariance = curve_fit(flognor, xcor[:], zcor[:], guess[:], maxfev=100000)
        guess = np.array([0.01, 0.5, max(zcor[:])*0.01])
        pars2, covariance2 = curve_fit(flognormal, xcor[:], zcor[:], guess[:], maxfev=100000)
        print("Obtained parameters of the fit: ", pars[:])
        truth = flognor(xcor[:], 0.04, 0.3, max(zcor[:]))
if(dname=='vcoll'): 
    ax.set_xlabel('velocity [m/s]')
    if fitdat:
        guess = np.array([100., 0.5, max(zcor[:])])
        pars, covariance = curve_fit(flognor, xcor[:], zcor[:], guess[:], maxfev=10000)
        pars2, covariance2 = curve_fit(flognormal, xcor[:], zcor[:], guess[:], maxfev=10000)
        truth = flognor(xcor[:], 46., 0.2, max(zcor[:]))
if (dname=='dcoll' or dname=='vcoll') and fitdat:
    fitdata = flognor(xcor[:], pars[0], pars[1], pars[2])
    fitdata2 = flognormal(xcor[:], pars2[0], pars2[1], pars2[2])
    ax.plot(xcor,truth,'1',ms=7,label='reference')
    ax.plot(xcor,fitdata,'2',ms=10,label=(r"$\exp\left(-\frac{\left[\log_{10}(x)-\log_{10}(%.3f)\right]^2}{2\,(%.3f)^2}\right)$" % (pars[0], abs(pars[1]))))
    ax.plot(xcor,fitdata2,label=(r"$\exp\left(-\frac{\left[\log(x)-\log(%.3f)\right]^2}{2\,(%.3f)^2}\right)/x$" % (pars2[0], abs(pars2[1]))))

if(dname=='angcoll'): 
    ax.set_xlabel(r"angle $\alpha$ [deg]")
    if fitdat:
        guess = np.array([1., 1.e6])
        pars, covariance = curve_fit(fnormal, xcor[:], zcor[:], guess[:], maxfev=10000)
        fitdata = fnormal(xcor[:], pars[0], pars[1])
        truth = fnormal(xcor[:], 31., 1.)
        #ax.plot(xcor,truth*max(zcor[:]),'2',label='ref run')
        ax.plot(xcor,fitdata,label=(r"$\exp(-\frac{\alpha^2}{2\,(%.3f)^2})*\sin(\alpha)$" % (abs(pars[0]))))

ax.set_ylabel('count')

plt.legend()
plt.savefig("../figs/"+dname+".png", bbox_inches='tight')

