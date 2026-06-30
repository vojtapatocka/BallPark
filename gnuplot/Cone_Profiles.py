import numpy as np
import matplotlib.pyplot as plt
from pylab import figure, show, legend, xlabel
#from scipy.optimize import curve_fit
import os.path
import sys

fit = True
volfit = True
L1norm = True
volumes = True
ref = 2
cutsummit = False
left = True
dirp = "../../"
dirr = "../run/"
dirf = "../figs/"
# Manual reading of the cone center from the respective Profiles.csv
icenl = [1210, 1090, 820] #[1264, 1100, 940]
#icenl = [0, 0, 0]

stuff = sys.argv[1].split(',')
for item in stuff:
    exec(item)

ixmin = 0
if fit: volumes=True

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = [20, 9]

fig1 = figure()
ax = fig1.add_subplot(111)

thirtyx = np.linspace(0,500,500)
thirtyy = -(thirtyx-500.)*np.tan(np.deg2rad(30.))
ax.plot(thirtyx[:],thirtyy[:],c='black')

sfile = dirp+"ProfileNECalorisA.csv"
print("Reading "+sfile)
xdatfull = np.loadtxt(sfile, usecols=(1), skiprows=1, delimiter=',')
zdatfull = np.loadtxt(sfile, usecols=(2), skiprows=1, delimiter=',')
icen = icenl[0]
print("Center point chosen at the distance: ", xdatfull[icen]/1.e3,' km')
if left:
    xdat1 = np.flip(-xdatfull[80:icen+1] + xdatfull[icen])
    zdat1 = np.flip(zdatfull[80:icen+1] - zdatfull[-1])
else:
    xdat1 = xdatfull[icen:] - xdatfull[icen]
    zdat1 = zdatfull[icen:] - zdatfull[-1]
if cutsummit:
    if left:
        icut = icen - 1000
        xdat1 = np.flip(-xdatfull[80:(icen-icut)] + xdatfull[icen])
        zdat1 = np.flip(zdatfull[80:(icen-icut)] - zdatfull[-1])        
    else:
        icut = 1470 - icen
        xdat1 = xdatfull[(icen+icut):] - xdatfull[icen]
        zdat1 = zdatfull[(icen+icut):] - zdatfull[-1]    
if(ref==1): 
    ztopref = zdat1[0]
    xmaxref = xdat1[-1]

if volumes:
    voldat = 0.
    for ix in range(len(xdat1)-1):
        voldat += 0.5*np.pi*(xdat1[ix]+xdat1[ix+1])*(zdat1[ix]+zdat1[ix+1])*(xdat1[ix+1]-xdat1[ix])
    print("Volume of the profile [km3]: ", voldat/1.e9)
    if(ref==1): volref = voldat

ax.plot(xdat1,zdat1,label='NE Caloris A')
#ax.axvline(x=xdatfull[1264],c='C0',ymin=0.,ymax=2000.)
#ax.axvline(x=xdatfull[1470],c='C0',linestyle='--',ymin=0.,ymax=2000.)
#ax.axvline(x=xdatfull[928],c='C0',linestyle=':',ymin=0.,ymax=2000.)

sfile = dirp+"ProfileNWCaloris.csv"
print("Reading "+sfile)
xdatfull = np.loadtxt(sfile, usecols=(1), skiprows=1, delimiter=',')
zdatfull = np.loadtxt(sfile, usecols=(2), skiprows=1, delimiter=',')
icen = icenl[1]
xdat2 = xdatfull[icen:] - xdatfull[icen]
zdat2 = zdatfull[icen:] - zdatfull[-1]
if left:
    xdat2 = np.flip(-xdatfull[300:icen+1] + xdatfull[icen])
    zdat2 = np.flip(zdatfull[300:icen+1] - zdatfull[300] + 10.)
if(ref==2): 
    ztopref = zdat2[0]
    xmaxref = xdat2[-1]

if volumes:
    voldat = 0.
    for ix in range(len(xdat2)-1):
        voldat += 0.5*np.pi*(xdat2[ix]+xdat2[ix+1])*(zdat2[ix]+zdat2[ix+1])*(xdat2[ix+1]-xdat2[ix])
    print("Volume of the profile [km3]: ", voldat/1.e9)  
    if(ref==2): volref = voldat

ax.plot(xdat2,zdat2,label='NW Caloris')
#ax.axvline(x=xdatfull[1100],c='C1',ymin=0.,ymax=2000.)

sfile = dirp+"ProfileHeaney.csv"
print("Reading "+sfile)
xdatfull = np.loadtxt(sfile, usecols=(1), skiprows=1, delimiter=',')
zdatfull = np.loadtxt(sfile, usecols=(2), skiprows=1, delimiter=',')
icen = icenl[2]
xdat3 = xdatfull[icen:] - xdatfull[icen]
zdat3 = zdatfull[icen:] - zdatfull[-1]
if left:
    xdat3 = np.flip(-xdatfull[:icen+1] + xdatfull[icen])
    zdat3 = np.flip(zdatfull[:icen+1] - 0.5*(zdatfull[0]+zdatfull[-1]))
if(ref==3): 
    ztopref = zdat3[0]
    xmaxref = xdat3[-1]

if volumes:
    voldat = 0.
    for ix in range(len(xdat3)-1):
        voldat += 0.5*np.pi*(xdat3[ix]+xdat3[ix+1])*(zdat3[ix]+zdat3[ix+1])*(xdat3[ix+1]-xdat3[ix])
    print("Volume of the profile [km3]: ", voldat/1.e9)  
    if(ref==3): volref = voldat

ax.plot(xdat3,zdat3,label='Heaney')
#ax.axvline(x=xdatfull[940],ymin=0.,ymax=2000.,c='C2')

sfile = dirr+"cone.dat"
print("Reading "+sfile)
xcor = np.loadtxt(sfile, usecols=(0))
zcor = np.loadtxt(sfile, usecols=(1))
if cutsummit and ref==1:
    for ix in range(len(xcor)):
        if(xcor[ix]>xdat1[0]): 
            ixmin=ix
            break

if volumes:
    voldat = 0.
    for ix in range(ixmin,len(xcor)-1):
        if (xcor[ix]<=xmaxref):
            voldat += 0.5*np.pi*(xcor[ix]+xcor[ix+1])*(zcor[ix]+zcor[ix+1])*(xcor[ix+1]-xcor[ix])
        else:
            break
    print("Volume of the simulated cone [km3]: ", voldat/1.e9)  

volscl = volref/voldat
hscl = ztopref/zcor[ixmin]
error = 0.
error2 = 0.
if fit:
    if (ref==1):
        fitx = xdat1[1:]
        fity = zdat1[1:]
    if (ref==2):
        fitx = xdat2[1:]
        fity = zdat2[1:]
    if (ref==3):
        fitx = xdat3[1:]
        fity = zdat3[1:]
    simy = np.interp(fitx[:], xcor[:], zcor[:]*volscl)
    simy2 = np.interp(fitx[:], xcor[:], zcor[:]*hscl)
    #ax.plot(fitx,fity,'1',label='fity')
    #ax.plot(fitx,simy,'2',ms='30',label='simy')
    #ax.plot(fitx,simy2,'3',ms='30',label='simy2')
    error = 0.
    error2 = 0.
    sqdat = 0.
    dx = (fitx[2]-fitx[1])
    power = 2
    if L1norm: power = 1
    for ix in range(len(simy)):
        faktor = 1.
        if volfit: faktor = fitx[ix]        
        error += np.power(np.abs(simy[ix]-fity[ix]),power)*dx*faktor
        error2 += np.power(np.abs(simy2[ix]-fity[ix]),power)*dx*faktor
        sqdat += np.power(np.abs(fity[ix]),power)*dx*faktor
    error /= sqdat
    error2 /= sqdat
    print("Volume_normalized_error ", error)
    print("Central_height_normalized_error ", error2)

if volumes: ax.plot(xcor,zcor*volscl,'o',label=("sim (volume) %.1f %%" % (100*error)))
ax.plot(xcor,zcor*hscl,'.',label=("sim (height) %.1f %%" % (100*error2)))
ax.axhline(y=0.,xmin=0.,xmax=15000.,c='black',linestyle='--')

ax.set_xlabel('distance [m]')
ax.set_ylabel('height [m]')
ax.set_ylim(bottom=-200., top=1200.)
ax.set_xlim(right=20000., left=0.)

plt.legend()
plt.savefig(dirf+"Cone.png", bbox_inches='tight')

