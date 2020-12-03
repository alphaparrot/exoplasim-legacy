import numpy as np
import matplotlib.pyplot as plt
import sys

cc = 299792458.0 #Speed of light

def readspec(sfile,cgs=False):
    f=open(sfile,"r")
    ftxt = f.read()
    if ftxt[0]=='#':
       newstyle=True
    else:
       newstyle=False
    f.close()
    if newstyle:
       lines = ftxt.split('\n')[:-2] #chop off timing line
       nhl = 0
       for entry in lines:
          if entry[0]=='#':
             nhl+=1
          else:
             break
       units = lines[nhl-2].split()[-1]
       if lines[nhl-2].split()[1]=="(ANGSTROM)":
          cgs=False
       nhi = 0
       for entry in lines[nhl:]:
          if float(entry.split()[0])<0.1: #angstroms
             nhi+=1
          else:
             break
       print "Cutting at %d: "%(nhl+nhi),lines[nhl+nhi]
       lines = lines[nhl+nhi:]
       wvls = np.zeros(len(lines))
       fluxes = np.zeros(len(lines))
       n=0
       for l in lines:
          wvls[n] = float(l.split()[0])
          fluxes[n] = float(l.split()[1])
          n+=1
          #print l.split()
    else:      
       lines = ftxt.split('\n')[:-1]
       wvls = np.zeros(len(lines))
       fluxes = np.zeros(len(lines))
       n=0
       for l in lines:
           wvls[n] = float(l.split()[0])
           fluxes[n] = 10**(float(l.split()[1].replace('D','E'))-8.0)
           n+=1
       units='ergs/sec/cm^2/angstrom'
    if not cgs:
        units='W/m^2/micron'
        wvls[:] = wvls[:]*1.0e-4 #microns
        fluxes[:] = fluxes[:]*10.0 #W/m^2/micron
    return wvls,fluxes,units

def coarsen(wvls,fluxes,w1,w2,num=5000):
    inds = np.arange(len(wvls))
    nlow = (inds[wvls>w1])[0]
    nhigh= (inds[wvls>w2])[0]
    
    inc = int(nhigh-nlow)/num
    
    waves = np.zeros(num)
    energ = wvls*fluxes
    flxs2 = np.zeros(num)
    
    for n in range(0,num):
        v0 = nlow+n*inc-1
        v1 = nlow+n*inc
        intg = 0.0
        for v in range(0,inc):
            intg+=0.5*(wvls[v1+v]-wvls[v0+v])*(fluxes[v1+v]+fluxes[v0+v])
        dnu = wvls[v1+inc]-wvls[v0]
        waves[n] = 0.5*(wvls[v1+inc]+wvls[v0])
        flxs2[n] = intg/dnu
        
    return waves,flxs2

def integrate_trap(x,y):
    net = 0
    for n in range(1,len(y)):
        net += 0.5*(x[n]-x[n-1])*(y[n]+y[n-1])
    return net

def normalize(wvls,fluxes):
    netf = np.trapz(fluxes,x=wvls)
    factor = 1366.941858/netf
    return fluxes*factor

def writedat(wvls,fluxes,name):
    f=open(name+".dat","w")
    sdat = ' Wavelength    Flux  \n'
    for n in range(0,len(wvls)):
        sdat+=str(wvls[n])+' '+str(fluxes[n])+'\n'
    f.write(sdat)
    f.close()
    
if __name__=="__main__":
    sfile = sys.argv[1]
    name = sys.argv[2]
    #wave0 = float(sys.argv[3]) #microns
    #wave1 = float(sys.argv[4]) #microns
    try:
        numw = int(sys.argv[3])
    except:
        numw = 2048
    try:
        norm = sys.argv[4]
    except:
        norm = False
    w,f,u = readspec(sfile)
    print u,w.min(),w.max()
    plt.plot(w,f)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    #wc,fc = coarsen(w,f,wave0,wave1,num=numw*4)
    #plt.plot(wc,fc)
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.xlabel("$\lambda$ [$\mu$m]")
    #plt.ylabel("$F_\lambda$ [W/m$^2$/$\mu$m]")
    #plt.show()
    w2 = np.concatenate([np.geomspace(0.2,0.75,num=numw/2+1)[:-1],np.geomspace(0.75,100.0,num=numw/2)])
    plt.plot(w2)
    plt.yscale('log')
    plt.show()
    f2 = np.interp(w2,w,f)
    E0 = np.trapz(f,x=w)
    E1 = np.trapz(f2,x=w2)
    print abs(E0-E1)/E0
    if norm:
        factor = 1366.941858/E1
        f2*=factor
        print np.trapz(f2,x=w2)
        print np.trapz(normalize(w,f),x=w)
    #print f2[np.argwhere(w2>40.0)],w[-10:]
    plt.plot(w2,f2)
    #plt.plot(wc,fc)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("$\lambda$ [$\mu$m]")
    plt.ylabel("$F_\lambda$ [W/m$^2$/$\mu$m]")
    plt.show()
    wvref = np.loadtxt("wvref.txt")
    f3 = np.interp(wvref,w2,f2)
    writedat(w2,f2/cc,name+"_hr")
    writedat(wvref,f3/cc,name)
    plt.plot(wvref,f3)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("$\lambda$ [$\mu$m]")
    plt.ylabel("$\lambda F_\lambda$ [W/m$^2$]")
    plt.show()
