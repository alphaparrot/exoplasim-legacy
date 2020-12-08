import numpy as np
import os
import sys
import netCDF4 as nc
import glob
import time

gplasim = True
TIMELIMIT = 1.44e5

#This version lets the model relax.

def spatialmath(lt,ln,variable,mean=True,radius=6.371e6):
    lt1 = np.zeros(len(lt)+1)
    lt1[0] = 90
    for n in range(0,len(lt)-1):
        lt1[n+1] = 0.5*(lt[n]+lt[n+1])
    lt1[-1] = -90
    ln1 = np.zeros(len(ln)+1)
    ln1[0] = -2.8125
    for n in range(0,len(ln)-1):
        ln1[n+1] = 0.5*(ln[n]+ln[n+1])
    ln1[-1] = 360.0-2.8125
    
    lt1*=np.pi/180.0
    ln1*=np.pi/180.0
    
    darea = np.zeros((len(lt),len(ln)))
    for jlat in range(0,len(lt)):
        for jlon in range(0,len(ln)):
            dln = ln1[jlon+1]-ln1[jlon]
            darea[jlat,jlon] = (np.sin(lt1[jlat])-np.sin(lt1[jlat+1]))*dln
    
    svar = variable*darea
    if mean:
        outvar = np.sum(svar)/np.sum(darea)
    else:
        outvar = np.sum(svar) * radius**2
    
    return outvar

def isflat(key="ts",mean=True,radius=6.371e6,baseline=13,threshhold=0.05):
    #Key is the netCDF variable to evaluate, mean toggles whether to track the average or total,
    #radius is the planet radius in meters, and baseline is the number of years over which to measure
    #slope. Default is to track surface temperature. Threshhold is the maximum slope we'll allow.
  files = sorted(glob.glob("*.nc"))
  nfiles = len(files)
  prior=False
  if len(glob.glob("thistory.ps*"))>0:
      thistory = np.loadtxt("thistory.pso")
      nfiles += len(thistory)
      prior=True
  dd = np.zeros(nfiles)
  nstart=0
  if prior:
      dd[:len(thistory)] = thistory[:]
      nstart=len(thistory)
  if len(files) < baseline+2:
    return False
  else:
    for n in range(0,len(files)):
        ncd = nc.Dataset(files[n],"r")
        variable = ncd.variables[key][:]
        if len(variable.shape)>3:
            variable = variable[:,-1,:,:]
        for m in range(0,variable.shape[0]):
            dd[n+nstart] += spatialmath(ncd.variables['lat'][:],ncd.variables['lon'][:],variable[m,:,:],
                                 mean=mean,radius=radius)
        dd[n+nstart] /= variable.shape[0] #Monthly mean
        ncd.close()
    n=len(dd)-3
    tt=np.arange(baseline)+1
    linfits=[]
    for n in range(len(dd)-3,len(dd)):
      sample=dd[n-(baseline-1):n+1]
      linfit=np.polyfit(tt,sample,1)[0]
      linfits.append(abs(linfit))
      
    avglinfit = (linfits[-3]+linfits[-2]+linfits[-1])/3.0
    if avglinfit <= 0.05:
      return np.mean(dd[-(baseline-2):])
    else:
      return False
  
def gethistory(key="ts",mean=True,radius=6.371e6):
    files = sorted(glob.glob("*.nc"))
    dd=np.zeros(len(files))
    for n in range(0,len(files)):
        ncd = nc.Dataset(files[n],"r")
        variable = ncd.variables[key][:]
        if len(variable.shape)>3:
            variable = variable[:,-1,:,:]
        for m in range(0,variable.shape[0]):
            dd[n] += spatialmath(ncd.variables['lat'][:],ncd.variables['lon'][:],variable[m,:,:],
                                 mean=mean,radius=radius)
        dd[n] /= variable.shape[0] #Monthly mean
        ncd.close()
    return dd
  
def hasnans():
    files = sorted(glob.glob("*.nc"))
    print("NetCDF  files:",files)
    if type(files)!=type([1,2,3]):
        files = [files,]
    ncd = nc.Dataset(files[-1],"r") #Should be most recent
    if np.sum(1.0*np.isnan(ncd.variables['ts'][-1,:]))>0.5:
        return True
    return False

def energybalanced(threshhold = 1.0e-4,baseline=50): #Takes an average of 200 years
    files = sorted(glob.glob("*.nc"))
    nfiles = len(files)
    prior=False
    if len(glob.glob("toahistory.ps*"))>0:
        toahistory = np.loadtxt("toahistory.pso")
        nfiles+=len(toahistory)
        shistory = np.loadtxt("shistory.pso")
        prior=True
    sbalance = np.zeros(nfiles)
    toabalance=np.zeros(nfiles)
    nstart=0
    if prior:
        sbalance[:len(toahistory)] = shistory[:]
        toabalance[:len(toahistory)] = toahistory[:]
        nstart = len(toahistory)
    if len(files) < baseline: #Run for minimum of baseline years
        return False
    else:
        for n in range(0,len(files)):
            ncd = nc.Dataset(files[n],"r")
            ntr = ncd.variables['ntr'][:]
            hfns = ncd.variables['hfns'][:]
            lat = ncd.variables['lat'][:]
            lon = ncd.variables['lon'][:]
            ncd.close()
            topt = np.zeros(12)
            bott = np.zeros(12)
            for m in range(0,12):
                topt[m] = spatialmath(lat,lon,ntr[m,:,:])
                bott[m] = spatialmath(lat,lon,hfns[m,:,:])
            sbalance[n+nstart] = np.mean(bott)
            toabalance[n+nstart] = np.mean(topt)
        savgs = []
        tavgs = []
        for n in range(9,len(sbalance)):
            savgs.append(abs(np.mean(sbalance[n-9:n+1]))) #10-year average energy balance
            tavgs.append(abs(np.mean(toabalance[n-9:n+1])))
        sslopes = []
        tslopes = []
        for n in range(4,len(savgs)): #5-baseline slopes in distance from energy balance
            sslopes.append(np.polyfit(np.arange(5)+1,savgs[n-4:n+1],1)[0])
            tslopes.append(np.polyfit(np.arange(5)+1,tavgs[n-4:n+1],1)[0])
        savgslope = abs(np.mean(sslopes[-30:])) #30-year average of 5-year slopes  
        tavgslope = abs(np.mean(tslopes[-30:]))
        os.system("echo '%02.8f  %02.8f'>>slopes.log"%(savgslope,tavgslope))
        print("%02.8f %02.8f"%(savgslope,tavgslope))
        if savgslope<threshhold and tavgslope<threshhold: #Both TOA and Surface are changing at average 
            return True                                  # of <0.1 mW/m^2/yr on 45-year baselines
        else:
            return False
        
def getbalance():
    files = sorted(glob.glob("*.nc"))
    ncd = nc.Dataset(files[-1],"r")
    ntr = ncd.variables['ntr'][:]
    hfns = ncd.variables['hfns'][:]
    lat = ncd.variables['lat'][:]
    lon = ncd.variables['lon'][:]
    ncd.close()
    topt = np.zeros(12)
    bott = np.zeros(12)
    for m in range(0,12):
        topt[m] = spatialmath(lat,lon,ntr[m,:,:])
        bott[m] = spatialmath(lat,lon,hfns[m,:,:])
    return (np.mean(bott),np.mean(topt))
    
        

if __name__=="__main__":
  if gplasim:
    wf=open("weathering.pso","w")
    wf.write("     CO2       AVG SURF T   WEATHERING    OUTGASSING      DpCO2       NEW CO2\n")
    wf.close()
  EXP="MOST"
  os.system("rm keepgoing")
  tstart = time.clock()
  NCPU=int(sys.argv[1])
  nlevs = int(sys.argv[2])
  #os.system("rm -f plasim_restart") #Uncomment for a fresh run when you haven't cleaned up beforehand
  os.system("rm -f Abort_Message")
  os.system("echo 'SURFACE      TOA'>balance.log")
  os.system("echo 'SURFACE      TOA'>slopes.log")
  exfiles = glob.glob("*DIAG*")
  year=len(exfiles)
  minyears=75
  maxyears=year+300
  relaxed=False
  while (year < minyears or not energybalanced(threshhold=4.0e-4)) and year<maxyears and (time.clock()-tstart)<=TIMELIMIT:
    year+=1
    dataname=EXP+".%04d"%year
    snapname=EXP+"_SNAP.%04d"%year
    diagname=EXP+"_DIAG.%04d"%year
    restname=EXP+"_REST.%03d"%year
    snowname=EXP+"_SNOW_%1d"%(year%5)
    os.system("mpiexec -np "+str(NCPU)+" most_plasim_t21_l"+str(nlevs)+"_p"+str(NCPU)+".x")
    os.system("[ -e restart_dsnow ] && rm restart_dsnow")
    os.system("[ -e restart_xsnow ] && rm restart_xsnow")
    os.system("[ -e Abort_Message ] && exit 1")
    os.system("[ -e plasim_output ] && mv plasim_output "+dataname)
    os.system("[ -e plasim_snapshot ] && mv plasim_snapshot "+snapname)
    os.system("[ -e plasim_diag ] && mv plasim_diag "+diagname)
    os.system("[ -e plasim_status ] && cp plasim_status plasim_restart")
    os.system("[ -e plasim_status ] && mv plasim_status "+restname)
    os.system("[ -e restart_snow ] && mv restart_snow "+snowname)
    os.system("[ -e "+dataname+" ] && ./burn7.x -n <example.nl>burnout "+dataname+" "+dataname+".nc")
    os.system("[ -e "+snapname+" ] && ./burn7.x -n <snapshot.nl>snapout "+snapname+" "+snapname+".nc")
    os.system("[ -e "+dataname+" ] && cp "+dataname+" "+EXP+"_OUT.%04d"%year)
    os.system("[ -e "+dataname+".nc ] && rm "+dataname)
    os.system("[ -e "+snapname+".nc ] && rm "+snapname)
    os.system("[ -e "+snapname+".nc ] && mv "+snapname+".nc snapshots/")
    if hasnans():
        os.system("echo 'NAN ENCOUNTERED'>>weathering.pso")
        break
    sb,tb = getbalance()
    os.system("echo '%02.6f  %02.6f'>>balance.log"%(sb,tb))
  os.system("rm keepgoing")
  if not hasnans() and not energybalanced(threshhold=4.0e-4):
    os.system("touch keepgoing")
    bott = gethistory(key="hfns")
    topt = gethistory(key="ntr")
    with open("shistory.pso","a+") as f:
        text='\n'+'\n'.join(bott.astype(str))
        f.write(text)
    with open("toahistory.pso","a+") as f:
        text='\n'+'\n'.join(topt.astype(str))
        f.write(text)
        
