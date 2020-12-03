import numpy as np
import os
import math
import struct

#This version lets the model relax.

def isflat(weatheringfile):
  f=open(weatheringfile,'r')
  r=f.read()
  f.close()
  r=r.split('\n')[1:-1]
  if len(r) < 13+2:
    return False
  else:
    d=0
    while d<len(r):
      r[d]=r[d].split()
      d+=1
    temps=[]
    n=0
    while n<len(r):
      temps.append(float(r[n][1]))
      n+=1
    n=len(temps)-3
    tt=np.arange(13)+1
    linfits=[]
    while n<len(temps):
      sample=temps[n-(13-1):n+1]
      linfit=np.polyfit(tt,sample,1)[0]
      linfits.append(abs(linfit))
      n+=1
    avglinfit = (linfits[-3]+linfits[-2]+linfits[-1])/3.0
    if avglinfit <= 0.05:
      return np.mean(temps[-(13-2):])
    else:
      return False
  


if __name__=="__main__":
  wf=open("weathering.pso","w")
  wf.write("     CO2       AVG SURF T   WEATHERING    OUTGASSING      DpCO2       NEW CO2\n")
  wf.close()
  EXP="MOST"
  NCPU=8
  #os.system("rm -f plasim_restart") #Uncomment for a fresh run when you haven't cleaned up beforehand
  os.system("rm -f Abort_Message")
  year=0
  minyear=13
  relaxed=False
  while year < minyear or not relaxed:
    year+=1
    dataname=EXP+".%04d"%year
    diagname=EXP+"_DIAG.%04d"%year
    restname=EXP+"_REST.%03d"%year
    snowname=EXP+"_SNOW_%1d"%(year%5)
    os.system("mpiexec -np "+str(NCPU)+" most_plasim_t21_l10_p"+str(NCPU)+".x")
    os.system("[ -e restart_dsnow ] && rm restart_dsnow")
    os.system("[ -e restart_xsnow ] && rm restart_xsnow")
    os.system("[ -e Abort_Message ] && exit 1")
    os.system("[ -e plasim_output ] && mv plasim_output "+dataname)
    os.system("[ -e plasim_diag ] && mv plasim_diag "+diagname)
    os.system("[ -e plasim_status ] && cp plasim_status plasim_restart")
    os.system("[ -e plasim_status ] && mv plasim_status "+restname)
    os.system("[ -e restart_snow ] && mv restart_snow "+snowname)
    os.system("[ -e "+dataname+" ] && ./burn7.x -n <example.nl>burnout "+dataname+" "+dataname+".nc")
    os.system("[ -e "+dataname+" ] && cp "+dataname+" "+EXP+"_OUT.%04d"%n)
    os.system("[ -e "+dataname+".nc ] && rm "+dataname)
    relaxed=isflat("weathering.pso")
    
