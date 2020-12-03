import numpy as np
import os
import math
import struct

#This version couples both a carbon cycle and a glacier model


def getmaxdsnow(filename1,filename2):
  f1=open(filename1,'rb')
  r1=f1.read()
  f1.close()
  
  f2=open(filename2,'rb')
  r2=f2.read()
  f2.close()
  
  dd1=np.array(struct.unpack('2048d',r1[4:-4])).reshape((32,64))
  dd2=np.array(struct.unpack('2048d',r2[4:-4])).reshape((32,64))
  
  dsnow = dd2-dd1
  maxdsnow = np.amax(np.abs(dsnow))
  
  return maxdsnow

def changeCO2(pCO2):
  nl=open("radmod_namelist","r")
  nltxt=nl.read().split('\n')
  nl.close()
  l=0
  while l<len(nltxt):
    line=nltxt[l].split()
    if line[0]=='CO2':
      line[2]=str(pCO2)
      nltxt[l]=' '+' '.join(line)
      break
    elif line[0]=='CO2=':
      line[1]=str(pCO2)
      nltxt[l]=' '+' '.join(line)
      break
    l+=1
  nltxt='\n'.join(nltxt)
  nl=open('radmod_namelist','w')
  nl.write(nltxt)
  nl.close()
  
def changep(psurf):
  nl=open("plasim_namelist","r")
  nltxt=nl.read().split('\n')
  nl.close()
  l=0
  while l<len(nltxt)-1:
    line=nltxt[l].split()
    if line[0]=='PSURF':
      line[2]=str(psurf)
      nltxt[l]=' '+' '.join(line)
      break
    elif line[0]=='PSURF=':
      line[1]=str(psurf)
      nltxt[l]=' '+' '.join(line)
      break
    l+=1
  nltxt='\n'.join(nltxt)
  nl=open('plasim_namelist','w')
  nl.write(nltxt)
  nl.close()  
  
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
  
def getct(weatheringfile):
  f=open(weatheringfile,'r')
  r=f.read().split('\n')[1:-1]
  f.close()
  l=r[-1].split()
  ct=float(l[1])
  return ct

def getndco2(weatheringfile):
  f=open(weatheringfile,'r')
  r=f.read().split('\n')[1:-1]
  f.close()
  l1=r[1].split()
  l2=r[-1].split()
  sco2=float(l1[0])*1e6
  tco2=float(l2[5])*1e6
  return (sco2,tco2)

def getdco2(weatheringfile):
  f=open(weatheringfile,'r')
  r=f.read().split('\n')[1:-1]
  f.close()
  l=r[-1].split()
  dco2=float(l[4])*1e6
  return dco2

if __name__=="__main__":
  start=True
  pCO2=330.0
  p0 = 1010670.0
  cyear=0
  n=0
  eCO2s=[]
  etemps=[]
  #changep(psurf*0.1)
  #co2s=getndco2('weathering.pso')
  #dco2=getdco2('weathering.pso')
  #ct=getct('weathering.pso')
  #eCO2s.append(co2s[1])
  #etemps.append(ct)
  #pCO2 = eCO2s[-1]
  psurf = p0 + pCO2
  changep(psurf*0.1)
  #os.system("rm -f plasim_restart")
  while n<300:
    changeCO2(pCO2/psurf*1e6) #Change pCO2 to CO2 ppmv
    wf=open("weathering.pso","w")
    wf.write("     CO2       AVG SURF T   WEATHERING    OUTGASSING      DpCO2       NEW CO2\n")
    wf.close()
    EXP="MOST"
    NCPU=16
    #os.system("rm -f plasim_restart") #Uncomment for a fresh run when you haven't cleaned up beforehand
    os.system("rm -f Abort_Message")
    year=0
    minyear=13
    relaxed=False
    while year < minyear or not relaxed:
      cyear+=1
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
    os.system("cp "+dataname+".nc "+EXP+"_OUT.%04d.nc"%n)
    os.system("cp "+restname+" "+EXP+"_REST.%04d"%n)
    os.system("cp "+snowname+" "+EXP+"_SNOW.%04d"%n)
    ct=getct('weathering.pso')
    co2s=getndco2('weathering.pso')
    dco2=getdco2('weathering.pso')
    etemps.append(ct)
    eCO2s.append(co2s[1])
    os.system("mv weathering.pso weathering"+str(n)+".pso")
    oldCO2 = pCO2
    if start:
      pCO2=eCO2s[-1]
      if np.sign(dco2) <= 0.0:
	pCO2 = pCO2*0.9
      else:
	pCO2 = pCO2*1.1
      start=False
      dtdco2 = 0.
    else:
      direction=np.sign(dco2)
      dtdco2 = abs((etemps[-1]-etemps[-2])/(eCO2s[-1]-eCO2s[-2]))
      if ct < 255.0:
	if ct > 235.0:
	  if ct > 240.0:
	    new_co2 = pCO2 + direction*min(0.5/(dtdco2+1.0e-14),2.5*pCO2) #Try to change equilibrium temperature by 0.5K
	  else:
	    new_co2 = pCO2 + direction*min(1.0/(dtdco2+1.0e-14),2.5*pCO2) #Try to change equilibrium temperature by 1.0K
	else:
	  new_co2 = pCO2 + direction*min(2.0/(dtdco2+1.0e-14),2.5*pCO2) #Try to change equilibrium temperature by 4.0K
      else:
	if ct < 273.0:
	  if ct > 268.0:
	    new_co2 = pCO2 + direction*min(1.0/(dtdco2+1.0e-14),2.5*pCO2) #Try to change equilibrium temperature by 1.0K
	  else:
	    new_co2 = pCO2 + direction*min(0.5/(dtdco2+1.0e-14),2.5*pCO2) #Try to change equilibrium temperature by 0.5K
	else:
	  new_co2 = pCO2 + direction*min(2.0/(dtdco2+1.0e-14),2.5*pCO2) #Try to change equilibrium temperature by 4.0K
      pCO2 = new_co2
    psurf = p0 + pCO2
    f=open('cyclelog.txt','a')
    f.write('Adjusted to '+str(pCO2)+'ubars CO2 and '+str(psurf*0.1)+' Pa atmosphere. dco2 was '+str(dco2)+", end co2 was "+str(co2s[1])+". Rate of temperature change was "+str(dtdco2)+" K/bar.\n")
    f.close()
    changep(psurf*0.1)
    n+=1
    changeCO2(pCO2/psurf*1e6) #Change pCO2 to CO2 ppmv
    #Adjust snowpack
    sfile1 = EXP+"_SNOW_0"
    sfile2 = EXP+"_SNOW_1"
    sfile3 = EXP+"_SNOW_2"
    sfile4 = EXP+"_SNOW_3"
    sfile5 = EXP+"_SNOW_4"
    os.system("cp newdsnow newdsnow_old")
    os.system("cp newxsnow newxsnow_old")
    deltat =  abs(pCO2-oldCO2)/(abs(dco2)+1.0e-14)
    os.system("./newsnow.x "+sfile1+" "+sfile2+" "+sfile3+" "+sfile4+" "+sfile5+" "+str(deltat))
    maxdsnow = getmaxdsnow("newdsnow_old","newdsnow")
    if maxdsnow>255.0: #Cap elevation changes to 300 meters. 2 km elevation changes may break PlaSim
      fraction=255.0/maxdsnow
      deltat*=fraction
      if ((dco2<0.) and (dco2*deltat>0.)):
	    f=open("cyclelog.txt","a")
	    f.write("WARNING TIME TRAVEL ISN'T OKAY\n")
	    f.close()
	    break
      new_co2 = oldCO2 + dco2*deltat
      pCO2 = new_co2
      psurf = p0+pCO2
      changep(psurf*0.1)
      changeCO2(pCO2/psurf*1e6)
      os.system("cp newdsnow_old newdsnow")
      os.system("cp newxsnow_old newxsnow")
      os.system("./newsnow.x "+sfile1+" "+sfile2+" "+sfile3+" "+sfile4+" "+sfile5+" "+str(deltat))
      f=open('cyclelog.txt','a')
      f.write('Adjusted to '+str(pCO2)+'ubars CO2 and '+str(psurf*0.1)+' Pa atmosphere. dco2 was '+str(dco2)+", end co2 was "+str(co2s[1])+". Rate of temperature change was "+str(dtdco2)+" K/bar. Maximum snow change was "+str(maxdsnow)+", so we capped at 85 meters.\n")
      f.close()
    os.system("mv newdsnow restart_dsnow")
    os.system("mv newxsnow restart_xsnow")
  #changeCO2(2.0e5)
