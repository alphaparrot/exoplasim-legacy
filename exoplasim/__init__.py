global sourcedir
sourcedir = None


import os
import sys
import numpy as np
import glob
import netCDF4 as nc
import exoplasim.gcmt

smws = {'mH2': 2.01588,
        'mHe': 4.002602,
        'mN2': 28.0134,
        'mO2': 31.9988,
        'mCO2':44.01,
        'mAr': 39.948,
        'mNe': 20.1797,
        'mKr': 83.798,
        'mH2O':18.01528}

gases_default = {'pH2': 0.0,
                   'pHe': 5.24e-6,
                   'pN2': 0.78084,
                   'pO2': 0.20946,
                   'pCO2':330.0e-6,
                   'pAr': 9.34e-3,
                   'pNe': 18.18e-6,
                   'pKr': 1.14e-6,
                   'pH2O':0.01}

def noneparse(text,dtype):
    if text=="None" or text=="none":
        return None
    else:
        return dtype(text)
    
#def readsourcepath():
    #with open("sourcepath","r") as sf:
        #spth = sf.read()
        #if spth.strip()=="None":
            #sourcedir=None
        #else:
            #sourcedir=spth.strip()
            
    #return sourcedir

class Model(object):
    def __init__(self,resolution="T21",layers=10,ncpus=4,precision=8,debug=False,inityear=0,
                 recompile=False,optimization=None,mars=False,workdir="most",source=None,
                 modelname="MOST_EXP"):
        '''Initialize an ExoPlaSim model configuration. By default, the current directory
           is assumed to be the working directory.'''
        
        global sourcedir
        
        if not sourcedir: #This means we haven't run yet, and have some post-install work to do
            os.system('spth=$(python3 -c "import exoplasim as exo; print(exo.__path__)") && echo $spth>sourcepath')
            with open("sourcepath","r") as spf:
                sourcedir = spf.read().strip()
                if sourcedir[0]=="[":
                    sourcedir=sourcedir[1:]
                if sourcedir[-1]=="]":
                    sourcedir=sourcedir[:-1]
                if sourcedir[0]=="'":
                    sourcedir=sourcedir[1:]
                if sourcedir[-1]=="'":
                    sourcedir=sourcedir[:-1]
            os.system("rm sourcepath")
            with open("%s/__init__.py"%sourcedir,"r") as sourcef:
                sourcecode = sourcef.read().split('\n')
            sourcecode[1] = 'sourcedir = "%s"'%sourcedir
            sourcecode = '\n'.join(sourcecode)
            os.system("cp %s/__init__.py %s/preinit.py"%(sourcedir,sourcedir))
            with open("%s/__init__.py"%sourcedir,"w") as sourcef:
                sourcef.write(sourcecode)
            cwd = os.getcwd()
            os.chdir(sourcedir)
            os.system("./configure.sh")
            os.chdir(cwd)
        
        self.runscript=None
        self.otherargs = []
        self.pgases = {}
        self.modelname=modelname
        self.cleaned=False
        
        if debug or optimization: #There is no need to set these for precompiled binaries
            recompile=True
        
        self.ncpus = ncpus
        if self.ncpus>1:
            self._exec = "mpiexec -np %d "%self.ncpus
        else:
            self._exec = "./"
        self.layers = layers
        
        self.odir = os.getcwd()
        if workdir[0]!="/":
            workdir = self.odir+"/"+workdir
        self.workdir = workdir
        os.system("mkdir "+self.workdir)
        self.currentyear=inityear
        
        # Depending on how the user has entered the resolution, set the appropriate number
        # of spectral modes and latitudes
        if resolution=="T21" or resolution=="t21" or resolution==21 or resolution==32:
            self.nsp=21
            self.nlats=32
        elif resolution=="T42" or resolution=="t42" or resolution==42 or resolution==64:
            self.nsp=42
            self.nlats=64
        elif resolution=="T63" or resolution=="t63" or resolution==63 or resolution==96:
            self.nsp=63
            self.nlats=96
        elif resolution=="T85" or resolution=="t85" or resolution==85 or resolution==128:
            self.nsp=85
            self.nlats=128
        elif resolution=="T106" or resolution=="T106" or resolution==106 or resolution==160:
            self.nsp=106
            self.nlats=160
        elif resolution=="T127" or resolution=="t127" or resolution==127 or resolution==192:
            self.nsp=127
            self.nlats=192
        elif resolution=="T170" or resolution=="t170" or resolution==170 or resolution==256:
            self.nsp=170
            self.nlats=256
        else:
            raise ValueError("Resolution unsupported. ExoPlaSim supports T21, T42, T63, T85, "+
                             "T106, T127, and T170 (32, 64, 96, 128, 160, 192, and 256 "+
                             "latitudes respectively")
        
        # If the executable does not exist, then regardless of whether we've been asked
        # to recompile, we'll have to recompile
        
        self.executable = sourcedir+"/plasim/run/most_plasim_t%d_l%d_p%d.x"%(self.nsp,self.layers,ncpus)
        
        if not source:
           source = "%s/plasim/run"%sourcedir
        
        burnsource = "%s/postprocessor"%sourcedir
        
        if recompile or not os.path.exists(self.executable):
            extraflags = ""
            if debug:
                extraflags+= "-d "
            if optimization:
                extraflags+= "-O %s"%optimization
            if mars:
                extraflags+= "-m "
            os.system("cwd=$(pwd) && "+
                      "cd %s && ./compile.sh -n %d -p %d -r T%d -v %d "%(sourcedir,self.ncpus,
                                                                        precision,self.nsp,
                                                                        self.layers)+
                      extraflags+" &&"+
                      "cd $cwd")
        
        os.system("cp %s/* %s/"%(source,self.workdir))
        os.system("cp %s/burn7.x %s/"%(burnsource,self.workdir))
        
        #Copy the executable to the working directory, and then CD there
        os.system("cp %s %s"%(self.executable,self.workdir))
        os.chdir(self.workdir)
        
        self.executable = self.executable.split("/")[-1] #Strip off all the preceding path
    
        
    def run(self,**kwargs):
        if not self.runscript:
            self._run(**kwargs)
        else:
            try:
                self.runscript(self,**kwargs) #runscript MUST accept a Model object as the first arg
            except Exception as e:
                print(e)
                self.crash()
        
    def runtobalance(self,threshhold = 4.0e-4,baseline=50,maxyears=300,minyears=75):
        runlimit = self.currentyear+maxyears
        if os.getcwd()!=self.workdir:
            os.chdir(self.workdir)
        os.system("mkdir snapshots")
        if self.highcadence["toggle"]:
            os.system("mkdir highcadence")
        while not self._isbalanced(threshhold=threshhold,baseline=baseline) \
                and self.currentyear<minyears and self.currentyear<runlimit:
            dataname="MOST.%05d"%self.currentyear
            snapname="MOST_SNAP.%05d"%self.currentyear
            hcname  ="MOST_HC.%05d"%self.currentyear
            diagname="MOST_DIAG.%05d"%self.currentyear
            restname="MOST_REST.%05d"%self.currentyear
            snowname="MOST_SNOW.%05d"%self.currentyear
            
            #Run ExoPlaSim
            os.system(self._exec+self.executable)
            
            #Sort, categorize, and arrange the various outputs
            os.system("[ -e restart_dsnow ] && rm restart_dsnow")
            os.system("[ -e restart_xsnow ] && rm restart_xsnow")
            os.system("[ -e Abort_Message ] && exit 1")
            os.system("[ -e plasim_output ] && mv plasim_output "+dataname)
            os.system("[ -e plasim_snapshot ] && mv plasim_snapshot "+snapname)
            if self.highcadence["toggle"]:
                os.system("[ -e plasim_hcadence ] && mv plasim_hcadence "+hcname)
            os.system("[ -e plasim_diag ] && mv plasim_diag "+diagname)
            os.system("[ -e plasim_status ] && cp plasim_status plasim_restart")
            os.system("[ -e plasim_status ] && mv plasim_status "+restname)
            os.system("[ -e restart_snow ] && mv restart_snow "+snowname)
            
            #Do any additional work
            timeavg=0
            snapsht=0
            highcdn=0
            try:
                timeavg=self.postprocess(dataname,"example.nl",
                                        log="burnout",crashifbroken=crashifbroken)
                snapsht=self.postprocess(snapname,"snapshot.nl",
                                        log="snapout",crashifbroken=crashifbroken)
                os.system("mv %s.nc snapshots/"%snapname)
                if self.highcadence["toggle"]:
                    highcdn=self.postprocess(hcname  ,"snapshot.nl",
                                            log="hcout"  ,crashifbroken=crashifbroken)
                    os.system("mv %s.nc highcadence/"%hcname)
            except Exception as e:
                print(e)
                self.crash()
            if clean:
                if timeavg:
                    os.system("rm %s"%dataname)
                if snapsht:
                    os.system("rm %s"%snapname)
                if highcdn:
                    os.system("rm %s"%hcname)
                    
            if crashifbroken: #Check to see that we aren't throwing NaNs
                try:
                    check=self.integritycheck(dataname+".nc")
                except Exception as e:
                    print(e)
                    self.crash()
                
            self.currentyear += 1
            sb,tb = self.getbalance()
            os.system("echo %02.6f  %02.6f'>>%s/balance.log"%(sb,tb,self.workdir))
        bott = self.gethistory(key="hfns")
        topt = self.gethistory(key="ntr")
        with open("%s/shistory.pso"%self.workdir,"a+") as f:
            text='\n'+'\n'.join(bott.astype(str))
            f.write(text)
        with open("%s/toahistory.pso"%self.workdir,"a+") as f:
            text='\n'+'\n'.join(topt.astype(str))
            f.write(text)
        finished = self._isbalanced(threshhold=threshhold,baseline=baseline)
        if not finished:
            return False
        return True
    
            
    def getbalance(self,year=-1):
        topt = self.inspect("ntr",savg=True,tavg=True,year=year)
        bott = self.inspect("hfns",savg=True,tavg=True,year=year)
        return topt,bott
    
    def gethistory(self,key="ts",mean=True,layer=-1):
        files = sorted(glob.glob("%s/MOST*.nc"%self.workdir))
        dd=np.zeros(len(files))
        for n in range(0,len(files)):
            ncd = self.get(n)
            variable = ncd.variables[key][:]
            lon = ncd.variables['lon'][:]
            lat = ncd.variables['lat'][:]
            if len(variable.shape)>3:
                variable = variable[:,layer,:,:]
            dd[n] = gcmt.spatialmath(variable,lon=lon,lat=lat,
                                     mean=mean,radius=self.radius)
            ncd.close()
        return dd
    
    
    def _isbalanced(self,threshold = 1.0e-4,baseline=50):
        nfiles = len((glob.glob("%s/MOST*.nc"%self.workdir)))
        prior=False
        if len(glob.glob(self.workdir+"/toahistory.ps*"))>0:
            toahistory = np.loadtxt(self.workdir+"/toahistory.pso")
            nfiles+=len(toahistory)
            shistory = np.loadtxt(self.workdir+"/shistory.pso")
            prior=True
        sbalance = np.zeros(nfiles)
        toabalance=np.zeros(nfiles)
        nstart=0
        if prior:
            sbalance[:len(toahistory)] = shistory[:]
            toabalance[:len(toahistory)] = toahistory[:]
            nstart = len(toahistory)
        if self.currentyear < baseline: #Run for minimum of baseline years
            return False
        else:
            for n in range(0,self.currentyear):
                topt,bott = self.getbalance(year=n)
                sbalance[n+nstart] = bott
                toabalance[n+nstart] = topt
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
            os.system("echo '%02.8f  %02.8f'>>%s/slopes.log"%(savgslope,tavgslope,self.workdir))
            if savgslope<threshhold and tavgslope<threshhold: #Both TOA and Surface are changing at average 
                return True                                  # of <0.1 mW/m^2/yr on 45-year baselines
            else:
                return False
        
    def _run(self,years=1,postprocess=True,crashifbroken=False,clean=True):
        if os.getcwd()!=self.workdir:
            os.chdir(self.workdir)
        os.system("mkdir snapshots")
        if self.highcadence["toggle"]:
            os.system("mkdir highcadence")
        for year in range(years):
            dataname="MOST.%05d"%self.currentyear
            snapname="MOST_SNAP.%05d"%self.currentyear
            hcname  ="MOST_HC.%05d"%self.currentyear
            diagname="MOST_DIAG.%05d"%self.currentyear
            restname="MOST_REST.%05d"%self.currentyear
            snowname="MOST_SNOW.%05d"%self.currentyear
            
            #Run ExoPlaSim
            os.system(self._exec+self.executable)
            
            #Sort, categorize, and arrange the various outputs
            os.system("[ -e restart_dsnow ] && rm restart_dsnow")
            os.system("[ -e restart_xsnow ] && rm restart_xsnow")
            os.system("[ -e Abort_Message ] && exit 1")
            os.system("[ -e plasim_output ] && mv plasim_output "+dataname)
            os.system("[ -e plasim_snapshot ] && mv plasim_snapshot "+snapname)
            if self.highcadence["toggle"]:
                os.system("[ -e plasim_hcadence ] && mv plasim_hcadence "+hcname)
            os.system("[ -e plasim_diag ] && mv plasim_diag "+diagname)
            os.system("[ -e plasim_status ] && cp plasim_status plasim_restart")
            os.system("[ -e plasim_status ] && mv plasim_status "+restname)
            os.system("[ -e restart_snow ] && mv restart_snow "+snowname)
            
            #Do any additional work
            timeavg=0
            snapsht=0
            highcdn=0
            if postprocess:
                try:
                    timeavg=self.postprocess(dataname,"example.nl",
                                            log="burnout",crashifbroken=crashifbroken)
                    snapsht=self.postprocess(snapname,"snapshot.nl",
                                            log="snapout",crashifbroken=crashifbroken)
                    os.system("mv %s.nc snapshots/"%snapname)
                    if self.highcadence["toggle"]:
                        highcdn=self.postprocess(hcname  ,"snapshot.nl",
                                                log="hcout"  ,crashifbroken=crashifbroken)
                        os.system("mv %s.nc highcadence/"%hcname)
                except Exception as e:
                    print(e)
                    self.crash()
            if clean:
                if timeavg:
                    os.system("rm %s"%dataname)
                if snapsht:
                    os.system("rm %s"%snapname)
                if highcdn:
                    os.system("rm %s"%hcname)
                    
            if crashifbroken: #Check to see that we aren't throwing NaNs
                try:
                    check=self.integritycheck(dataname+".nc")
                except Exception as e:
                    print(e)
                    self.crash()
                
            self.currentyear += 1
    
    def postprocess(self,inputfile,namelist,log="postprocess.log",crashifbroken=False):
        stat=os.system("./burn7.x -n<%s>%s %s %s.nc"%(namelist,log,inputfile,inputfile))
        if stat==0:
            print("NetCDF output written to %s.nc; log written to %s"%(inputfile,log))
            return 1
        else:
            if crashifbroken:
                raise RuntimeError("Error writing output to %s.nc; "%inputfile +
                                   "log written to %s"%log)
            else:
                print("Error writing output to %s.nc; log written to %s"%(inputfile,log))
            return 0
        
    def integritycheck(self,ncfile): #MUST pass a NetCDF file that contains surface temperature
        if os.getcwd()!=self.workdir:
            os.chdir(self.workdir)
        ioe=1
        if not os.path.exists(ncfile): #If the specified output file does not exist, create it
            ioe = self.postprocess(ncfile[:-3],"example.nl",crashifbroken=True)
        if ioe:
            ncd = nc.Dataset(ncfile,"r")
            try:
                ts = ncd.variables["ts"][:]
            except:
                raise RuntimeError("Output is missing surface temperature; check logs for errors")
            if np.sum(np.isnan(ts))+np.sum(np.isinf(ts)) > 0.5:
                raise RuntimeError("Non-finite values found in surface temperature")
            return 1
        else:
            return 0
                
    def finalize(self,outputdir,allyears=False,keeprestarts=False,clean=True):
        if outputdir[0]!="/":
            outputdir = os.getcwd()+"/"+outputdir
        if not os.path.isdir(outputdir):
            os.system("mkdir %s"%outputdir)
        if allyears:
            os.chdir(outputdir)
            os.system("mkdir %s"%self.modelname)
            os.system("cp %s/MOST*.nc %s/"%(self.workdir,self.modelname))
            if self.snapshots:
                os.system("cp -r %s/snapshots %s/snapshots"%(self.workdir,self.modelname))
            if self.highcadence['toggle']:
                os.system("cp -r %s/highcadence %s/highcadence"%(self.workdir,self.modelname))
            os.system("cp %s/MOST*DIAG* %s/"%(self.workdir,self.modelname))
            if keeprestarts:
                os.system("cp %s/MOST_REST* %s/"%(self.workdir,self.modelname))
            #else:
            #    restarts = sorted(glob.glob("%s/MOST_REST*"%self.workdir))
            #    os.system("cp %s %s/%s_restart"%(restarts[-1],
            #                                        self.modelname,self.modelname))
            newworkdir = os.getcwd()+"/"+self.modelname
        else:
            outputs = sorted(glob.glob("%s/MOST*.nc"%self.workdir))
            os.chdir(outputdir)
            os.system("cp %s %s.nc"%(outputs[-1],self.modelname))
            diags = sorted(glob.glob("%s/MOST*DIAG*"%self.workdir))
            os.system("cp %s %s.DIAG"%(diags[-1],self.modelname))
            if self.snapshots:
                snps = sorted(glob.glob("%s/snapshots/*.nc"%self.workdir))
                os.system("cp %s %s_snapshot.nc"%(snps[-1],self.modelname))
            if self.highcadence["toggle"]:
                hcs = sorted(glob.glob("%s/highcadence/MOST*.nc"%self.workdir))
                os.system("cp %s %s_highcadence.nc"%(hcs[-1],self.modelname))
            if keeprestarts:
                rsts = sorted(glob.glob("%s/MOST_REST*"%self.workdir))
                os.system("cp %s %s_restart"%(rsts[-1],self.modelname))
            if clean:
                newworkdir = os.getcwd()
                self.cleaned=True
        if clean:
            os.system("rm -rf %s"%self.workdir)
            self.workdir = newworkdir
                
    
    def get(self,year,snapshot=False,highcadence=False):
        #Note: if the work directory has been cleaned out, only the final year will be returned.
        if snapshot and not highcadence:
            name = "snapshots/MOST_SNAP.%05d.nc"%year
        elif highcadence and not snapshot:
            name = "highcadence/MOST_HC.%05d.nc"%year
        else:
            name = "MOST.%05d.nc"%year
        if self.cleaned:
            if snapshot and not highcadence:
                name = "%s_snapshot.nc"%self.modelname
            elif highcadence and not snapshot:
                name = "%s_highcadence.nc"%self.modelname
            else:
                name = "%s.nc"%self.modelname
        if os.path.exists(self.workdir+"/"+name):
            ncd = nc.Dataset(self.workdir+"/"+name,"r")
            return ncd
        else:
            raise RuntimeError("Output file %s not found."%(self.workdir+"/"+name))
    
    def inspect(self,variable,year=-1,ignoreNaNs=True,snapshot=False,
                highcadence=False,savg=False,tavg=False,layer=None):
        #Note: if the work directory has been cleaned out, only the final year will be returned.
        if snapshot and not highcadence:
            pattern = "snapshots/MOST_SNAP"
        elif highcadence and not snapshot:
            pattern = "highcadence/MOST_HC"
        else:
            pattern = "MOST"
        
        if ignoreNaNs:
            meanop = np.nanmean
        else:
            meanop = np.mean
        
        if year<0:
            nfiles = len(glob.glob(self.workdir+"/"+pattern+"*.nc"))
            year = nfiles+year
    
        ncd = self.get(year,snapshot=snapshot,highcadence=highcadence)
        
        var = ncd.variables[variable][:]
        
        if variable!="lat" and variable!="lon" and variable!="lev" and variable!="time":
            lon = ncd.variables['lon'][:]
            lat = ncd.variables['lat'][:]
            lev = ncd.variables['lev'][:]
            if not savg and not tavg:
                return var
            elif tavg and not savg:
                return meanop(var,axis=0)
            elif tavg and savg:
                if type(layer)!=type(None) and len(var.shape)==4: #3D spatial array, plus time
                #We're going to get a scalar; user has specified a level
                    return gcmt.spatialmath(var,lat=lat,lon=lon,ignoreNaNs=ignoreNaNs,lev=layer)
                elif type(layer)==type(None) and len(var.shape)==4:
                    #We're going to get a 1D vertical profile where each layer is a spatial avg
                    output = np.zeros(lev.shape)
                    for l in range(len(lev)):
                        output[l] = gcmt.spatialmath(var,lat=lat,lon=lon,
                                                     ignoreNaNs=ignoreNaNs,lev=l)
                    return output
                elif len(var.shape)==3: #2D spatial array, plus time
                    return gcmt.spatialmath(var,lat=lat,lon=lon,ignoreNaNs=ignoreNaNs)
            else:
                if type(layer)!=type(None) and len(var.shape)==4: #3D spatial array, plus time
                #We're going to get a 1D array; user has specified a level
                    output = np.zeros(var.shape[0])
                    for t in range(len(output)):
                        output[t] = gcmt.spatialmath(var,lat=lat,lon=lon,ignoreNaNs=ignoreNaNs,
                                                     lev=layer,time=t)
                    return output
                elif type(layer)==type(None) and len(var.shape)==4:
                    #We're going to get a 1D vertical profile plus time
                    output = np.zeros((var.shape[0],len(lev)))
                    for t in range(var.shape[0]):
                        for l in range(len(lev)):
                            output[t,l] = gcmt.spatialmath(var,lat=lat,lon=lon,time=t,
                                                     ignoreNaNs=ignoreNaNs,lev=l)
                    return output
                elif len(var.shape)==3: #2D spatial array, plus time
                    #We're going to get a 1D array
                    output = np.zeros(var.shape[0])
                    for t in range(len(output)):
                        output[t] = gcmt.spatialmath(var,lat=lat,lon=lon,
                                                     time=t,ignoreNaNs=ignoreNaNs)
                    return output
                else:
                    return -1
        else:
            return var
                
    
    def crash(self):
        os.chdir(self.workdir)
        os.chdir("..")
        os.system("mv %s %s_crashed"%(self.workdir,os.getcwd()+"/"+self.modelname))
        raise RuntimeError("ExoPlaSim has crashed or begun producing garbage. All working files have been moved to %s_crashed/"%(os.getcwd()+"/"+self.modelname))
        
    def configure(self,noutput=True,flux=1367.0,startemp=None,starspec=None,pH2=None,
               pHe=None,pN2=None,pO2=None,pCO2=None,pAr=None,pNe=None,
               pKr=None,pH2O=None,gascon=None,pressure=None,pressurebroaden=True,
               vtype=0,rotationperiod=24.0,synchronous=False,substellarlon=180.0,
               restartfile=None,gravity=9.80665,radius=1.0,eccentricity=None,
               obliquity=None,lonvernaleq=None,fixedorbit=False,orography=None,
               seaice=True,co2weathering=False,evolveco2=False,physicsfilter=None,
               filterkappa=8.0,filterpower=8,filterLHN0=15.0,diffusionwaven=None,
               qdiffusion=None,tdiffusion=None,zdiffusion=None,ddiffusion=None,
               diffusionpower=None,erosionsupplylimit=None,outgassing=50.0,snowicealbedo=None,
               twobandalbedo=False,maxsnow=None,soilalbedo=None,oceanalbedo=None,
               oceanzenith="ECHAM-3",wetsoil=False,soilwatercap=None,aquaplanet=False,
               desertplanet=False,soilsaturation=None,drycore=False,ozone=True,
               cpsoil=None,soildepth=1.0,mldepth=50.0,
               writefrequency=None,modeltop=None,stratosphere=False,
               tropopause=None,timestep=45.0,runscript=None,columnmode=None,
               highcadence={"toggle":0,"start":320,"end":576,"interval":4},
               snapshots=None,resources=[],landmap=None,stormclim=False,nstorms=4,
               stormcapture={"VITHRESH":0.145,"GPITHRESH":0.37,"VMXTHRESH":33.0,
                             "LAVTHRESH":1.2e-5,"VRMTHRESH":0.577,"MINSURFTEMP":298.15,
                             "MAXSURFTEMP":373.15,"WINDTHRESH":33.0,"SWINDTHRESH":20.5,
                             "SIZETHRESH":30,"ENDTHRESH":16,"MINSTORMLEN":256,
                             "MAXSTORMLEN":1024,"NKTRIGGER":0,"toggle":0},
               topomap=None,otherargs={}):
        
        self._edit_namelist("plasim_namelist","NOUTPUT",str(noutput*1))
        self.noutput = noutput
        self._edit_namelist("planet_namelist","GSOL0",str(flux))
        self.flux = flux
        if startemp:
            self._edit_namelist("radmod_namelist","NSTARTEMP","1")
            self._edit_namelist("radmod_namelist","STARBBTEMP",str(startemp))
        self.startemp = startemp
        if starspec:
            self._edit_namelist("radmod_namelist","NSTARFILE","1")
            self._edit_namelist("radmod_namelist","STARFILE","'%s'"%starspec)
            self._edit_namelist("radmod_namelist","STARFILEHR","'%s_hr.dat'"%(starspec)[:-4])
        self.starspec = starspec
        
        if pH2:
            self.pgases["pH2"]=pH2
        if pHe:
            self.pgases["pHe"]=pHe
        if pN2:
            self.pgases["pN2"]=pN2
        if pO2:
            self.pgases["pO2"]=pO2
        if pCO2:
            self.pgases["pCO2"]=pCO2
        if pAr:
            self.pgases["pAr"]=pAr
        if pNe:
            self.pgases["pNe"]=pNe
        if pKr:
            self.pgases["pKr"]=pKr
        if pH2O:
            self.pgases["pH2O"]=pH2O
        
        if len(self.pgases)==0:
            if not pressure:
                #self.pgases=gases_default
                #self.pressure=0.0
                #for gas in self.pgases:
                    #self.pressure+=self.pgases[gas] #in bars
                self.pressure = 1.011
            else:
                self.pressure = pressure
        
        else:
            self.pressure=0.0
            for gas in self.pgases:
                self.pressure+=self.pgases[gas]
        
        gasesvx = {}
        for gas in self.pgases:
            gasesvx[gas[1:]] = self.pgases[gas]/self.pressure
        self.mmw = 0
        for gas in gasesvx:
            self.mmw += gasesvx[gas]*smws['m'+gas]
        if self.mmw==0:
            self.mmw = 28.970253
        self.gascon = 8314.46261815324 / self.mmw
        
        if gascon:
            self.gascon=gascon
            self.mmw = 8314.46261815324/self.gascon
        
        print('Mean Molecular Weight set to %1.4f g/mol'%self.mmw)
        
        if pressure:
            if pressure != self.pressure:  #User has specified a different pressure than sum of gas pressures
                self.pressure = pressure
        
        self.CO2ppmv = 0.0
        if 'pCO2' in self.pgases:
            self.CO2ppmv = self.pgases['pCO2']/self.pressure * 1.0e6 #ppmv
        
        self._edit_namelist("radmod_namelist","CO2",str(self.CO2ppmv))
        self._edit_namelist("plasim_namelist","PSURF",str(self.pressure*1.0e5)) #Pa
        self._edit_namelist("planet_namelist","GASCON",str(self.gascon))
        self._edit_namelist("radmod_namelist","NPBROADEN",str(pressurebroaden*1))
        self.pressurebroaden=pressurebroaden
        self._edit_namelist("plasim_namelist","NEQSIG",str(vtype))
        self.vtype=vtype
        if rotationperiod!=1.0:
            self._edit_namelist("planet_namelist","ROTSPD",str(1.0/float(rotationperiod)))
            self._edit_namelist("plasim_namelist","N_DAYS_PER_YEAR",
                                str(max(int(360.0/float(rotationperiod)/12+0.5),1)*12))
        self.rotationperiod=rotationperiod
        if synchronous:
            self._edit_namelist("radmod_namelist","NFIXED","1")
            self._edit_namelist("radmod_namelist","FIXEDLON",str(substellarlon))
        self.synchronous=synchronous
        self.substellarlon=substellarlon
        
        if restartfile:
            os.system("cp %s %s/plasim_restart"%(restartfile,self.workdir))
        else:
            os.system("rm %s/plasim_restart"%self.workdir)
        self.restartfile=restartfile
        
        self._edit_namelist("planet_namelist","GA",str(gravity))
        self._edit_postnamelist("example.nl","gravity",str(gravity))
        self._edit_postnamelist("snapshot.nl","gravity",str(gravity))
        self.gravity=gravity
        
        self._edit_namelist("planet_namelist","PLARAD",str(radius*6371220.0))
        self._edit_postnamelist("example.nl","radius",str(radius*6371220.0))
        self._edit_postnamelist("snapshot.nl","radius",str(radius*6371220.0))
        self.radius=radius
        
        if type(eccentricity)!=type(None):
            self._edit_namelist("planet_namelist","ECCEN",str(eccentricity))
        self.eccentricity=eccentricity
        
        if type(obliquity)!=type(None):
            self._edit_namelist("planet_namelist","OBLIQ",str(obliquity))
        self.obliquity=obliquity
        
        if type(lonvernaleq)!=type(None):
            self._edit_namelist("planet_namelist","MVELP",str(lonvernaleq))
        self.lonvernaleq=lonvernaleq
        
        self._edit_namelist("planet_namelist","NFIXORB",str(fixedorbit*1))
        self.fixedorbit=fixedorbit
        
        if type(orography)!=type(None):
            self._edit_namelist("landmod_namelist","OROSCALE",str(orography))
            self._edit_namelist("glacier_namelist","NGLACIER","1")
        self.orography=orography
                
        self._edit_namelist("radmod_namelist","NRADICE",str(seaice*1))
        self.seaice=seaice
        
        self._edit_namelist("carbonmod_namelist","NCARBON",str(co2weathering*1))
        self._edit_namelist("carbonmod_namelist","NCO2EVOLVE",str(evolveco2*1))
        self.co2weathering=co2weathering
        self.evolveco2=evolveco2
        if type(erosionsupplylimit)!=type(None):
            self._edit_namelist("carbonmod_namelist","NSUPPLY","1")
            self._edit_namelist("carbonmod_namelist","WMAX",str(erosionsupplylimit))
        self.erosionsupplylimit=erosionsupplylimit
        self._edit_namelist("carbonmod_namelist","VOLCANCO2",str(outgassing))
        self.outgassing=outgassing
        
        if physicsfilter:
            vals = physicsfilter.split("|")
            if "gp" in vals:
                self._edit_namelist("plasim_namelist","NGPTFILTER","1")
            if "sp" in vals:
                self._edit_namelist("plasim_namelist","NSPVFILTER","1")
            if "none" in vals:
                self._edit_namelist("plasim_namelist","NFILTER","0")
            if "cesaro" in vals:
                self._edit_namelist("plasim_namelist","NFILTER","1")
            if "exp" in vals:
                self._edit_namelist("plasim_namelist","NFILTER","2")
            if "lh" in vals:
                self._edit_namelist("plasim_namelist","NFILTER","3")
        self.physicsfilter=physicsfilter
        self._edit_namelist("plasim_namelist","FILTERKAPPA",str(filterkappa))
        self.filterkappa=filterkappa
        self._edit_namelist("plasim_namelist","NFILTEREXP",str(filterpower))
        self.filterpower=filterpower
        self._edit_namelist("plasim_namelist","LANDHOSKN0",str(filterLHN0))
        self.filterLHN0=filterLHN0
        if diffusionwaven:
            self._edit_namelist("plasim_namelist","NHDIFF",str(diffusionwaven))
        self.diffusionwaven=diffusionwaven
        if qdiffusion:
            self._edit_namelist("plasim_namelist","TDISSQ","%d*%f"%(self.layers,qdiffusion))
        self.qdiffusion=qdiffusion
        if tdiffusion:
            self._edit_namelist("plasim_namelist","TDISST","%d*%f"%(self.layers,tdiffusion))
        self.tdiffusion=tdiffusion
        if zdiffusion:
            self._edit_namelist("plasim_namelist","TDISSZ","%d*%f"%(self.layers,zdiffusion))
        self.zdiffusion=zdiffusion
        if ddiffusion:
            self._edit_namelist("plasim_namelist","TDISSD","%d*%f"%(self.layers,ddiffusion))
        self.ddiffusion=ddiffusion
        if diffusionpower:
            self._edit_namelist("plasim_namelist","NDEL","%d*%d"%(self.layers,diffusionpower))
        self.diffusionpower=diffusionpower
        
        if type(snowicealbedo)!=type(None):
            alb = str(snowicealbedo)
            self._edit_namelist("seamod_namelist","ALBICE",alb)
            self._edit_namelist("seamod_namelist","DICEALBMN","%s,%s"%(alb,alb))
            self._edit_namelist("seamod_namelist","DICEALBMX","%s,%s"%(alb,alb))
            self._edit_namelist("landmod_namelist","DSNOWALBMN","%s,%s"%(alb,alb))
            self._edit_namelist("landmod_namelist","DSNOWALBMN","%s,%s"%(alb,alb))
            self._edit_namelist("landmod_namelist","DGLACALBMN","%s,%s"%(alb,alb))
            self._edit_namelist("landmod_namelist","DSNOWALB","%s,%s"%(alb,alb))
        self.snowicealbedo=snowicealbedo
        
        self._edit_namelist("radmod_namelist","NSIMPLEALBEDO",str((not twobandalbedo)*1))
        self.twobandalbedo=twobandalbedo
        
        if maxsnow:
            self._edit_namelist("landmod_namelist","DSMAX",str(maxsnow))
        self.maxsnow=maxsnow
        
        if type(soilalbedo)!=type(None):
            alb = str(soilalbedo)
            os.system("rm %s/*0174.sra"%self.workdir)
            self._edit_namelist("landmod_namelist","ALBLAND",alb)
            self._edit_namelist("landmod_namelist","DGROUNDALB","%s,%s"%(alb,alb))
        self.soilalbedo=soilalbedo
        
        if type(oceanalbedo)!=type(None):
            alb = str(oceanalbedo)
            self._edit_namelist("seamod_namelist","ALBSEA",alb)
            self._edit_namelist("seamod_namelist","DOCEANALB","%s,%s"%(alb,alb))
        self.oceanalbedo=oceanalbedo
        
        if oceanzenith=="lambertian" or oceanzenith=="Lambertian" or oceanzenith=="uniform":
            self._edit_namelist("radmod_namelist","NECHAM","0")
            self._edit_namelist("radmod_namelist","NECHAM6","0")
        if oceanzenith=="default" or oceanzenith=="plasim" or oceanzenith=="ECHAM-3":
            self._edit_namelist("radmod_namelist","NECHAM","1")
            self._edit_namelist("radmod_namelist","NECHAM6","0")
        if oceanzenith=="ECHAM-6":
            self._edit_namelist("radmod_namelist","NECHAM","0")
            self._edit_namelist("radmod_namelist","NECHAM6","1")
        self.oceanzenith=oceanzenith
        
        self._edit_namelist("landmod_namelist","NWETSOIL",str(wetsoil*1))
        self.wetsoil=wetsoil
        if soilwatercap:
            self._edit_namelist("landmod_namelist","WSMAX",str(soilwatercap))
            os.system("rm %s/*0229.sra"%self.workdir)
        if desertplanet:
            self._edit_namelist("plasim_namelist","NDESERT","1")
            self._edit_namelist("landmod_namelist","NWATCINI","1")
            self._edit_namelist("landmod_namelist","DWATCINI","0.0")
            os.system("rm %s/*.sra"%self.workdir)
        if type(soilsaturation)!=type(None):
            self._edit_namelist("landmod_namelist","NWATCINI","1")
            self._edit_namelist("landmod_namelist","DWATCINI",str(soilsaturation))
            os.system("rm %s/*0229.sra"%self.workdir)
        if aquaplanet:
            self._edit_namelist("plasim_namelist","NAQUA","1")
            os.system("rm %s/*.sra"%self.workdir)
        self.soilwatercap=soilwatercap
        self.soilsaturation=soilsaturation
        self.desertplanet=desertplanet
        self.aquaplanet=aquaplanet
        
        if drycore:
            self._edit_namelist("fluxmod_namelist","NEVAP","0")
        self.drycore=drycore
        
        if columnmode:
            parts = columnmode.split("|")
            if "static" in parts:
                self._edit_namelist("plasim_namelist","NADV","0")
            if "clear" in parts:
                self._edit_namelist("radmod_namelist","NCLOUDS","0")
                self._edit_namelist("radmod_namelist","ACLLWR","0.0")
        self.columnmode=columnmode
        
        self._edit_namelist("radmod_namelist","NO3",str(ozone*1))
        self.ozone=ozone
        
        if cpsoil:
            self._edit_namelist("landmod_namelist","SOILCAP",str(cpsoil))
        self.cpsoil = cpsoil
        
        self.dzsoils = np.array([0.4, 0.8, 1.6, 3.2, 6.4])*soildepth
        self._edit_namelist("landmod_namelist","DSOILZ",",".join(self.dzsoils.astype(str)))
        self.soildepth=soildepth
        
        self._edit_namelist("oceanmod_namelist","MLDEPTH",str(mldepth))
        self.mldepth=mldepth
        
        if writefrequency:
            self._edit_namelist("plasim_namelist","NWPD",str(writefrequency))
        self.writefrequency=writefrequency
        
        if stratosphere:
            self._edit_namelist("plasim_namelist","NEQSIG","5")
            if modeltop:
                self._edit_namelist("plasim_namelist","PTOP2",str(modeltop*100.0)) #convert hPa->Pa
            if tropopause:
                self._edit_namelist("plasim_namelist","PTOP",str(tropopause*100.0))
        else:
            if modeltop:
                self._edit_namelist("plasim_namelist","PTOP",str(modeltop*100.0)) #convert hPa->Pa
        self.stratosphere=stratosphere
        self.modeltop=modeltop
        self.tropopause=tropopause
        
        self._edit_namelist("plasim_namelist","MPSTEP",str(timestep))
        self._edit_namelist("plasim_namelist","NSTPW",str(int(7200//int(timestep))))
        self.timestep=timestep
        
        self.runscript=runscript
        
        self._edit_namelist("plasim_namelist","NHCADENCE",str(highcadence["toggle"]))
        self._edit_namelist("plasim_namelist","HCSTARTSTEP",str(highcadence["start"]))
        self._edit_namelist("plasim_namelist","HCENDSTEP",str(highcadence["end"]))
        self._edit_namelist("plasim_namelist","HCINTERVAL",str(highcadence["interval"]))
        self.highcadence=highcadence
        
        if snapshots:
            self._edit_namelist("plasim_namelist","NSNAPSHOT","1")
            self._edit_namelist("plasim_namelist","NSTPS",str(snapshots))
        self.snapshots=snapshots
        
        if len(resources)>0:
            for res in resources:
                os.system("cp %s %s/"%(res,self.workdir))
        self.resources=resources
        
        if landmap or topomap:
            os.system("rm %s/*.sra"%self.workdir)
        if landmap:
            os.system("cp %s %s/N%03d_surf_0172.sra"%(landmap,self.workdir,self.nlats))
        if topomap:
            os.system("cp %s %s/N%03d_surf_0129.sra"%(topomap,self.workdir,self.nlats))
        self.landmap=landmap
        self.topomap=topomap
        
        if stormclim:
            self._edit_namelist("hurricane_namelist","NSTORMDIAG","1")
            self._add_postcodes("example.nl",[322,323,324,325,326,327,328,329])
            self._add_postcodes("snapshot.nl",[322,323,324,325,326,327,328,329])
        self.stormclim=stormclim
        self._edit_namelist("hurricane_namelist","NSTORMS",str(int(nstorms)))
        self.nstorms=nstorms
        if stormcapture["toggle"]:
            self._edit_namelist("hurricane_namelist","HC_CAPTURE","1")
            for param in stormcapture:
                if param!="toggle":
                    self._edit_namelist("hurricane_namelist",param,str(stormcapture[param]))
        self.stormcapture=stormcapture
        
        if len(otherargs)>0:
            for key in otherargs:
                parts = key.split("<")
                destination=parts[0].split("@")
                value=parts[1]
                field=destination[0]
                namelist=destination[1]
                self._edit_namelist(namelist,field,value)
                self.otherargs.append(key)
     
    def loadconfig(self,configfile):
        with open(configfile,"r") as cfgf:
            cfg = cfgf.read().split("\n")
        noutput=bool(int(cfg[0]))
        flux=float(cfg[1])
        startemp=noneparse(cfg[2],float)
        starspec=noneparse(cfg[3],str)
        gases = cfg[4].split("&")
        for gas in gases:
            species = gas.split("|")
            amt = float(species[1])
            species = species[0]
            self.pgases[species] = amt
        gascon = float(cfg[5])
        pressure = noneparse(cfg[6],float)
        pressurebroaden = bool(int(cfg[7]))
        vtype = int(cfg[8])
        rotationperiod = float(cfg[9])
        synchronous = bool(int(cfg[10]))
        substellarlon = float(cfg[11])
        restartfile = noneparse(cfg[12],str)
        gravity = float(cfg[13])
        radius = float(cfg[14])
        eccentricity = noneparse(cfg[15],float)
        obliquity = noneparse(cfg[16],float)
        lonvernaleq = noneparse(cfg[17],float)
        fixedorbit = bool(int(cfg[18]))
        orography = noneparse(cfg[19],float)
        seaice = bool(int(cfg[20]))
        co2weathering = bool(int(cfg[21]))
        evolveco2 = bool(int(cfg[22]))
        physicsfilter = noneparse(cfg[23],str)
        filterkappa = float(cfg[24])
        filterpower = int(cfg[25])
        filterLHN0 = float(cfg[26])
        diffusionwaven = noneparse(cfg[27],int)
        qdiffusion = noneparse(cfg[28],float)
        tdiffusion = noneparse(cfg[29],float)
        zdiffusion = noneparse(cfg[30],float)
        ddiffusion = noneparse(cfg[31],float)
        diffusionpower = noneparse(cfg[32],int)
        erosionsupplylimit = noneparse(cfg[33],float)
        outgassing = float(cfg[34])
        snowicealbedo = noneparse(cfg[35],float)
        twobandalbedo = bool(int(cfg[36]))
        maxsnow = noneparse(cfg[37],float)
        soilalbedo = noneparse(cfg[38],float)
        oceanalbedo = noneparse(cfg[39],float)
        oceanzenith = cfg[40]
        wetsoil = bool(int(cfg[41]))
        soilwatercap = noneparse(cfg[42],float)
        aquaplanet = bool(int(cfg[43]))
        desertplanet = bool(int(cfg[44]))
        soilsaturation = noneparse(cfg[45],float)
        drycore = bool(int(cfg[46]))
        ozone = bool(int(cfg[47]))
        cpsoil = noneparse(cfg[48],float)
        soildepth = float(cfg[49])
        mldepth = float(cfg[50])
        writefrequency = noneparse(cfg[51],int)
        modeltop = noneparse(cfg[52],float)
        stratosphere = bool(int(cfg[53]))
        tropopause = noneparse(cfg[54],float)
        timestep = float(cfg[55])
        runscript = noneparse(cfg[56],str)
        columnmode = noneparse(cfg[57],str)
        hcdict = cfg[58].split("&")
        highcadence = {}
        for hc in hcdict:
            parts = hc.split("|")
            highcadence[parts[0]] = int(parts[1])
        snapshots = noneparse(cfg[59],int)
        resources = []
        reslist = cfg[60].split("&")
        if len(reslist)>0 and reslist[0]!='':
            for res in reslist:
                resources.append(res)
        landmap = noneparse(cfg[61],str)
        topomap = noneparse(cfg[62],str)
        stormclim = bool(int(cfg[63]))
        nstorms = int(cfg[64])
        stormcapture = {}
        stormdict = cfg[65].split("&")
        for item in stormdict:
            parts = item.split("|")
            if parts[0]=="toggle" or parts[0]=="NKTRIGGER" or parts[0]=="SIZETHRESH" \
             or parts[0]=="ENDTHRESH" or parts[0]=="MINSTORMLEN" or parts[0]=="MAXSTORMLEN":
                stormcapture[parts[0]] = int(parts[1])
            else:
                stormcapture[parts[0]] = float(parts[1])
        otherargs = {}
        otherdict = cfg[66].split("&")
        if len(otherdict)>1 or otherdict[0]!='':
            print(otherdict)
            for item in otherdict:
                parts = item.split("~")
                if parts[1]=="f":
                    dtype=float
                elif parts[1]=="i":
                    dtype=int
                elif parts[1]=="s":
                    dtype=str
                parts = parts[0].split("|")
                otherargs[parts[0]] = dtype(parts[1])
            
        self.configure(noutput=noutput,flux=flux,startemp=startemp,starspec=starspec,
                       gascon=gascon,pressure=pressure,pressurebroaden=pressurebroaden,
                       vtype=vtype,rotationperiod=rotationperiod,synchronous=synchronous,
                       substellarlon=substellarlon,restartfile=restartfile,gravity=gravity,
                       radius=radius,eccentricity=eccentricity,obliquity=obliquity,
                       lonvernaleq=lonvernaleq,fixedorbit=fixedorbit,orography=orography,
                       seaice=seaice,co2weathering=co2weathering,evolveco2=evolveco2,
                       physicsfilter=physicsfilter,filterkappa=filterkappa,
                       filterpower=filterpower,filterLHN0=filterLHN0,diffusionwaven=diffusionwaven,
                       qdiffusion=qdiffusion,tdiffusion=tdiffusion,zdiffusion=zdiffusion,
                       ddiffusion=ddiffusion,diffusionpower=diffusionpower,
                       erosionsupplylimit=erosionsupplylimit,outgassing=outgassing,
                       snowicealbedo=snowicealbedo,twobandalbedo=twobandalbedo,maxsnow=maxsnow,
                       soilalbedo=soilalbedo,oceanalbedo=oceanalbedo,oceanzenith=oceanzenith,
                       wetsoil=wetsoil,soilwatercap=soilwatercap,aquaplanet=aquaplanet,
                       desertplanet=desertplanet,soilsaturation=soilsaturation,
                       drycore=drycore,ozone=ozone,cpsoil=cpsoil,soildepth=soildepth,
                       mldepth=mldepth,writefrequency=writefrequency,modeltop=modeltop,
                       stratosphere=stratosphere,tropopause=tropopause,timestep=timestep,
                       runscript=runscript,columnmode=columnmode,highcadence=highcadence,
                       snapshots=snapshots,resources=resources,landmap=landmap,stormclim=stormclim,
                       nstorms=nstorms,stormcapture=stormcapture,topomap=topomap,
                       otherargs=otherargs)       
    
    def modify(self,**kwargs):
        setgas=False
        setgascon=False
        changeatmo=False
        changeland=False
        oldpressure = 0.0
        for gas in self.pgases:
            oldpressure += self.pgases[gas]
        for key,value in kwargs.items():
            if key=="noutput":
                self._edit_namelist("plasim_namelist","NOUTPUT",str(value*1))
                self.noutput = value
            if key=="flux":
                self._edit_namelist("planet_namelist","GSOL0",str(value))
                self.flux = value
            if key=="startemp":
                startemp=value
                if startemp:
                    self._edit_namelist("radmod_namelist","NSTARTEMP","1")
                    self._edit_namelist("radmod_namelist","STARBBTEMP",str(startemp))
                self.startemp = startemp
            if key=="starspec":
                starspec=value
                if starspec:
                    self._edit_namelist("radmod_namelist","NSTARFILE","1")
                    self._edit_namelist("radmod_namelist","STARFILE","'%s'"%starspec)
                    self._edit_namelist("radmod_namelist","STARFILEHR","'%s_hr.dat'"%(starspec)[:-4])
                self.starspec = starspec
            if key=="pH2":
                setgas=True
                if type(value)!=type(None):
                    self.pgases["pH2"]=value
                else:
                    self.pgases["pH2"]=0.0
            if key=="pHe":
                setgas=True
                if type(value)!=type(None):
                    self.pgases["pHe"]=value
                else:
                    self.pgases["pHe"]=0.0
            if key=="pN2":
                setgas=True
                if type(value)!=type(None):
                    self.pgases["pN2"]=value
                else:
                    self.pgases["pN2"]=0.0
            if key=="pO2":
                setgas=True
                if type(value)!=type(None):
                    self.pgases["pO2"]=value
                else:
                    self.pgases["pO2"]=0.0
            if key=="pCO2":
                setgas=True
                if type(value)!=type(None):
                    self.pgases["pCO2"]=value
                else:
                    self.pgases["pCO2"]=0.0
            if key=="pAr":
                setgas=True
                if type(value)!=type(None):
                    self.pgases["pAr"]=value
                else:
                    self.pgases["pAr"]=0.0
            if key=="pNe":
                setgas=True
                if type(value)!=type(None):
                    self.pgases["pNe"]=value
                else:
                    self.pgases["pNe"]=0.0
            if key=="pKr":
                setgas=True
                if type(value)!=type(None):
                    self.pgases["pKr"]=value
                else:
                    self.pgases["pKr"]=0.0
            if key=="pH20":
                setgas=True
                if type(value)!=type(None):
                    self.pgases["pH2O"]=value
                else:
                    self.pgases["pH2O"]=0.0
            if key=="pressure":
                pressure=value
                setpressure=True
            if key=="gascon":
                setgascon=True
                gascon=value
            if key=="pressurebroaden":
                self.pressurebroaden=value
                self._edit_namelist("radmod_namelist","NPBROADEN",str(self.pressurebroaden*1))
            if key=="vtype":
                self.vtype=value
                self._edit_namelist("plasim_namelist","NEQSIG",str(self.vtype))
            if key=="rotationperiod":
                self.rotationperiod=value
                if self.rotationperiod!=1.0:
                    self._edit_namelist("planet_namelist","ROTSPD",
                                        str(1.0/float(self.rotationperiod)))
                    self._edit_namelist("plasim_namelist","N_DAYS_PER_YEAR",
                                        str(max(int(360.0/float(self.rotationperiod)/12+0.5),1)*12))
            if key=="synchronous":
                self.synchronous=value
                self._edit_namelist("radmod_namelist","NFIXED",str(self.synchronous*1))
            if key=="substellarlon":
                self.substellarlon=value
                self._edit_namelist("radmod_namelist","FIXEDLON",str(self.substellarlon))
            if key=="restartfile":
                self.restartfile=value
                if self.restartfile:
                    os.system("cp %s %s/plasim_restart"%(restartfile,self.workdir))
                else:
                    os.system("rm %s/plasim_restart"%self.workdir)
            if key=="gravity":
                self.gravity=value
                self._edit_namelist("planet_namelist","GA",str(self.gravity))
                self._edit_postnamelist("example.nl","gravity",str(self.gravity))
                self._edit_postnamelist("snapshot.nl","gravity",str(self.gravity))
            if key=="radius":
                self.radius=value
                self._edit_namelist("planet_namelist","PLARAD",str(self.radius*6371220.0))
                self._edit_postnamelist("example.nl","radius",str(self.radius*6371220.0))
                self._edit_postnamelist("snapshot.nl","radius",str(self.radius*6371220.0))
            if key=="eccentricity":
                self.eccentricity=value
                self._edit_namelist("planet_namelist","ECCEN",str(self.eccentricity))
            if key=="obliquity":
                self.obliquity=value
                self._edit_namelist("planet_namelist","OBLIQ",str(self.obliquity))
            if key=="lonvernaleq":
                self.lonvernaleq=value
                self._edit_namelist("planet_namelist","MVELP",str(self.lonvernaleq))
            if key=="fixedorbit":
                self.fixedorbit=value
                self._edit_namelist("planet_namelist","NFIXORB",str(self.fixedorbit*1))
                
            if key=="orography":
                self.orography=value
                self._edit_namelist("landmod_namelist","OROSCALE",str(self.orography))
                self._edit_namelist("glacier_namelist","NGLACIER",str((self.orography!=1)*1))
            if key=="seaice":
                self.seaice=value
                self._edit_namelist("radmod_namelist","NRADICE",str(self.seaice*1))
            if key=="co2weathering":
                self.co2weathering=value
                self._edit_namelist("carbonmod_namelist","NCARBON",str(self.co2weathering*1))
            if key=="evolveco2":
                self.evolveco2=value
                self._edit_namelist("carbonmod_namelist","NCO2EVOLVE",str(self.evolveco2*1))
            if key=="erosionsupplylimit":
                self.erosionsupplylimit=value
                flag = bool(self.erosionsupplylimit)*1
                self._edit_namelist("carbonmod_namelist","NSUPPLY",str(flag))
                self._edit_namelist("carbonmod_namelist","WMAX",
                                        str(self.erosionsupplylimit*flag+1.0*(1-flag)))
            if key=="outgassing":
                self.outgassing=value
                self._edit_namelist("carbonmod_namelist","VOLCANCO2",str(self.outgassing))
                
            if key=="physicsfilter":
                self.physicsfilter=value
                if self.physicsfilter:
                    vals = self.physicsfilter.split("|")
                    if "gp" in vals:
                        self._edit_namelist("plasim_namelist","NGPTFILTER","1")
                    if "sp" in vals:
                        self._edit_namelist("plasim_namelist","NSPVFILTER","1")
                    if "none" in vals:
                        self._edit_namelist("plasim_namelist","NFILTER","0")
                    if "cesaro" in vals:
                        self._edit_namelist("plasim_namelist","NFILTER","1")
                    if "exp" in vals:
                        self._edit_namelist("plasim_namelist","NFILTER","2")
                    if "lh" in vals:
                        self._edit_namelist("plasim_namelist","NFILTER","3")
                else:
                    self._edit_namelist("plasim_namelist","NGPTFILTER","0")
                    self._edit_namelist("plasim_namelist","NSPVFILTER","0")
            if key=="filterkappa":
                self.filterkappa=value
                self._edit_namelist("plasim_namelist","FILTERKAPPA",str(self.filterkappa))
            if key=="filterpower":
                self.filterpower=value
                self._edit_namelist("plasim_namelist","NFILTEREXP",str(self.filterpower))
            if key=="filterLHN0":
                self.filterLHN0=value
                self._edit_namelist("plasim_namelist","LANDHOSKN0",str(self.filterLHN0))
            if key=="diffusionwaven":
                self.diffusionwaven=value
                if value:
                    self._edit_namelist("plasim_namelist","NHDIFF",str(self.diffusionwaven))
                else:
                    self._rm_namelist_param("plasim_namelist","NHDIFF")
            if key=="qdiffusion":
                self.qdiffusion=value
                if value:
                    self._edit_namelist("plasim_namelist","TDISSQ","%d*%f"%(self.layers,
                                                                        self.qdiffusion))
                else:
                    self._rm_namelist_param("plasim_namelist","TDISSQ")
            if key=="tdiffusion":
                self.tdiffusion=value
                if value:
                    self._edit_namelist("plasim_namelist","TDISST","%d*%f"%(self.layers,
                                                                        self.tdiffusion))
                else:
                    self._rm_namelist_param("plasim_namelist","TDISST")
            if key=="zdiffusion":
                self.zdiffusion=value
                if value:
                    self._edit_namelist("plasim_namelist","TDISSZ","%d*%f"%(self.layers,
                                                                        self.zdiffusion))
                else:
                    self._rm_namelist_param("plasim_namelist","TDISSZ")
            if key=="ddiffusion":
                self.ddiffusion=value
                if value:
                    self._edit_namelist("plasim_namelist","TDISSD","%d*%f"%(self.layers,
                                                                        self.ddiffusion))
                else:
                    self._rm_namelist_param("plasim_namelist","TDISSD")
            if key=="diffusionpower":
                self.diffusionpower=value
                if value:
                    self._edit_namelist("plasim_namelist","NDEL","%d*%d"%(self.layers,
                                                                      self.diffusionpower))
                else:
                    self._rm_namelist_param("plasim_namelist","NDEL")
            if key=="snowicealbedo":
                self.snowicealbedo=value
                if type(self.snowicealbedo)!=type(None):
                    alb = str(self.snowicealbedo)
                    self._edit_namelist("seamod_namelist","ALBICE",alb)
                    self._edit_namelist("seamod_namelist","DICEALBMN","%s,%s"%(alb,alb))
                    self._edit_namelist("seamod_namelist","DICEALBMX","%s,%s"%(alb,alb))
                    self._edit_namelist("landmod_namelist","DSNOWALBMN","%s,%s"%(alb,alb))
                    self._edit_namelist("landmod_namelist","DSNOWALBMN","%s,%s"%(alb,alb))
                    self._edit_namelist("landmod_namelist","DGLACALBMN","%s,%s"%(alb,alb))
                    self._edit_namelist("landmod_namelist","DSNOWALB","%s,%s"%(alb,alb))
                else:
                    self._rm_namelist_param("seamod_namelist","ALBICE")
                    self._rm_namelist_param("seamod_namelist","DICEALBMN")
                    self._rm_namelist_param("seamod_namelist","DICEALBMX")
                    self._rm_namelist_param("landmod_namelist","DSNOWALBMN")
                    self._rm_namelist_param("landmod_namelist","DSNOWALBMN")
                    self._rm_namelist_param("landmod_namelist","DGLACALBMN")
                    self._rm_namelist_param("landmod_namelist","DSNOWALB")
            if key=="twobandalbedo":
                self.twobandalbedo=value
                self._edit_namelist("radmod_namelist","NSIMPLEALBEDO",
                                    str((not self.twobandalbedo)*1))
            if key=="maxsnow":
                self.maxsnow=value
                if maxsnow:
                    self._edit_namelist("landmod_namelist","DSMAX",str(self.maxsnow))
                else:
                    self._rm_namelist_param("landmod_namelist","DSMAX")
            if key=="soilalbedo":
                self.soilalbedo=value
                if type(self.soilalbedo)!=type(None):
                    alb = str(self.soilalbedo)
                    os.system("rm %s/*0174.sra"%self.workdir)
                    self._edit_namelist("landmod_namelist","ALBLAND",alb)
                    self._edit_namelist("landmod_namelist","DGROUNDALB","%s,%s"%(alb,alb))
                else:
                    self._rm_namelist_param("landmod_namelist","ALBLAND")
                    self._rm_namelist_param("landmod_namelist","DGROUNDALB")
            
            if key=="oceanalbedo":
                self.oceanalbedo=value
                if type(self.oceanalbedo)!=type(None):
                    alb = str(self.oceanalbedo)
                    self._edit_namelist("seamod_namelist","ALBSEA",alb)
                    self._edit_namelist("seamod_namelist","DOCEANALB","%s,%s"%(alb,alb))
                else:
                    self._rm_namelist_param("seamod_namelist","ALBSEA")
                    self._rm_namelist_param("seamod_namelist","DOCEANALB")
                
            if key=="oceanzenith":
                self.oceanzenith=value
                if self.oceanzenith=="lambertian" or self.oceanzenith=="Lambertian" or self.oceanzenith=="uniform":
                    self._edit_namelist("radmod_namelist","NECHAM","0")
                    self._edit_namelist("radmod_namelist","NECHAM6","0")
                if self.oceanzenith=="default" or self.oceanzenith=="plasim" or self.oceanzenith=="ECHAM-3":
                    self._edit_namelist("radmod_namelist","NECHAM","1")
                    self._edit_namelist("radmod_namelist","NECHAM6","0")
                if self.oceanzenith=="ECHAM-6":
                    self._edit_namelist("radmod_namelist","NECHAM","0")
                    self._edit_namelist("radmod_namelist","NECHAM6","1")
                
            if key=="wetsoil":
                self.wetsoil=value
                self._edit_namelist("landmod_namelist","NWETSOIL",str(self.wetsoil*1))
            if key=="soilwatercap":
                self.soilwatercap=value
                if self.soilwatercap:
                    self._edit_namelist("landmod_namelist","WSMAX",str(self.soilwatercap))
                    os.system("rm %s/*0229.sra"%self.workdir)
                else:
                    self._rm_namelist_param("landmod_namelist","WSMAX")
            if key=="soilsaturation":
                self.soilsaturation=value
                if type(self.soilsaturation)!=type(None):
                    self._edit_namelist("landmod_namelist","NWATCINI","1")
                    self._edit_namelist("landmod_namelist","DWATCINI",str(self.soilsaturation))
                    os.system("rm %s/*0229.sra"%self.workdir)
                else:
                    self._edit_namelist("landmod_namelist","NWATCINI","0")
                    self._rm_namelist_param("landmod_namelist","DWATCINI")
            if key=="desertplanet":
                self.desertplanet=value
                if self.desertplanet:
                    self._edit_namelist("plasim_namelist","NDESERT","1")
                    self._edit_namelist("landmod_namelist","NWATCINI","1")
                    self._edit_namelist("landmod_namelist","DWATCINI","0.0")
                    os.system("rm %s/*.sra"%self.workdir)
                else:
                    self._edit_namelist("landmod_namelist","NDESERT","0")
                    self._rm_namelist_param("landmod_namelist","NWATCINI")
                    self._rm_namelist_param("landmod_namelist","DWATCINI")
            if key=="aquaplanet":
                self.aquaplanet=value
                if self.aquaplanet:
                    self._edit_namelist("plasim_namelist","NAQUA","1")
                    os.system("rm %s/*.sra"%self.workdir)
                else:
                    self._edit_namelist("plasim_namelist","NAQUA","0")

            if key=="drycore":
                self.drycore=value
                self._edit_namelist("fluxmod_namelist","NEVAP",str((not self.drycore)*1))
        
            if key=="columnmode":
                self.columnmode=value
                if self.columnmode:
                    parts = self.columnmode.split("|")
                    if "static" in parts:
                        self._edit_namelist("plasim_namelist","NADV","0")
                    else:
                        self._edit_namelist("plasim_namelist","NADV","1")
                    if "clear" in parts:
                        self._edit_namelist("radmod_namelist","NCLOUDS","0")
                        self._edit_namelist("radmod_namelist","ACLLWR","0.0")
                    else:
                        self._edit_namelist("radmod_namelist","NCLOUDS","1")
                        self._rm_namelist_param("radmod_namelist","ACLLWR")
                else:
                    self._rm_namelist_param("plasim_namelist","NADV")
                    self._rm_namelist_param("radmod_namelist","NCLOUDS")
                    self._rm_namelist_param("radmod_namelist","ACLLWR")
                    
            if key=="ozone":
                self.ozone=value
                self._edit_namelist("radmod_namelist","NO3",str(self.ozone*1))
                
            if key=="cpsoil":
                self.cpsoil=value
                if self.cpsoil:
                    self._edit_namelist("landmod_namelist","SOILCAP",str(self.cpsoil))
                else:
                    self._rm_namelist_param("landmod_namelist","SOILCAP")
                
            if key=="soildepth":
                self.soildepth=value
                self.dzsoils = np.array([0.4, 0.8, 1.6, 3.2, 6.4])*soildepth
                self._edit_namelist("landmod_namelist","DSOILZ",",".join(self.dzsoils.astype(str)))
                
            if key=="mldepth":
                self.mldepth=value
                self._edit_namelist("oceanmod_namelist","MLDEPTH",str(self.mldepth))
                
            if key=="writefrequency":
                self.writefrequency=value
                if self.writefrequency:
                    self._edit_namelist("plasim_namelist","NWPD",str(self.writefrequency))
                else:
                    self._rm_namelist_param("plasim_namelist","NWPD")
                
            if key=="modeltop":
                changeatmo=True
                self.modeltop=value
            if key=="stratosphere":
                changeatmo=True
                self.stratosphere=value
            if key=="tropopause":
                changeatmo=True
                self.tropopause=value
             
            if key=="timestep":
                self.timestep=value
                self._edit_namelist("plasim_namelist","MPSTEP",str(self.timestep))
                self._edit_namelist("plasim_namelist","NSTPW",str(int(7200//int(self.timestep))))
            
            if key=="runscript":
                self.runscript=value
                
            if key=="highcadence":
                self.highcadence=value
                self._edit_namelist("plasim_namelist","NHCADENCE",str(self.highcadence["toggle"]))
                self._edit_namelist("plasim_namelist","HCSTARTSTEP",str(self.highcadence["start"]))
                self._edit_namelist("plasim_namelist","HCENDSTEP",str(self.highcadence["end"]))
                self._edit_namelist("plasim_namelist","HCINTERVAL",str(self.highcadence["interval"]))
                
            if key=="snapshots":
                self.snapshots=value
                if self.snapshots:
                    self._edit_namelist("plasim_namelist","NSNAPSHOT","1")
                    self._edit_namelist("plasim_namelist","NSTPS",str(self.snapshots))
                else:
                    self._rm_namelist_param("plasim_namelist","NSTPS")
                    self._edit_namelist("plasim_namelist","NSNAPSHOT","0")
                    
            if key=="resources":
                self.resources=value
                if len(self.resources)>0:
                    for res in self.resources:
                        os.system("cp %s %s/"%(res,self.workdir))
                
            if key=="landmap":
                self.landmap=value
                changeland=True
            if key=="topomap":
                self.topomap=value
                changeland=True
                
            if key=="stormclim":
                self.stormclim=value
                if self.stormclim:
                    self._edit_namelist("hurricane_namelist","NSTORMDIAG","1")
                    self._add_postcodes("example.nl",[322,323,324,325,326,327,328,329])
                    self._add_postcodes("snapshot.nl",[322,323,324,325,326,327,328,329])
                else:
                    self._edit_namelist("hurricane_namelist","NSTORMDIAG","0")
                    self._rm_postcodes("example.nl",[322,323,324,325,326,327,328,329])
                    self._rm_postcodes("snapshot.nl",[322,323,324,325,326,327,328,329])
            if key=="nstorms":
                self.nstorms=value
                self._edit_namelist("hurricane_namelist","NSTORMS",str(self.int(nstorms)))
            if key=="stormcapture":
                self.stormcapture=value
                if self.stormcapture["toggle"]:
                    self._edit_namelist("hurricane_namelist","HC_CAPTURE","1")
                    for param in self.stormcapture:
                        if param!="toggle":
                            self._edit_namelist("hurricane_namelist",param,
                                                str(self.stormcapture[param]))
                else:
                    self._edit_namelsit("hurricane_namelist","HC_CAPTURE","0")
                
            if key=="otherargs":
                otherargs=value
                if len(otherargs)>0:
                    for key in otherargs:
                        parts = key.split("<")
                        destination=parts[0].split("@")
                        value=parts[1]
                        field=destination[0]
                        namelist=destination[1]
                        self._edit_namelist(namelist,field,value)
                        self.otherargs.append(key)
                
        if setgas:
            
            if setpressure:
                newpressure=pressure
            
            pressure=0.0
            for gas in self.pgases:
                pressure+=self.pgases[gas]
            
            pscalef = pressure/oldpressure
            
            gasesvx = {}
            for gas in self.pgases:
                gasesvx[gas[1:]] = self.pgases[gas]/pressure
            self.mmw = 0
            for gas in gasesvx:
                self.mmw += gasesvx[gas]*smws['m'+x]
            self.gascon = 8314.46261815324 / self.mmw
            
            if setgascon:
                self.gascon=gascon
                self.mmw = 8314.46261815324/self.gascon
            
            print('Mean Molecular Weight set to %1.4f g/mol'%self.mmw)
            if setpressure:
                self.pressure=newpressure
            else:
                self.pressure *= pscalef
            print("Surface Pressure set to %1.6f bars"%self.pressure)
            self._edit_namelist("plasim_namelist","PSURF",str(self.pressure))
            self._edit_namelist("planet_namelist","GASCON",str(self.gascon))
            
                
            if 'pCO2' in self.pgases:
                self.CO2ppmv = self.pgases['pCO2']/self.pressure * 1.0e6 #ppmv
                self._edit_namelist("radmod_namelist","CO2",self.CO2ppmv)
        
        else:
            if setpressure:
                self.pressure=pressure
                self._edit_namelist("plasim_namelist","PSURF",str(self.pressure))
            
        if changeatmo:
               
            if self.stratosphere:
                self._edit_namelist("plasim_namelist","NEQSIG","5")
                if self.modeltop:
                    self._edit_namelist("plasim_namelist","PTOP2",str(self.modeltop*100.0)) #convert hPa->Pa
                if self.tropopause:
                    self._edit_namelist("plasim_namelist","PTOP",str(self.tropopause*100.0))
            else:
                if self.modeltop:
                    self._edit_namelist("plasim_namelist","PTOP",str(self.modeltop*100.0)) #convert hPa->Pa
                    
        if changeland:
            if self.landmap or self.topomap:
                os.system("rm %s/*.sra"%self.workdir)
            if self.landmap:
                os.system("cp %s %s/N%03d_surf_0172.sra"%(self.landmap,self.workdir,self.nlats))
            if self.topomap:
                os.system("cp %s %s/N%03d_surf_0129.sra"%(self.topomap,self.workdir,self.nlats))
    
    def save(self,filename):
        np.save(filename,self,allow_pickle=True)
        
    def exportcfg(self,filename=None):
        '''export model configuration to a text file that can be used as configuration input'''
        if not filename:
            filename = self.modelname+".cfg"
        cfg = []
        
        cfg.append(str(self.noutput*1))
        cfg.append(str(self.flux))
        cfg.append(str(self.startemp))
        cfg.append(str(self.starspec))#=noneparse(cfg[3],str)
        gases = []
        for gas in self.pgases:
            gases.append(gas+"|"+str(self.pgases[gas]))
        gases = "&".join(gases)
        cfg.append(gases)# = cfg[4].split("&")
        cfg.append(str(self.gascon))# = float(cfg[5])
        cfg.append(str(self.pressure))# = float(cfg[6])
        cfg.append(str(self.pressurebroaden*1))# = bool(int(cfg[7]))
        cfg.append(str(self.vtype))# = int(cfg[8])
        cfg.append(str(self.rotationperiod))# = float(cfg[9])
        cfg.append(str(self.synchronous*1))# = bool(int(cfg[10]))
        cfg.append(str(self.substellarlon))# = float(cfg[11])
        cfg.append(str(self.restartfile))# = noneparse(cfg[12],str)
        cfg.append(str(self.gravity))# = float(cfg[13])
        cfg.append(str(self.radius))# = float(cfg[14])
        cfg.append(str(self.eccentricity))# = noneparse(cfg[15],float)
        cfg.append(str(self.obliquity))# = noneparse(cfg[16],float)
        cfg.append(str(self.lonvernaleq))# = noneparse(cfg[17],float)
        cfg.append(str(self.fixedorbit*1))# = bool(int(cfg[18]))
        cfg.append(str(self.orography))# = noneparse(cfg[19],float)
        cfg.append(str(self.seaice*1))# = bool(int(cfg[20]))
        cfg.append(str(self.co2weathering*1))# = bool(int(cfg[21]))
        cfg.append(str(self.evolveco2*1))# = bool(int(cfg[22]))
        cfg.append(str(self.physicsfilter))# = noneparse(cfg[23],str)
        cfg.append(str(self.filterkappa))# = float(cfg[24])
        cfg.append(str(self.filterpower))# = int(cfg[25])
        cfg.append(str(self.filterLHN0))# = float(cfg[26])
        cfg.append(str(self.diffusionwaven))# = noneparse(cfg[27],int)
        cfg.append(str(self.qdiffusion))# = noneparse(cfg[28],float)
        cfg.append(str(self.tdiffusion))# = noneparse(cfg[29],float)
        cfg.append(str(self.zdiffusion))# = noneparse(cfg[30],float)
        cfg.append(str(self.ddiffusion))# = noneparse(cfg[31],float)
        cfg.append(str(self.diffusionpower))# = noneparse(cfg[32],int)
        cfg.append(str(self.erosionsupplylimit))# = noneparse(cfg[33],float)
        cfg.append(str(self.outgassing))# = float(cfg[34])
        cfg.append(str(self.snowicealbedo))# = noneparse(cfg[35],float)
        cfg.append(str(self.twobandalbedo*1))# = bool(int(cfg[36]))
        cfg.append(str(self.maxsnow))# = noneparse(cfg[37],float)
        cfg.append(str(self.soilalbedo))# = noneparse(cfg[38],float)
        cfg.append(str(self.oceanalbedo))# = noneparse(cfg[39],float)
        cfg.append(str(self.oceanzenith))# = cfg[40]
        cfg.append(str(self.wetsoil*1))# = bool(int(cfg[41]))
        cfg.append(str(self.soilwatercap))# = noneparse(cfg[42],float)
        cfg.append(str(self.aquaplanet*1))# = bool(int(cfg[43]))
        cfg.append(str(self.desertplanet*1))# = bool(int(cfg[44]))
        cfg.append(str(self.soilsaturation))# = noneparse(cfg[45],float)
        cfg.append(str(self.drycore*1))# = bool(int(cfg[46]))
        cfg.append(str(self.ozone*1))# = bool(int(cfg[47]))
        cfg.append(str(self.cpsoil))# = noneparse(cfg[48],float)
        cfg.append(str(self.soildepth))# = float(cfg[49])
        cfg.append(str(self.mldepth))# = float(cfg[50])
        cfg.append(str(self.writefrequency))# = noneparse(cfg[51],int)
        cfg.append(str(self.modeltop))# = noneparse(cfg[52],float)
        cfg.append(str(self.stratosphere*1))# = bool(int(cfg[53]))
        cfg.append(str(self.tropopause))# = noneparse(cfg[54],float)
        cfg.append(str(self.timestep))# = float(cfg[55])
        cfg.append(str(self.runscript))# = noneparse(cfg[56],str)
        cfg.append(str(self.columnmode))# = noneparse(cfg[57],str)
        hcdict = []
        for hc in self.highcadence:
            hcdict.append(hc+"|"+str(self.highcadence[hc]))
        hcdict = "&".join(hcdict)
        cfg.append(hcdict)
        cfg.append(str(self.snapshots))# = noneparse(cfg[59],int)
        cfg.append("&".join(self.resources))
        cfg.append(str(self.landmap))# = noneparse(cfg[61],str)
        cfg.append(str(self.topomap))# = noneparse(cfg[62],str)
        cfg.append(str(self.stormclim*1))# = bool(int(cfg[63]))
        cfg.append(str(self.nstorms))# = int(cfg[64])
        stormdict = []
        for arg in self.stormcapture:
            stormdict.append(arg+"|"+str(self.stormcapture[arg]))
        stormdict = "&".join(stormdict)
        cfg.append(stormdict)
        otherdict = []
        for arg in self.otherargs:
            item = arg+"|"+str(self.otherargs[arg])
            if type(self.otherargs[arg])==int:
                item+="~i"
            if type(self.otherargs[arg])==str:
                item+="~s"
            if type(self.otherargs[arg])==float:
                item+="~f"
            otherdict.append(item)
        cfg.append("&".join(otherdict))
        
        with open(filename,"w") as cfgf:
            cfgf.write("\n".join(cfg))
        
    def _rm_namelist_param(self,namelist,arg,val):
        '''Remove an argument from a namelist'''
        
        f=open(self.workdir+"/"+namelist,"r")
        fnl=f.read().split('\n')
        f.close()
        found=False
        fnl1=fnl[1].split(' ')
        if '=' in fnl1:
            mode='EQ'
        else:
            mode='CM'
        #print fnl1
        item = fnl1[-1]
        if item=='':
            item = fnl1[-2]
        if item.strip()[-1]!=",":
            mode='EQ'
        
        for l in range(1,len(fnl)-2):
            fnl[l]=fnl[l].split(' ')
            if arg in fnl[l]:
                found=True
            elif (arg+'=') in fnl[l]:
                found=True
            if found:
                fnl.pop(l)
                break
                
        f=open(self.workdir+"/"+namelist,"w")
        f.write('\n'.join(fnl))
        f.close()
              
    def _edit_namelist(self,namelist,arg,val):
        '''Either edit or add argument/value pair to a namelist'''
        
        f=open(self.workdir+"/"+namelist,"r")
        fnl=f.read().split('\n')
        f.close()
        found=False
        fnl1=fnl[1].split(' ')
        if '=' in fnl1:
            mode='EQ'
        else:
            mode='CM'
        #print fnl1
        item = fnl1[-1]
        if item=='':
            item = fnl1[-2]
        if item.strip()[-1]!=",":
            mode='EQ'
        
        for l in range(1,len(fnl)-2):
            fnl[l]=fnl[l].split(' ')
            if arg in fnl[l]:
                fnl[l]=['',arg,'','=','',str(val),'']
                found=True
            elif (arg+'=') in fnl[l]:
                tag = ','
                item = fnl[l][-1]
                if item=='':
                    item = fnl[l][-2]
                if item.strip()[-1]!=',':
                    tag = ''
                fnl[l]=['',arg+'=','',str(val),'',tag]
                found=True
            fnl[l]=' '.join(fnl[l])
        if not found:
            if mode=='EQ':
                fnl.insert(-3,' '+arg+' = '+val+' ')
            else:
                fnl.insert(-3,' '+arg+'= '+val+' ,')
            
        f=open(self.workdir+"/"+namelist,"w")
        f.write('\n'.join(fnl))
        f.close()
        
    def _edit_postnamelist(self,namelist,arg,val):
        '''Edit postprocessing namelist'''
        
        with open(self.workdir+"/"+namelist,"r") as f:
            pnl = f.read().split('\n')
            
        flag=False
        pnl = [y for y in pnl if y!='']
        for n in range(len(pnl)):
            if pnl[n].split('=')[0].strip()==arg:
                pnl[n]=arg+"="+val
                flag=True
                break
        if not flag:
            pnl.append(arg+'='+val)
        pnl.append('')
        
        with open(self.workdir+"/"+namelist,"w") as f:
            f.write('\n'.join(pnl))
            
    def _add_postcodes(self,namelist,newcodes):
        '''Add postprocessor codes to postprocessor namelist'''
        
        with open(self.workdir+"/"+namelist,"r") as f:
            pnl = f.read().split('\n')
        pnl = [y for y in pnl if y!='']
        for n in range(len(pnl)):
            if pnl[n].split('=')[0].strip()=="code":
                codes = pnl[n].split('=')[1].strip().split(',')
                lineno=n
                break
        ncodes = [int(n) for n in codes]
        for n in newcodes:
            if n not in ncodes:
                ncodes.append(n)
        pnl[lineno]+=','+','.join([str(n) for n in ncodes])
        #print "Writing to %s/%s: \n"%(home,filename)+'\n'.join(pnl)+"\n"
        with open(self.workdir+"/"+namelist,"w") as f:
            f.write('\n'.join(pnl)+"\n")

    def _rm_postcodes(self,namelist,rmcodes):
        '''Add postprocessor codes to postprocessor namelist'''
        
        with open(self.workdir+"/"+namelist,"r") as f:
            pnl = f.read().split('\n')
        pnl = [y for y in pnl if y!='']
        for n in range(len(pnl)):
            if pnl[n].split('=')[0].strip()=="code":
                codes = pnl[n].split('=')[1].strip().split(',')
                lineno=n
                break
        ncodes = [int(n) for n in codes]
        
        newcodes = []
        for n in ncodes:
            if n not in rmcodes:
                newcodes.append(n)
        pnl[lineno]+=','+','.join([str(n) for n in newcodes])
        #print "Writing to %s/%s: \n"%(home,filename)+'\n'.join(pnl)+"\n"
        with open(self.workdir+"/"+namelist,"w") as f:
            f.write('\n'.join(pnl)+"\n")


class TLaquaplanet(Model):
    def configure(self,timestep=30.0,snapshots=720,eccentricity=0.0,ozone=False,
                  obliquity=0.0,physicsfilter="gp|exp|sp",**kwargs):
        super(TLaquaplanet,self).configure(synchronous=True,fixedorbit=True,aquaplanet=True,
                          eccentricity=eccentricity,obliquity=obliquity,timestep=timestep,
                          snapshots=snapshots,physicsfilter=physicsfilter,ozone=ozone,
                          **kwargs)
        
class TLlandplanet(Model): #Default will be ZERO soil water; set soilsaturation if you want any
    def configure(self,timestep=30.0,snapshots=720,eccentricity=0.0,ozone=False,
                  obliquity=0.0,physicsfilter="gp|exp|sp",**kwargs):
        super(TLlandplanet,self).configure(synchronous=True,fixedorbit=True,desertplanet=True,
                          eccentricity=eccentricity,obliquity=obliquity,timestep=timestep,
                          snapshots=snapshots,physicsfilter=physicsfilter,ozone=ozone,
                          **kwargs)
        
class Earthlike(Model):
    def configure(self,timestep=45.0,snapshots=480,**kwargs):
        super(Earthlike,self).configure(vtype=4,modeltop=50.0,timestep=timestep,
                          snapshots=snapshots,**kwargs)

