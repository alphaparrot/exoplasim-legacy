import numpy as np
import netCDF4 as nc
import sys
import glob

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

if __name__=="__main__":
    
    prefix = sys.argv[1]
    
    start = int(sys.argv[2])
    try:
        end = int(sys.argv[3])    
    except:
        end = len(glob.glob(prefix+"*.nc"))
    
    keys = ["clt" ,
            "hfls",
            "hfns",
            "hfss",
            "nbr" ,
            "prw" ,
            "rls" ,
            "rlut",
            "rss" ,
            "rst" ,
            "rsut",
            "ssru",
            "stru",
            "sic",
            "sit",
            "ts"  ]
    
    nrecs = end - start + 1
    
    history = {}
    
    for k in keys:
        history[k] = np.zeros(nrecs)
    
    for n in range(start,end+1):
        ncd = nc.Dataset(prefix+".%04d.nc"%n,"r")
        lat = ncd.variables['lat'][:]
        lon = ncd.variables['lon'][:]
        for k in keys:
            try:
                varx = ncd.variables[k][:]
                nt = varx.shape[0]
                x = 0
                for t in range(0,nt):
                    x += spatialmath(lat,lon,varx[t,:,:])
                x /= float(nt)
                history[k][n-start] = x
            except:
                pass
        ncd.close()
        print(n)
        
    np.save(prefix+"_history.npy",history)
    
    
    