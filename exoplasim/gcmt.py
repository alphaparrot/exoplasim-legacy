import numpy as np
import netCDF4 as nc

class _Dataset:
    def __init__(self,filename):
        if filename[-3:]==".nc":
            self.body = nc.Dataset(filename,"r")
            self.variables=self.body.variables
        elif filename[-4:]==".npy":
            self.body = np.load(filename)
            self.variables = self.body.item()
        else:
            raise DatafileError("Unknown dataset format")
        
    def close(self):
        try:
            self.body.close()
        except:
            pass

class DimensionError(Exception):
    pass

class UnitError(Exception):
    pass

class DatafileError(Exception):
    pass

def parse(file,variable,lat=None,lon=None):
    """Retrieve a variable from a NetCDF file

    Parameters
    ----------
    file : str
        Path to a NetCDF file
    variable : str
        Name of the variable to extract
    lat,lon : str, optional
        If the latitude and longitude arrays have non-standard names, specify them here.

    Returns
    -------
    numpy.ndarray
        Requested output field

    """
    ncd=_Dataset(file)
    variable = ncd.variables[variable][:]
    
    if lat:
        lt = ncd.variables[lat][:]
    elif "lat" in ncd.variables:
        lt = ncd.variables['lat'][:]
    elif "lt" in ncd.variables:
        lt = ncd.variables["lt"][:]
    elif "lats" in ncd.variables:
        lt = ncd.variables['lats'][:]
    elif "latitude" in ncd.variables:
        lt = ncd.variables['latitude'][:]
    elif "latitudes" in ncd.variables:
        lt = ncd.variables['latitudes'][:]
    else:
        raise DatafileError("Unknown datafile format; unsure how to extract latitude")
    
    if lon:
        ln = ncd.variables[lon][:]
    elif "lat" in ncd.variables:
        ln = ncd.variables['lon'][:]
    elif "lt" in ncd.variables:
        ln = ncd.variables["ln"][:]
    elif "lats" in ncd.variables:
        ln = ncd.variables['lons'][:]
    elif "latitude" in ncd.variables:
        ln = ncd.variables['longitude'][:]
    elif "latitudes" in ncd.variables:
        ln = ncd.variables['longitudes'][:]
    else:
        raise DatafileError("Unknown datafile format; unsure how to extract longitude")
    
    ncd.close()
    
    return ln,lt,variable

def make2d(variable,lat=None,lon=None,time=None,lev=None,ignoreNaNs=True):
    """Compress a variable in two dimensions by slicing or averaging.

    Parameters
    ----------
    variable : numpy.ndarray
        The variable to operate on
    lat,lon,lev : int, str, optional
        Either an index on which to slice, or either of "sum" or "mean", indicating what
        should be done along that axis.
    time : int, optional
        The time index on which to slice. If unspecified, a time average will be returned.
    ignoreNaNs : bool, optional
        If set, will use NaN-safe numpy operators.
        
    Returns
    -------
    numpy.ndarray
        A 2-D array
    """
    if ignoreNaNs:
        sumop = np.nansum
        meanop = np.nanmean
    else:
        sumop = np.sum
        meanop = np.mean
    if len(variable.shape)==2:
        return variable
    if time:
        try:
            variable=variable[time,:]
        except:
            raise UnitError("You have probably passed a float time to a variable with no "+
                            "information about what that means. You should pass an integer "+
                            "time index instead")
    elif time==None and len(variable.shape)>2:
        variable=meanop(variable,axis=0)
    elif time==0:
        variable=variable[time,:]
    if len(variable.shape)>2:
        if lev!=None:
            if type(lev)==int:
                variable=variable[lev,:]
            elif lev=="sum":
                variable=sumop(variable,axis=0)
            elif lev=="mean":
                variable=meanop(variable,axis=0)
            else:
                raise UnitError("Unknown level specification")
        elif lat!=None and lon==None:
            if type(lat)==int:
                variable=variable[:,lat,:]
            elif lat=="sum":
                variable=sumop(variable,axis=1)
            elif lat=="mean":
                variable=meanop(variable,axis=1)
            else:
                raise UnitError("Unknown latitude specification")
        elif lon!=None and lat==None:
            if type(lon)==int:
                variable=variable[:,:,lon]
            elif lon=="sum":
                variable=sumop(variable,axis=2)
            elif lon=="mean":
                variable=meanop(variable,axis=2)
            else:
                raise UnitError("Unknown longitude specification")
        else:
            raise DimensionError("Inappropriate or insufficient dimensional constraints")
    
    return variable
    

def spatialmath(variable,lat=None,lon=None,file=None,mean=True,time=None,
               ignoreNaNs=True,lev=None,radius=6.371e6):
    """Compute spatial means or sums of data

    Parameters
    ----------
    variable : str, numpy.ndarray
        The variable to operate on. Can either be a data array, or the name of a variable. If the latter, file must be specified.
    lat,lon : numpy.ndarray, optional
        Latitude and longitude arrays. If file is provided and lat and lon are not, they will be
        extracted from the file.
    file : str, optional
        Path to a NetCDF output file to open and extract data from.
    mean : bool, optional
        If True, compute a global mean. If False, compute a global sum.
    time : int, optional
        The time index on which to slice. If unspecified, a time average will be returned.
    ignoreNaNs : bool, optional
        If True, use NaN-safe numpy operators.
    lev : int, optional
        If set, slice a 3D spatial array at the specified level.
    radius : float, optional
        Radius of the planet in meters. Only used if mean=False.
        
    Returns
    -------
    float

    """
    
    if ignoreNaNs:
        sumop = np.nansum
        meanop = np.nanmean
    else:
        sumop = np.sum
        meanop = np.mean
        
    if file:
        ln,lt,variable = parse(file,variable,lat=lat,lon=lon)
        
    else:
        if type(lat)==type(None) or type(lon)==type(None):
            raise DimensionError("Need to provide latitude and longitude data")
        ln=lon
        lt=lat
    variable = make2d(variable,time=time,lev=lev,ignoreNaNs=ignoreNaNs)
    
    lt1 = np.zeros(len(lt)+1)
    lt1[0] = 90
    for n in range(0,len(lt)-1):
        lt1[n+1] = 0.5*(lt[n]+lt[n+1])
    lt1[-1] = -90
    dln = np.diff(ln)[0]
    ln1 = np.zeros(len(ln)+1)
    ln1[0] = -dln
    for n in range(0,len(ln)-1):
        ln1[n+1] = 0.5*(ln[n]+ln[n+1])
    ln1[-1] = 360.0-dln
    
    lt1*=np.pi/180.0
    ln1*=np.pi/180.0
    
    darea = np.zeros((len(lt),len(ln)))
    for jlat in range(0,len(lt)):
        for jlon in range(0,len(ln)):
            dln = ln1[jlon+1]-ln1[jlon]
            darea[jlat,jlon] = (np.sin(lt1[jlat])-np.sin(lt1[jlat+1]))*dln
    
    svar = variable*darea
    if mean:
        outvar = sumop(svar)/sumop(darea)
    else:
        outvar = sumop(svar) * radius**2
    
    return outvar


def cspatialmath(variable,lat=None,lon=None,file=None,mean=True,time=None,
               ignoreNaNs=True,lev=None,radius=6.371e6,poles=False):
    """Compute spatial means or sums of data, but optionally don't go all the way to the poles.

    Sometimes, saying that the latitudes covered go all the way to :math:`\pm90^\circ` results in
    errors, and accurate accounting requires excluding the poles themselves. This function
    is identical to spatialmath, except that it provides that option.

    Parameters
    ----------
    variable : str, numpy.ndarray
        The variable to operate on. Can either be a data array, or the name of a variable. If the latter, file must be specified.
    lat,lon : numpy.ndarray, optional
        Latitude and longitude arrays. If file is provided and lat and lon are not, they will be
        extracted from the file.
    file : str, optional
        Path to a NetCDF output file to open and extract data from.
    mean : bool, optional
        If True, compute a global mean. If False, compute a global sum.
    time : int, optional
        The time index on which to slice. If unspecified, a time average will be returned.
    ignoreNaNs : bool, optional
        If True, use NaN-safe numpy operators.
    lev : int, optional
        If set, slice a 3D spatial array at the specified level.
    radius : float, optional
        Radius of the planet in meters. Only used if mean=False.
    poles : bool, optional
        If False (default), exclude the poles.
        
    Returns
    -------
    float

    """
    
    if ignoreNaNs:
        sumop = np.nansum
        meanop = np.nanmean
    else:
        sumop = np.sum
        meanop = np.mean
        
    if file:
        ln,lt,variable = parse(file,variable,lat=lat,lon=lon)
        
    else:
        if type(lat)==type(None) or type(lon)==type(None):
            raise DimensionError("Need to provide latitude and longitude data")
        ln=lon
        lt=lat
    variable = make2d(variable,time=time,lev=lev,ignoreNaNs=ignoreNaNs)
    
    lt1 = np.zeros(len(lt)+1)
    dlt1 = abs(np.diff(lt)[0])
    dlt2 = abs(np.diff(lt)[-1])
    lt1[0] = lt[0]+0.5*dlt1
    for n in range(0,len(lt)-1):
        lt1[n+1] = 0.5*(lt[n]+lt[n+1])
    lt1[-1] = lt[-1]-0.5*dlt2
    if poles:
        lt1[0] = 90.0
        lt1[-1] = -90.0
    
    dln = np.diff(ln)[0]
    ln1 = np.zeros(len(ln)+1)
    ln1[0] = ln[0]-dln*0.5
    for n in range(0,len(ln)-1):
        ln1[n+1] = 0.5*(ln[n]+ln[n+1])
    ln1[-1] = ln[-1]+dln*0.5
    
    lt1*=np.pi/180.0
    ln1*=np.pi/180.0
    
    darea = np.zeros((len(lt),len(ln)))
    for jlat in range(0,len(lt)):
        for jlon in range(0,len(ln)):
            dln = ln1[jlon+1]-ln1[jlon]
            darea[jlat,jlon] = (np.sin(lt1[jlat])-np.sin(lt1[jlat+1]))*dln
    
    svar = variable*darea
    if mean:
        outvar = sumop(svar)/sumop(darea)
    else:
        outvar = sumop(svar) * radius**2
    
    return outvar

def wrap2d(var):
    '''Add one element to the longitude axis to allow for wrapping'''
    newvar = np.zeros(np.array(var.shape)+np.array((0,1)))
    newvar[:,:-1] = var[:,:]
    newvar[:,-1] = var[:,0]
    return newvar

    
def streamfxn(file,time=None):
    '''Return the streamfunction

    Parameters
    ----------
    file : str
        Path to an ExoPlaSim NetCDF output file.
        
    Returns
    -------
    numpy.ndarray
        The streamfunction for the given file.
    '''
    ln,lt,levs=parse(file,"lev")
    plevs = levs*spatialmath("ps",file=file,time=time)
    ln,lt,va = parse(file,"va")
    va = make2d(va,lon="mean",time=time)
    strf = np.zeros(va.shape)
    pref = 2*np.pi*6.371e6*np.cos(lt*np.pi/180.0)/9.81
    ps = np.array([0.0,]+list(plevs))
    vas = np.zeros(np.array(va.shape)+np.array((1,0)))
    vas[1:,:] = va[:,:]
    for k in range(0,len(plevs)):
        strf[k,:] = pref[:]*np.trapz(vas[0:k+1,:],x=ps[0:k+1],axis=0)
    return strf

    