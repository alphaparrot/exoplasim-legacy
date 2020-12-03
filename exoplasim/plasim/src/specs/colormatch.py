import numpy as np
from scipy import interpolate

#From http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
widegamut = np.array([[ 1.4628067, -0.1840623, -0.2743606],
                      [-0.5217933,  1.4472381,  0.0677227],
                      [ 0.0349342, -0.0968930,  1.2884099]])




def _loadcie():
    wvl,xx,yy,zz = np.loadtxt("cmf.csv",unpack=True,delimiter=',')
    return wvl,xx,yy,zz
    
def interpolant(wvl):
    w0,xx,yy,zz = _loadcie()
    fx = interpolate.interp1d(w0,xx)
    fy = interpolate.interp1d(w0,yy)
    fz = interpolate.interp1d(w0,zz)
    imin = np.where(wvl>np.amin(w0))[0][0]
    imax = np.where(wvl<np.amax(w0))[0][-1]
    wn = wvl[imin:imax+1]
    xn = fx(wn)
    yn = fy(wn)
    zn = fz(wn)
    return (xn,yn,zn)
    
def makexyz(wvl,spec,interpolant=None):
    if np.amin(wvl)<1.0e-3: #probably meters not nanometers
        wvl*=1.0e9
    w0,xx,yy,zz = _loadcie()
    imin = np.where(wvl>np.amin(w0))[0][0]
    imax = np.where(wvl<np.amax(w0))[0][-1]
    wn = wvl[imin:imax+1]
    specn = spec[imin:imax+1]
    if not interpolant:
        fx = interpolate.interp1d(w0,xx)
        fy = interpolate.interp1d(w0,yy)
        fz = interpolate.interp1d(w0,zz)
        
        xn = fx(wn)
        yn = fy(wn)
        zn = fz(wn)
    else:
        xn = interpolant[0]
        yn = interpolant[1]
        zn = interpolant[2]
    
    XI = np.trapz(xn*specn,x=wn)
    YI = np.trapz(yn*specn,x=wn)
    ZI = np.trapz(zn*specn,x=wn)
    xyzmin = np.amin((XI,YI,ZI))
    if xyzmin<0:
        XI -=xyzmin
        YI -=xyzmin
        ZI -=xyzmin
    if (XI+YI+ZI)>0:
        xnu = XI/(XI+YI+ZI)
        ynu = YI/(XI+YI+ZI)
        znu = 1.0-(xnu+ynu)
    else:
        xnu=0
        ynu=0
        znu=0
    
    return xnu,ynu,YI

def xyz2rgb(x,y,normalization):
    
    z = 1-(x+y)
    
    r = np.sum(widegamut[0,:]*x)*normalization
    g = np.sum(widegamut[1,:]*y)*normalization
    b = np.sum(widegamut[2,:]*z)*normalization
    cmax = np.amax((r,g,b))
    
    #return r/cmax,g/cmax,b/cmax
    return r,g,b


def spec2rgb(wvl,spec,normalization=None):
    x,y,I = makexyz(wvl,spec)
    if normalization:
        norm = normalization
    else:
        norm = 1.0
    
    r,g,b = xyz2rgb(x,y,norm)
    return r,g,b


def specs2rgb(wvl,specs):
    interpol = interpolant(wvl)
    intensities = np.zeros((len(specs),3))
    for n in range(0,len(specs)):
        intensities[n,:] = makexyz(wvl,specs[n,:],interpolant=interpol)
    norms = intensities[:,2]/np.amax(intensities[:,2])
    colors = np.zeros((len(specs),3))
    for n in range(0,len(specs)):
        colors[n,:] = xyz2rgb(intensities[n,0],intensities[n,1],norms[n])
    
    return colors
    
    