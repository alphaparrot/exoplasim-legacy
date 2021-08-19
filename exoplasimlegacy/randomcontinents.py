import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import os, sys
import argparse as ag
import matplotlib.colors as colors
from pathlib import Path
    


def _wrap2d(datd,vals):
    modf=np.zeros(datd.ndim,dtype=int)
    modf[-1]=1
    dd=np.zeros(datd.shape+modf)
    dd[:,0:datd.shape[-1]]=datd
    dd[:,datd.shape[-1]]=vals
    return dd

def writeSRA(name,kcode,field,NLAT,NLON):
    """Write a lat-lon field to a formatted .sra file
    
    Parameters
    ----------
    name : str
        The name with which to label this map
    kcode : int
        The integer map code for specifying what kind of boundary file this is (see the PlaSim documentation for more details)
    field : numpy.ndarray
        The map to write to file. Should have the dimensions (NLAT,NLON).
    NLAT : int
        The number of latitudes
    NLON : int
        The number of longitudes
        
    """
    label=name+'_surf_%04d.sra'%kcode
    header=[kcode,0,20170927,0,NLON,NLAT,0,0]
    fmap = field.reshape((int(NLAT*NLON//8),8))
    sheader = ''
    for h in header:
        sheader+=" %11d"%h
    
    lines=[]
    i=0
    while i<NLAT*NLON/8:
        l=''
        for n in fmap[i,:]:
            l+=' %9.3f'%n
        lines.append(l)
        i+=1

    text=sheader+'\n'+'\n'.join(lines)+'\n' 

    f=open(label,'w')
    f.write(text)
    f.close()
    
    
def writePGM(name,heightfield):    
    """Write a lat-lon field to a .pgm image file (usually topo field)
    
    Parameters
    ----------
    name : str
        The name with which to label this map
    heightfield : numpy.ndarray
        The 2-D map to write to file.
        
    """
    shape=heightfield.shape
    filetext = ("P2\n"+
                "# Heightfield map for %s planet\n"%name+
                "%d %d\n"%(shape[1],shape[0])+
                "65535\n")
    img = ((np.round((heightfield/heightfield.max())*65535)).astype(int)).astype(str)
    for k in range(shape[0]):
        filetext+=' '.join(img[k,:])+'\n'
    with open(name+".pgm","w") as fw:
        fw.write(filetext)

def main():
    """Command-line tool to randomly generate continents up to specified land fraction. Topography optional.
    
    Do not invoke as an imported function; must run directly.

**Options**
        -z,--topo   
            Generate topographical geopotential map
        -c,--continents   
            Number of continental cratons
        -f,--landfraction   
            Land fraction
        -n,--name   
            Assign a name for the planet
        -m,--maxz   
            Maximum elevation in km assuming Earth gravity
        --nlats   
            Number of latitudes (evenly-spaced)--will also set longitudes (twice as many). If unset, PlaSim latitudes and longitudes will be used (T21 resolution)"
        -l,--hemispherelongitude   
            Confine land to a hemisphere centered on a given longitude
        -o,--orthographic   
            Plot orthographic projections centered on hemispherelongitude 

    Yields
    ------
    name_surf_0172.sra
        Land mask SRA file
    name_surf_0129.sra (optional)
        Topography geopotential SRA file (if requested)

    """
    parser = ag.ArgumentParser(description="Randomly generate continents up to a specified land-fraction. Topography optional.")
    parser.add_argument("-z","--topo",action="store_true",help="Generate topographical geopotential map",)
    parser.add_argument("-c","--continents",type=int,default=7,help="Number of continental cratons")
    parser.add_argument("-f","--landfraction",type=float,default=0.29,help="Land fraction")
    parser.add_argument("-n","--name",default="Alderaan",help="Assign a name for the planet")
    parser.add_argument("-m","--maxz",default=10.0,type=float,help="Maximum elevation in km assuming Earth gravity")
    if os.path.exists("T21.nc"):
        parser.add_argument("--nlats",type=int,help="Number of latitudes (evenly-spaced)--will also set longitudes (twice as many). If unset, PlaSim latitudes and longitudes will be used (T21 resolution)")
    else:
        parser.add_argument("--nlats",default=32,type=int,help="Number of latitudes (evenly-spaced)--will also set longitudes (twice as many).")
    parser.add_argument("-l","--hemispherelongitude",type=float,default=np.nan,help="Confine land to a hemisphere centered on a given longitude")
    parser.add_argument("-o","--orthographic",action="store_true",help="Plot orthographic projections centered on hemispherelongitude")
    args = parser.parse_args()
    
    if not args.nlats:
        import netCDF4 as nc
        dims = nc.Dataset("T21.nc","r")
    
        lts = dims.variables['lat'][:]
        lns = dims.variables['lon'][:]
    else:
        lts = np.linspace(90,-90,num=args.nlats+2)[1:-1]
        lns = np.linspace(0,360,num=args.nlats*2+1)[:-1]
    lons, lats = np.meshgrid(lns,lts)
    
    minlon=0.0
    maxlon=360.0
    lonrange=360.0
    l0 = 0.5*(minlon+maxlon)
    wraplon=False
    if np.isfinite(args.hemispherelongitude):
        l0 = args.hemispherelongitude
        minlon=args.hemispherelongitude-90.0
        maxlon=args.hemispherelongitude+90.0
        lonrange=180.0
        if minlon<0:
            minlon+=360.0
            wraplon=True
        if maxlon>360.0:
            maxlon-=360.0
            wraplon=True
        if wraplon:
            z1 = minlon
            z2 = maxlon
            minlon = min(z1,z2)
            maxlon = max(z1,z2)
    
    
    darea = np.zeros((len(lts),len(lns)))
    NLAT = len(lts)
    lts1 = np.zeros(NLAT)
    lts2 = np.zeros(NLAT)
    lts1[0] = 0.5*np.pi
    lts1[NLAT-1] = 0.5*(lts[NLAT-2]+lts[NLAT-1])*np.pi/180.0
    lts2[0] = 0.5*(lts[0]+lts[1])*np.pi/180.0
    lts2[NLAT-1] = -0.5*np.pi
    for jlat in range(1,NLAT-1):
        lts1[jlat] = 0.5*(lts[jlat-1]+lts[jlat])*np.pi/180.0
        lts2[jlat] = 0.5*(lts[jlat+1]+lts[jlat])*np.pi/180.0
    NLON = len(lns)
    for jlat in range(0,NLAT):
        darea[jlat,:] = 0.5/NLON*(np.sin(lts1[jlat])-np.sin(lts2[jlat]))
    
    globalarea = np.sum(darea)
    
    hlns = np.zeros(NLON+1)
    hlts = np.zeros(NLAT+1)
    hlts[0] = 90.0
    hlts[-1] = -90.0
    hlts[1:-1] = 0.5*(lts[:-1]+lts[1:])
    hlns[0] = lns[0]-0.5*(lns[1]-lns[0])
    hlns[-1] = lns[-1]+0.5*(lns[-1]-lns[-2])
    hlns[1:-1] = 0.5*(lns[:-1]+lns[1:])
    hlons, hlats = np.meshgrid(hlns,hlts)
    
    rlats = hlats*np.pi/180.0
    rlons = hlons*np.pi/180.0
    rl0 = l0*np.pi/180.0
    rcoord = 2*np.arcsin(np.sqrt(np.sin(0.5*rlats)**2+np.cos(rlats)*np.sin(0.5*(rlons-rl0))**2))
    thetacoord = np.arctan2(np.cos(rlats)*np.sin(rlons-rl0),np.sin(rlats))
    thetacoord[thetacoord<0] += 2*np.pi
    
    latsw=_wrap2d(lats,lats[:,0])
    lonsw=_wrap2d(lons,360.0)
    hlatsw=_wrap2d(hlats,hlats[:,0])
    hlonsw=_wrap2d(hlons,360.0)
    
    ncontinents = args.continents #There's no guarantee you actually get this many if it's >1, since continents can merge
    
    landfraction = args.landfraction
    
    name=args.name
    
    grid = np.zeros((NLAT,NLON))
    wgrid = np.zeros((NLAT+2,NLON+2))
    cratons = np.zeros((NLAT+2,NLON+2))+np.nan
    seams = np.zeros((NLAT,NLON))
    landarea = 0.0
    
    ilons,ilats = np.meshgrid(range(NLON),range(NLAT))
    ilns1 = ilons.flatten()
    ilts1 = ilats.flatten()
    
    for c in range(ncontinents):
        while True:
            #theta = np.random.uniform()*360.0
            #phi = np.arcsin(1-2*np.random.uniform())*180.0/np.pi
            idx = np.random.choice(range(NLAT*NLON),p=darea.flatten())
            #ilt = np.argmin(abs(phi-lts))
            #iln = np.argmin(abs(theta-lns))
            ilt = ilts1[idx]
            iln = ilns1[idx]
            keeppoint = True
            if wraplon:
                if lns[iln]>minlon and lns[iln]<maxlon:
                    keeppoint=False
            else:
                if lns[iln]<minlon or lns[iln]>maxlon:
                    keeppoint=False
            if grid[ilt,iln]<0.5 and keeppoint:
                grid[ilt,iln] = 1.0
                cratons[ilt+1,iln+1] = c
                seams[ilt,iln] = 1.0
                landarea += darea[ilt,iln]
                break
        wgrid[1:-1,1:-1] = grid[:,:]
        wgrid[0,1:-1] = grid[0,:][::-1]
        wgrid[-1,1:-1] = grid[-1,:][::-1]
        wgrid[:,0] = wgrid[:,-2]
        wgrid[:,-1] = wgrid[:,1]
        cratons[0,1:-1] = cratons[1,1:-1][::-1]
        cratons[-1,1:-1] = cratons[-2,1:-1][::-1]
        cratons[:,0] = cratons[:,-2]
        cratons[:,-1] = cratons[:,1]
        
    history=[]
    while landarea<=landfraction-np.nanmin(darea):
        while True:
            theta = np.random.uniform()*lonrange
            if wraplon:
                theta += maxlon
                if theta>360.0:
                    theta-=360.0
            else:
                theta += minlon
                if theta>360.0:
                    theta-=360.0
            phi = np.arcsin(1-2*np.random.uniform())*180.0/np.pi
            ilt = np.argmin(abs(phi-lts))
            iln = np.argmin(abs(theta-lns))
            if grid[ilt,iln]<0.5 and np.sum(wgrid[ilt:ilt+3,iln:iln+3])>0.5:
                goodtogo=True
                if np.sum(wgrid[ilt:ilt+3,iln:iln+3])<5.0: #Try to bias towards land-locked ocean cells
                    if np.random.uniform()<0.9:
                        goodtogo=False
                if goodtogo:
                    grid[ilt,iln] = 1.0
                    crsq = cratons[ilt:ilt+3,iln:iln+3]
                    cratons[ilt+1,iln+1] = np.nanmean(crsq[crsq>0.0])
                    #crsqm = crsq.min()
                    #if crsqm==0.0:
                    #    crsqm = crsq[crsq>0.0].min()
                    #if crsq.max()-crsqm > 0.0:
                        #seams[ilt,iln]=1.0
                    landarea += darea[ilt,iln]
                    history.append(landarea)
                    break
        wgrid[1:-1,1:-1] = grid[:,:]
        wgrid[0,1:-1] = grid[0,:][::-1]
        wgrid[-1,1:-1] = grid[-1,:][::-1]
        wgrid[:,0] = wgrid[:,-2]
        wgrid[:,-1] = wgrid[:,1]
        cratons[0,1:-1] = cratons[1,1:-1][::-1]
        cratons[-1,1:-1] = cratons[-2,1:-1][::-1]
        cratons[:,0] = cratons[:,-2]
        cratons[:,-1] = cratons[:,1]
    
    
    print(args.orthographic)
    
    if not args.orthographic:
        tm = plt.pcolormesh(hlons,hlats,grid,cmap='gist_earth',vmin=-0.3,vmax=2.5)
        plt.xlabel('Degrees Longitude')
        plt.ylabel('Degrees Latitude')
        plt.ylim(np.amin(lts),np.amax(lts))
        plt.xlim(np.amin(lns),np.amax(lns))
        plt.title("Continents")
        plt.savefig(name+"_lsm.png",bbox_inches='tight')
        plt.savefig(name+"_lsm.pdf",bbox_inches='tight')
        plt.close('all')  
    else:
        iln0 = np.argmin(abs(lns-minlon))
        iln1 = np.argmin(abs(lns-maxlon))
        shift = iln0
        if wraplon:
            shift = iln1
        x = np.roll(thetacoord,-shift,axis=1)[:,:NLON/2+1]
        y = np.roll(rcoord,-shift,axis=1)[:,:NLON/2+1]
        z = np.roll(grid,-shift,axis=1)[:,:NLON/2]
        fig,ax=plt.subplots(subplot_kw={"projection":"polar"},figsize=(9,9))
        ax.set_theta_zero_location('N')
        ax.pcolormesh(x,np.sin(y),z,cmap='gist_earth',vmin=-0.3,vmax=2.5)
        ax.set_rticks([])
        ax.set_thetagrids([])
        plt.title("Continents")
        plt.savefig(name+"_lsm.png",bbox_inches='tight')
        plt.savefig(name+"_lsm.pdf",bbox_inches='tight')
        plt.close('all')  
        
    
    writeSRA(name,172,grid,NLAT,NLON)
    
    
    plt.close('all')
    if not args.orthographic:
        t = plt.pcolormesh(hlons,hlats,cratons[1:-1,1:-1],cmap='gist_ncar',vmin=-0.3,vmax=ncontinents+2)
        plt.xlabel('Degrees Longitude')
        plt.ylabel('Degrees Latitude')
        plt.ylim(np.amin(lts),np.amax(lts))
        plt.xlim(np.amin(lns),np.amax(lns))
        plt.title("Continental Cratons")
        plt.savefig(name+"_cratons.png",bbox_inches='tight')
        plt.savefig(name+"_cratons.pdf",bbox_inches='tight')
    
        plt.close('all')
    else:
        iln0 = np.argmin(abs(lns-minlon))
        iln1 = np.argmin(abs(lns-maxlon))
        shift = iln0
        if wraplon:
            shift = iln1
        x = np.roll(thetacoord,-shift,axis=1)[:,:NLON/2+1]
        y = np.roll(rcoord,-shift,axis=1)[:,:NLON/2+1]
        z = np.roll(cratons[1:-1,1:-1],-shift,axis=1)[:,:NLON/2]
        fig,ax=plt.subplots(subplot_kw={"projection":"polar"},figsize=(9,9))
        ax.set_theta_zero_location('N')
        ax.pcolormesh(x,np.sin(y),z,cmap='gist_ncar',vmin=-0.3,vmax=ncontinents+2)
        ax.set_rticks([])
        ax.set_thetagrids([])
        plt.title("Continental Cratons")
        plt.savefig(name+"_cratons.png",bbox_inches='tight')
        plt.savefig(name+"_cratons.pdf",bbox_inches='tight')
        plt.close('all')  
    
    if args.topo:
        seeds = np.copy(seams)
        gcratonsx,gcratonsy = np.gradient(cratons[1:-1,1:-1],lts,lns)
        seams[:] = np.sqrt(gcratonsx**2+gcratonsy**2)
        seams[np.isnan(seams)]=0.0
        
        g0 = 9.80665
        if np.nanmax(seams)>0.0:
            geopotential = g0*(seams/np.nanmax(seams))*args.maxz*1000.0
        else:
            geopotential = g0*(grid*0.1+seeds*3.0)*args.maxz*1000.0
            seams[:] = seeds[:]
        geopotential[grid==0.0] = np.nan
        
        if not args.orthographic:
            tm = plt.pcolormesh(hlons,hlats,grid+seams*3.0,cmap='plasma')
            plt.title("Craton Seams")
            plt.xlabel("Degrees Longitude")
            plt.ylabel("Degrees Latitude")
            plt.ylim(np.amin(lts),np.amax(lts))
            plt.xlim(np.amin(lns),np.amax(lns))
            plt.savefig(name+"_seams.png",bbox_inches='tight')
            plt.savefig(name+"_seams.pdf",bbox_inches='tight')
            plt.close('all')
        else:
            iln0 = np.argmin(abs(lns-minlon))
            iln1 = np.argmin(abs(lns-maxlon))
            shift = iln0
            if wraplon:
                shift = iln1
            x = np.roll(thetacoord,-shift,axis=1)[:,:NLON/2+1]
            y = np.roll(rcoord,-shift,axis=1)[:,:NLON/2+1]
            z = np.roll(grid+seams*3.0,-shift,axis=1)[:,:NLON/2]
            fig,ax=plt.subplots(subplot_kw={"projection":"polar"},figsize=(9,9))
            ax.set_theta_zero_location('N')
            t=ax.pcolormesh(x,np.sin(y),z,cmap='plasma',vmin=-0.3,vmax=2.5)
            ax.set_rticks([])
            ax.set_thetagrids([])
            plt.colorbar(t)
            plt.title("Craton Seams")
            plt.savefig(name+"_seams.png",bbox_inches='tight')
            plt.savefig(name+"_seams.pdf",bbox_inches='tight')
            plt.close('all')  
    
        hlnsz = np.linspace(hlns[0],hlns[-1],num=max(400,NLON*2))
        hltsz = np.linspace(hlts[0],hlts[-1],num=max(200,NLAT*2))
        lnsz = 0.5*(hlnsz[:-1]+hlnsz[1:])
        ltsz = 0.5*(hltsz[:-1]+hltsz[1:])
    #lnsz = np.linspace(lns[0],lns[-1],num=128)
    #ltsz = np.linspace(lts[0],lts[-1],num=64)
        geo = np.copy(geopotential)
        geo[np.isnan(geopotential)] = 0.0
        
        kx=3
        ky=3
        if NLAT>128:
            kx=1
            ky=1
        
        print(lns.min(),lns.max(),lts.min(),lts.max())
        print(lnsz[0],lnsz[-1],ltsz[-1],ltsz[0])
        geozspline = interpolate.RectBivariateSpline(lns,lts[::-1],np.transpose(geo[::-1,:]),bbox=[lnsz[0],lnsz[-1],ltsz[-1],ltsz[0]],kx=kx,ky=ky)
        geoz = np.transpose(geozspline(lnsz,ltsz[::-1]))
        geoz = geoz[::-1,:]
        contzspline = interpolate.RectBivariateSpline(lns,lts[::-1],np.transpose(grid[::-1,:]),bbox=[lnsz[0],lnsz[-1],ltsz[-1],ltsz[0]],kx=kx,ky=ky)
        contz = np.transpose(contzspline(lnsz,ltsz[::-1]))
        contz = contz[::-1,:]
        
        topo = np.copy(geoz)+contz*g0*10.0 #Default lowlands of 10 meters above sea level
    
        maxiters = 10*int(np.sqrt(NLAT/32))
    
        NLATZ = len(ltsz)
        NLONZ = len(lnsz)
    
        dl1 = NLONZ/2+2
        dl2 = NLONZ/2
        dl3 = NLONZ/2-2
    
        for i in range(0,maxiters):
            topo = np.maximum(topo,geoz)
            topo[np.isnan(geoz)]=0.0
            nutop = np.zeros(contz.shape)
            for jlat in range(0,NLATZ):
                j1 = jlat-1
                j2 = jlat+1
                flip1=False
                flip2=False
                if j1<0:
                    j1=0
                    flip1=True
                if j2==NLATZ:
                    j2=NLATZ-1
                    flip2=True
                for jlon in range(0,NLONZ):
                    l1 = jlon-1
                    l2 = jlon+1
                    if l1<0:
                        l1=NLONZ-1
                    if l2==NLONZ:
                        l2=0
                    c1 = topo[j1,(l1  +flip1*dl1)%NLONZ]
                    c2 = topo[j1,(jlon+flip1*dl2)%NLONZ]
                    c3 = topo[j1,(l2  +flip1*dl3)%NLONZ]
                    c4 = topo[jlat,l1]
                    c5 = topo[jlat,l2]
                    c6 = topo[j2,(l1  +flip2*dl1)%NLONZ]
                    c7 = topo[j2,(jlon+flip2*dl2)%NLONZ]
                    c8 = topo[j2,(l2  +flip2*dl3)%NLONZ]
                    nutop[jlat,jlon] = 0.2*topo[jlat,jlon]+np.nanmean([c1,c2,c3,c4,c5,c6,c7,c8])*0.8
            topo[:] = nutop[:]
        topospline = interpolate.interp2d(lnsz,ltsz[::-1],(topo[::-1,:]),kind='linear')
        dtopo = (topospline(lns,lts[::-1]))
        dtopo = dtopo[::-1,:]
        dtopo[np.isnan(geopotential)] = 0.0
        
        writeSRA(name,129,dtopo,NLAT,NLON)
        
        if not args.orthographic:
            plt.close('all')
            t=plt.pcolormesh(hlons,hlats,dtopo,cmap='gist_earth',norm=colors.LogNorm(vmin=10.0))
            plt.xlabel('Degrees Longitude')
            plt.ylabel('Degrees Latitude')
            plt.ylim(np.amin(lts),np.amax(lts))
            plt.xlim(np.amin(lns),np.amax(lns))
            plt.title("Topography")
            plt.colorbar(t,label="Geopotential [m$^2$/s$^2$]")
            plt.savefig(name+"_geoz.png",bbox_inches='tight')
            plt.savefig(name+"_geoz.pdf",bbox_inches='tight')
        else:
            iln0 = np.argmin(abs(lns-minlon))
            iln1 = np.argmin(abs(lns-maxlon))
            shift = iln0
            if wraplon:
                shift = iln1
            x = np.roll(thetacoord,-shift,axis=1)[:,:NLON/2+1]
            y = np.roll(rcoord,-shift,axis=1)[:,:NLON/2+1]
            z = np.roll(dtopo,-shift,axis=1)[:,:NLON/2]
            fig,ax=plt.subplots(subplot_kw={"projection":"polar"},figsize=(9,9))
            ax.set_theta_zero_location('N')
            t=ax.pcolormesh(x,np.sin(y),z,cmap='gist_earth',norm=colors.LogNorm(vmin=10.0))
            ax.set_rticks([])
            ax.set_thetagrids([])
            plt.title("Topography")
            plt.colorbar(t,label="Geopotential [m$^2$/s$^2$]")
            plt.savefig(name+"_geoz.png",bbox_inches='tight')
            plt.savefig(name+"_geoz.pdf",bbox_inches='tight')
            plt.close('all')  
        
        hf = np.copy(dtopo)
        hf[grid>0.5] += 1000.0
        
        writePGM(name,hf)
        
if __name__=="__main__" and (Path(sys.argv[0]).name!="sphinx-build" and 
                             Path(sys.argv[0]).name!="build.py"):
    main()