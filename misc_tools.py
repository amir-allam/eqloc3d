import subprocess
import time
import struct
import os,sys
#import utm
sys.path.append('%s/data/python' % os.environ['ANTELOPE'])
from antelope.datascope import closing, dbopen
from antelope.stock import str2epoch, epoch2str
from numpy import * #Should probably be more specific
from matplotlib import pyplot as plt
from array import array #This is for fast io when writing binary 
from collections import defaultdict

#Read ASCII arrival time file output by FMM code

def num(s):
#Convert a string to a number, choosing int or float automatically
# SLOW, don't use for large lists
    try:
        return int(s)
    except ValueError:
        return float(s)

def load_faults():
#Load the california fault map
    fnam='cal_faults.dat'
    fid=open(fnam,'r')
    a=fid.readlines()
    faultlon=[];faultlat=[]
    for line in a:
        tmp=line.strip().split()
        if isnan(num(tmp[0])):
            faultlon.append(nan)
            faultlat.append(nan)
        else:
            faultlon.append(num(tmp[0]))
            faultlat.append(num(tmp[1]))
    print 'Faults loaded...'
    return faultlon,faultlat

class TomoDD_MOD():
    #Contains the velocity and metadata from a TomoDD MOD file which looks like
    #bld NX NY NZ
    # X1 X2... XN
    # Y1 Y2... YN
    # Z1 Z2... ZN
    #V111 V211 V311 ... VN11
    #V121 V221 V321 ... VN21
    #...
    #V1N1 V2N1 V3N1 ... VNN1
    #V112 V212 V312 ... VN12
    #...
    #V11N V21N V31N ... VN1N
    #...
    # THEN, the entire V matrix is repeated for Vp/Vs
    # so, total file length is NY*NZ*2 + 4
    #
    # NOTE: Left-handed coordinates; depth is positive
    def __init__(self):
        pass
    def read(self,fnam='MOD'):
        #Read a TomoDD MOD velocity file. It is an ascii file that looks like:
        #Open and read header
        self.fnam=fnam
        fid = open(fnam,'r')
        tmp=fid.readline().strip().split()
        self.bld=num(tmp[0])
        self.nx=num(tmp[1])
        self.ny=num(tmp[2])
        self.nz=num(tmp[3]) #Number of grid points
        tmp=fid.readline().strip().split()
        self.qx=asarray([float(ii) for ii in tmp]) 
        tmp=fid.readline().strip().split()
        self.qy=asarray([float(ii) for ii in tmp])
        tmp=fid.readline().strip().split()
        self.qz=asarray([float(ii) for ii in tmp])
        #Loop to create a numpy matrix
        self.vel=empty((self.nx,self.ny,self.nz)) #Check that this imports properly from numpy
        for iz in range(self.nz):
            for iy in range(self.ny):
                tmp=fid.readline().strip().split()
                for ix in range(self.nx):
                    self.vel[ix,iy,iz]=float(tmp[ix])
        self.velrat=empty((self.nx,self.ny,self.nz))
        for iz in range(self.nz):
            for iy in range(self.ny):
                tmp=fid.readline().strip().split()
                for ix in range(self.nx):
                    self.velrat[ix,iy,iz]=float(tmp[ix])
        fid.close()
    def write(self,outfnam='out.MOD'): #Write to MOD format
        self.outfnam=outfnam
        fid=open(outfnam,'w')
        outs=str(self.bld)+' '+str(self.nx)+' '+str(self.ny)+' '+str(self.nz)+'\n'
        fid.write(outs)
        outs=''.join(str(ii)+' ' for ii in self.qx)+'\n'
        fid.write(outs)
        outs=''.join(str(ii)+' ' for ii in self.qy)+'\n'
        fid.write(outs)
        outs=''.join(str(ii)+' ' for ii in self.qz)+'\n'
        fid.write(outs)
        for iz in range(self.nz):
            for iy in range(self.ny):
                outs=''.join('{1:.{0}f} '.format(3,ii)  for ii in self.vel[:,iy,iz])+'\n'
                fid.write(outs)
        for iz in range(self.nz):
            for iy in range(self.ny):
                outs=''.join('{1:.{0}f} '.format(3,ii)  for ii in self.velrat[:,iy,iz])+'\n'
                fid.write(outs)
        fid.close()

class Fmm_vgrids():
# vgrids.in is the FMM velocity file format. Looks like:
#    ngrids ntypes
#    nradius nlat nlon     !number of points
#    dradius dlat dlon     !grid spacings; uniform in each direction
#    oradius olat olon     !origin of the grid
#    V(1,1,1)
#    V(2,1,1)
#    V(1,2,1)
#    V(1,1,2)
#    etc...
#
#       NOTE: Right-handed; oradius is somewhere deep in the earth
#
# ngrids is the number of volumes. We generally use 1
# ntypes will be 1 for just P, 2 for P and S
    def __init__(self):
        pass
    def read(self,fnam='vgrids.in'): #Create a matrix of velocities
        self.fnam=fnam
        fid = open(fnam,'r')
        tmp=fid.readline().strip().split()
        self.ngrids,self.ntypes=num(tmp[0]),num(tmp[1])
        tmp=fid.readline().strip().split()
        self.nrad,self.nlat,self.nlon=num(tmp[0]),num(tmp[1]),num(tmp[2])
        tmp=fid.readline().strip().split()
        self.drad,self.dlat,self.dlon=num(tmp[0]),num(tmp[1]),num(tmp[2])
        tmp=fid.readline().strip().split()
        self.orad,self.olat,self.olon=num(tmp[0]),num(tmp[1]),num(tmp[2])
        #Loop to create a numpy matrix
        self.vel=empty((self.nlon,self.nlat,self.nrad))
        for irad in range(self.nrad):
            for ilat in range(self.nlat):
                for ilon in range(self.nlon):
                    tmp=fid.readline().strip().split()
                    self.vel[ilon,ilat,irad]=float(tmp[0])
        #THERE SHOULD BE MORE STATEMENTS HERE IN CASE NTYPES,NGRIDS!= 1
        fid.close()
    def write(self,outfnam='out.vgrids'): #Write to vgrids format
        self.outfnam=outfnam
        fid=open(outfnam,'w')
        outs=str(self.ngrids)+' '+str(self.ntypes)+'\n'
        fid.write(outs)
        outs=str(self.nrad)+' '+str(self.nlat)+' '+str(self.nlon)+'\n'
        fid.write(outs)
        outs=str(self.drad)+' '+str(self.dlat)+' '+str(self.dlon)+'\n'
        fid.write(outs)
        outs=str(self.orad)+' '+str(self.olat)+' '+str(self.olon)+'\n'
        fid.write(outs)
        for ilon in range(self.nlon):
            for ilat in range(self.nlat):
                for irad in range(self.nrad):
                    fid.write('{1:.{0}f}'.format(3,self.vel[ilon,ilat,irad])+'\n')
        fid.close()

class Interface():
#A 2D matrix containing interface information, e.g., topography or moho depth
    def __init__(self):
        pass
    def set_data(self,x,y,z):
    #This builds the class ad-hoc. x,y,z are all numpy arrays. x,y are vectors, z is 2D
        self.x,self.y,self.z=x,y,z
        deg_to_rad=math.acos(-1)/180.0  #math.acos(-1) is pi
        self.nx,self.ny=len(x),len(y) #total number of x and y coordinates
        self.ox,self.oy=min(x),min(y)   #origins in lat,lon
        self.dx=self.x[1]-self.x[0]     #Spacings in degrees
        self.dy=self.y[1]-self.y[0]
        self.ox_rad=self.ox*deg_to_rad  #Origins in radians
        self.oy_rad=self.oy*deg_to_rad
        self.dx_rad=self.dx*deg_to_rad  #Spacings in radians
        self.dy_rad=self.dy*deg_to_rad
    def read_topo_netcdf(self,fnam,subsamp=1):
    #Read a netcdf formatted topography file
        from scipy.io import netcdf_file as netcdf
        topo=netcdf(fnam,'r')
        x=topo.variables['x'][::subsamp] #vector
        y=topo.variables['y'][::subsamp] #vector
        z=topo.variables['z'][::subsamp,::subsamp] #matrix
        self.subsamp=subsamp
        self.set_data(x,y,z)
    def interp(self,xi,yi):
    #Find the values at points xi,yi by cubic spline
        import scipy.interpolate
        gnn=scipy.interpolate.RectBivariateSpline(self.x,self.y,self.z.transpose())
        self.x_unint=self.x #Remember non-interpolated values
        self.y_unint=self.y
        self.z_unint=self.z
        zi=gnn.__call__(xi,yi) #the actual interpolation
        self.set_data(xi,yi,zi)
    def write_fmm(self,fnam='out_interfaces',ifappend=False):
    #Write to the fmm ascii format. ifappend should be set to True for any interface after the first; this will only write the z values and not the header
        if not ifappend:
            fid=open(fnam,'w')
            fid.write('2\n') # THIS ASSUMES THAT THERE ARE ONLY 2 INTERFACES
            fid.write('%u %u\n'%(self.ny,self.nx) )
            fid.write('%f %f\n'%(self.dy_rad,self.dx_rad) )
            fid.write('%f %f\n'%(self.oy_rad,self.ox_rad) )
        else:
            fid=open(fnam,'a')
        for iy in range(self.ny):#Write the interface as a vector
            outs=''.join('{1:.{0}f}\n'.format(3,ii)  for ii in self.z[:,iy])
            fid.write(outs)
        fid.close()

def write_receivers_fmm(stalist,outfnam='tst_receivers.in'):
#Writes receivers.in file for the FMM code
#  This could be an attribute ouf the StationList class
#  Also, fmm receiver elevations should be negative for right-handed system
    fid = open(outfnam,'w')
    fid.write(str(len(stalist))+'\n')
    for sta in stalist:
        fid.write('{1:.{0}f} {2:.{0}f} {3:.{0}f}'.format(
                    4,sta.elev/1000*-1,sta.lat,sta.lon) ) #Replace 0.00 with elev
        fid.write('\n1\n1\n1\n')

def read_arrtimes(fnam='arrtimes.dat',ifplot=0): #default name
#Read arrival times from file fnam
    print 'Reading arrival times from '+fnam
    #Open and read header
    fid = open(fnam,'r')
    tmp=fid.readline().strip().split()
    nz=num(tmp[0]);nlat=num(tmp[1]);nlon=num(tmp[2]) #Number of grid points
    tmp=fid.readline().strip().split()
    dz=num(tmp[0]);dlat=num(tmp[1]);dlon=num(tmp[2]) #Spacing of grid points
    tmp=fid.readline().strip().split()
    oz=num(tmp[0]);olat=num(tmp[1]);olon=num(tmp[2]) #Origin of grid points
    tmp=fid.readline().strip().split()
    poo=num(tmp[0])                         #number of sets of arrival times
    tmp=fid.readline().strip().split()
    poo=num(tmp[0])                         #source and path for arrival time ???

    #Now read the traveltimes and sort into a matrix
    data=empty((nlat,nlon,nz))
    for ix in range(nlon):
        for iy in range(nlat):
            for iz in range(nz):
                tmp=float(fid.readline().strip())
                data[iy,ix,iz]=tmp
    fid.close()

    #Create vectors of geographic coordinates
    elon=olon+(dlon*nlon);
    elat=olat+(dlat*nlat);
    lonvec=linspace(olon,elon,nlon);
    latvec=linspace(olat,elat,nlat);

    #Now plot
    if ifplot: #adhoc
        flon,flat=load_faults()
        plt.figure(figsize=(6,6))
        plt.plot(flon,flat,'k')
        plt.axis([-118,-115,32,35])
        plt.pcolor(lonvec,latvec,data[:,:,6])
    return data

def find_containing_cube(px,py,pz,xvec,yvec,zvec):
#Find the 8 endpoints for the cell which contains point px,py
#  We take advantage of the regular grid
#  Assumes the point is inside the volume defined by xvec,yvec,zvec
#  Returns an array of size 8,3 where the rows contain x,y,z coordinates of the cubes endpoints
#  Also returns indexes of endpoints
    #Find the nearest node point and indexes <--"indices" sounds stupid to me
    xind,xnode=_find_nearest(px,xvec)
    yind,ynode=_find_nearest(py,yvec)
    zind,znode=_find_nearest(pz,zvec)
    #Now check if the 3 coordinates of p are greater or less than the node it is nearest
    if px>=xnode: #px is east of the nearest node
        xi=xind+1
        xn=xvec[xi]
    else:        #px is west of the nearest node
        xi=xind-1
        xn=xvec[xi]
    if py>=ynode: #px is north of the nearest node
        yi=yind+1
        yn=yvec[yi]
    else:        #px is south of the nearest node
        yi=yind-1
        yn=yvec[yi]
    if pz>=znode: #px is above the nearest node
        zi=zind+1
        zn=zvec[zi]
    else:        #px is below the nearest node
        zi=zind-1
        zn=zvec[zi]
    #Add new endpoints to define the cube
    endpoints=[]
    endpoints.append( [xnode,ynode,znode] )
    endpoints.append( [xn,ynode,znode] )
    endpoints.append( [xn,yn,znode]    )
    endpoints.append( [xnode,yn,znode] )
    endpoints.append( [xnode,ynode,zn] )
    endpoints.append( [xn,ynode,zn] )
    endpoints.append( [xn,yn,zn]    )
    endpoints.append( [xnode,yn,zn] )
    #Add indices
    indexes=[]
    indexes.append( [xind,yind,zind] )
    indexes.append( [xi,yind,zind] )
    indexes.append( [xi,yi,zind]    )
    indexes.append( [xind,yi,zind] )
    indexes.append( [xind,yind,zi] )
    indexes.append( [xi,yind,zi] )
    indexes.append( [xi,yi,zi]    )
    indexes.append( [xind,yi,zi] )

    return endpoints,indexes

def find_nearest(nparray,value):
    #Returns the nearest item in nparray to value
    idx= (abs(nparray-value)).argmin()
    return nparray.flat[idx]

def find_nearest_index(nparray,value):
    idx= (abs(nparray-value)).argmin()
    return idx

def _find_nearest(px,xvec):
#Find the nearest x in xvec 
#  returns index
    best_ind=0
    shortest=100000000.0;
    for ii in range(len(xvec)):
        if abs(xvec[ii]-px)<shortest:
            shortest=abs(xvec[ii]-px)
            best_ind=ii
    return best_ind,xvec[best_ind]

def read_binary_float(fid,n=0,precision='double'):
#read the nth float value from a binary file at the with the given precision
# following python conventions, 0 is the index of the first value
    if precision is 'single':
        numbytes=4;packstr='f'
    else:
        numbytes=8;packstr='d'
    offset=n*numbytes #Go to the right spot in the file
    fid.seek(offset)
    raw=fid.read(numbytes)
    tmp=struct.unpack(packstr,raw)
    val=tmp[0]
    return val

def grid_search_traveltimes_rms(arrsta,qx,qy,qz,arrvec,li):
#Find the minimum value of some criterion by performing a grid search
# We aren't necessarily searching the whole grid; we may be skipping values
#   sta         list of station names; strings
#   qx,qy,qz    vectors of indices to search through
#   arrvec      vector of arrivals in the same order as sta
#   li          Linear_index class for the entire traveltime grid
    from numpy import array #There is no reason this should be here, but the next line generated error messages if it wasn't. I'm confused
    rms=array([])
    search_inds=Linear_index(len(qx),len(qy),len(qz))
    for ix in qx:  #Loop over the three vectors, searching every point
        for iy in qy:
            for iz in qz:
                calctt=array([]); #initialize the calculated tt vector
                ind=li.get_1D(ix,iy,iz) #Find the vector index
                for sta in arrsta: #Build vector of calculated ttimes
                    #print 'bin.'+sta+'.traveltime'
                    fid = open('bin.'+sta+'.traveltime')
                    calctt=append(calctt, read_binary_float(fid,ind) )
                    #print read_binary_float(fid,ind)
                    fid.close()
                rms=append(rms, sqrt(mean( (arrvec-calctt)**2 )) ) #Root-mean-square, yo
    #Now find the 3D index of the best point so far
    min_ind=rms.argmin()
    (minx,miny,minz)=search_inds.get_3D(min_ind)
    minx=qx[minx]; miny=qy[miny]; minz=qz[minz];
    return minx,miny,minz

def grid_search_traveltimes_origin(arrsta,qx,qy,qz,arrvec,li):
#Find the minimum value of the origin time standard deviation following Ben-Zion et al., 1992 (JGR)
#  sta          list of station names; strings
#   qx,qy,qz    vectors of indices to search through
#   arrvec      vector of absolute arrivals in the same order as sta
#   li          Linear_index class for the entire traveltime grid
    from numpy import array #There is no reason this should be here, but the next line generated error messages if it wasn't. I'm confused
    origin_std=array([])
    origin_mean=array([])
    #origin_std=empty([ len(qy),len(qx),len(qz)] )+1000 #Give large starting values
    #origin_mean=empty([ len(qy),len(qx),len(qz)] )
    search_inds=Linear_index(len(qx),len(qy),len(qz))
    for ix in qx: #range(len(qx)):  #Loop over the three vectors, searching every point
        print 'On ix: ', ix,'/',len(qx),'\n'
        for iy in qy: #range(len(qy)):
            for iz in qz: #range(len(qz)):
                calctt=array([]); #initialize the calculated tt vector
                #ind=li.get_1D(qx[ix],qy[iy],qz[iz]) #Find the vector index
                ind=li.get_1D(ix,iy,iz)
                calctt=read_tt_vector(arrsta,ind) #Make a traveltime vector from calculated times
                orivec=arrvec-calctt #Take the difference
                origin_mean=append( origin_mean,orivec.mean()) #find the mean origin time
                if min(calctt)<0: #If the traveltime <0, this gridpoint is null
                    origin_std=append(origin_std,1000) #A large dummy value
                    continue
                origin_std=append(origin_std,orivec.std)
                #origin_std[iy,ix,iz]=orivec.std()
                #origin_mean[iy,ix,iz]=orivec.mean()
    #Now find the 3D index of the best point so far
    #new_origin=origin_mean.flatten()[oristd.argmin()] #The origin time with lowest std
    min_ind=origin_std.argmin()
    (minx,miny,minz)=search_inds.get_3D(min_ind)
    minx=qx[minx]; miny=qy[miny]; minz=qz[minz];
    return minx,miny,minz,origin_mean[min_ind]

def _grid_search_traveltimes_origin(arrsta,qx,qy,qz,arrvec,li):
#Find the minimum value of the origin time standard deviation following Ben-Zion et al., 1992 (JGR)
#  sta          list of station names; strings
#   qx,qy,qz    vectors of indices to search through
#   arrvec      vector of absolute arrivals in the same order as sta
#   li          Linear_index class for the entire traveltime grid
    from numpy import array #There is no reason this should be here, but the next line generated error messages if it wasn't. I'm confused
    rms=array([])
    #origin_std=array([])
    origin_std=empty([ len(qy),len(qx),len(qz)] )+1000 #Give large starting values
    origin_mean=empty([ len(qy),len(qx),len(qz)] )
    search_inds=Linear_index(len(qx),len(qy),len(qz))
    for ix in range(len(qx)):  #Loop over the three vectors, searching every point
        print 'On ix: ', ix,'/',len(qx),'\n'
        for iy in range(len(qy)):
            for iz in range(len(qz)):
                calctt=array([]); #initialize the calculated tt vector
                ind=li.get_1D(qx[ix],qy[iy],qz[iz]) #Find the vector index
                for sta in arrsta: #Build vector of calculated ttimes
                    fid = open('bin.'+sta+'.traveltime')
                    tmp=read_binary_float(fid,ind)
                    #if tmp<0: tmp=NaN #Kill negative travel times
                    calctt=append(calctt, tmp )
                    fid.close()
                orivec=arrvec-calctt
                origin_std[iy,ix,iz]=orivec.std()
                origin_mean[iy,ix,iz]=orivec.mean()
    #Now find the 3D index of the best point so far
    min_ind=origin_std.argmin()
    (minx,miny,minz)=search_inds.get_3D(min_ind)
    minx=qx[minx]; miny=qy[miny]; minz=qz[minz];
    return origin_std,origin_mean

def exp_grid_search(arrsta,qx,qy,qz,arrvec,li,ep):
#Find the minimum value of the origin time standard deviation following Ben-Zion et al., 1992 (JGR)
#  sta          list of station names; strings
#   qx,qy,qz    vectors of indices to search through
#   arrvec      vector of absolute arrivals in the same order as sta
#   li          Linear_index class for the entire traveltime grid
    from numpy import array #There is no reason this should be here, but the next line generated error messages if it wasn't. I'm confused
    rms=array([])
    origin_std=array([])
    origin_mean=array([])
    #origin_std=empty([ len(qy),len(qx),len(qz)] )+1000 #Give large starting values
    #origin_mean=empty([ len(qy),len(qx),len(qz)] )
    search_inds=Linear_index(len(qx),len(qy),len(qz))
    for ix in qx: #range(len(qx)):  #Loop over the three vectors, searching every point
        for iy in qy: #range(len(qy)):
            for iz in qz: #range(len(qz)):
                calctt=array([]); #initialize the calculated tt vector
                #ind=li.get_1D(qx[ix],qy[iy],qz[iz]) #Find the vector index
                ind=li.get_1D(ix,iy,iz)
                calctt=read_tt_vector(arrsta,ind) #Make a traveltime vector from calculated times
                orivec=arrvec-calctt #Take the difference
                origin_mean=append( origin_mean,ep-orivec.mean()) #find the mean origin time
                if min(calctt)<0: #If the traveltime <0, this gridpoint is null
                    origin_std=append(origin_std,1000) #A large dummy value
                    continue
                origin_std=append(origin_std,orivec.std)
                #origin_std[iy,ix,iz]=orivec.std()
                #origin_mean[iy,ix,iz]=orivec.mean()
    #Now find the 3D index of the best point so far
    #new_origin=origin_mean.flatten()[oristd.argmin()] #The origin time with lowest std
    min_ind=origin_mean.argmin()
    (minx,miny,minz)=search_inds.get_3D(min_ind)
    minx=qx[minx]; miny=qy[miny]; minz=qz[minz];
    return minx,miny,minz,origin_mean[min_ind]+ep

def locate_eq(ev):
    #Locate an earthquake based on the arrivals in ev, traveltime files which are already saved
    from antelope.stock import pfread, pfin
    from mtools import Phalist,Stalist
    params=pfin('eqloc3d.pf')
    loc_params = params['location_parameters']
    prop_params= params['propagation_grid']
    earth_rad=6371.0

    #Get Propagation grid paramters
    nlat=int(prop_params['nlat'])
    nlon=int(prop_params['nlon'])
    nz=int(prop_params['nr'])
    li = Linear_index(nlon, nlat, nz)
    olon=float(prop_params['minlon'])
    olat=float(prop_params['minlat'])
    oz=float(prop_params['minz'])
    dlon=float(prop_params['dlon'])
    dlat=float(prop_params['dlat'])
    dz=  float(prop_params['dr'])

    #Build vectors of geographic coordinates
    qlon = arange(olon, dlon * nlon + olon, dlon)
    qlat = arange(olat, dlat * nlat + olat, dlat)
    qdep = arange(earth_rad-dz*nz,earth_rad,dz)+oz+dz

    #Grid search for best location
    start_time=time.time()
    absvec=[]
    arrvec=[]
    arrsta=[]        #a list of station names
    for arrival in ev.arrivals:
        if arrival.phase is 'P':
            arrvec.append(arrival.time - ev.time)
            absvec.append(arrival.time)
            arrsta.append(arrival.sta)
        if not os.path.isfile(arrival.sta+'traveltime'):
            continue
    if len(arrvec)<6: #About this many phases are needed to get a decent result
        return None
    absvec=asarray(absvec)
    arrvec=asarray(arrvec)
    #print 'Number of phases used: ',len(arrvec)

    #Search coarsely
    #dstep should go in parameter file.
    dstep = int(loc_params['dstep2'])
    dx, dy, dz = nlon / dstep, nlat / dstep, nz / dstep
    #dx, dy, dz = 1,1,1 #Remove this later
    qx, qy, qz = range(1, nlon, dx), range(1, nlat, dy), range(1, nz, dz);
    minx, miny, minz, orgmin = grid_search_traveltimes_origin(arrsta, qx, qy,
                                                                qz, absvec, li)
    #minx, miny, minz, orgmin = exp_grid_search(arrsta,qx,qy,qz, absvec, li,ev.time)
    #Finer search
#    minx, miny, minz = grid_search_traveltimes_rms(arrsta, qx, qy, qz,
#                                                            arrvec,li)
    buff = int(loc_params['buff2'])
    qx = range(minx - buff, minx + buff)
    qy = range(miny - buff, miny + buff)
    qz = range(minz - buff, minz + buff);
    qx = fix_boundary_search(qx, li.nx)
    qy = fix_boundary_search(qy, li.ny)
    qz = fix_boundary_search(qz, li.nz)
    minx, miny, minz, orgmin = grid_search_traveltimes_origin(arrsta, qx, qy,
                                                                qz, absvec, li)
    #minx, miny, minz, orgmin = exp_grid_search(arrsta,qx,qy,qz, absvec, li,ev.time)
#    minx, miny, minz = grid_search_traveltimes_rms(arrsta, qx, qy, qz,
#                                                            arrvec,li)
    orgmin=0

    #Find the best subgrid location
    c,resid=get_subgrid_loc(minx,miny,minz,arrvec,arrsta,li)
    delta_x=qlon[1]-qlon[0]
    delta_y=qlat[1]-qlat[0]
    delta_z=qdep[1]-qdep[0]
    loc_change=c*[delta_x,delta_y,delta_z] #Subgrid location change in lon/lat/depth

    #Find the best-fit source location in geographic coordinates
    newloc=[newlon,newlat,newz]=[ qlon[minx],qlat[miny],qdep[minz] ] #+loc_change
    elapsed_time=time.time()-start_time
    print ev.orid,len(arrvec),newlon,newlat,earth_rad-newz,ev.lon,ev.lat,ev.depth,ev.time-orgmin,elapsed_time,resid
    #This function will return location, origin time, and error estimates for both

def get_subgrid_loc(ix,iy,iz,arrvec,arrsta,li):
    #Test least squares on real data
    import numpy as np
    from scipy import linalg
    import matplotlib.pyplot as plt

    #Get traveltime vectors for the closest point and its neighbors
    ind=li.get_1D(ix,iy,iz)
    tt000= read_tt_vector(arrsta,ind)
    ind=li.get_1D(ix+1,iy,iz)
    tt100= read_tt_vector(arrsta,ind)
    ind=li.get_1D(ix,iy+1,iz)
    tt010= read_tt_vector(arrsta,ind)
    ind=li.get_1D(ix,iy,iz+1)
    tt001= read_tt_vector(arrsta,ind)
    #backwards
    ind=li.get_1D(ix-1,iy,iz)
    btt100= read_tt_vector(arrsta,ind)
    ind=li.get_1D(ix,iy-1,iz)
    btt010= read_tt_vector(arrsta,ind)
    ind=li.get_1D(ix,iy,iz-1)
    btt001= read_tt_vector(arrsta,ind)
    #Calculate forward derivatives
    dt_dx=tt100-tt000
    dt_dy=tt010-tt000
    dt_dz=tt001-tt000
    #backwards
    bdt_dx=tt000-btt100
    bdt_dy=tt000-btt010
    bdt_dz=tt000-btt001
    #Build A matrix and r in r=Ax  (x is the spatial vector here [x,y,z]
    A=c_[dt_dx,dt_dy,dt_dz]
    r=arrvec-tt000
    c,resid,rank,sigma=linalg.lstsq(A,r)
    #backwards
    bA=c_[bdt_dx,bdt_dy,bdt_dz]
    bc,resid,rank,sigma=linalg.lstsq(bA,r)
    return c,resid

def read_tt_vector(stanames,ind):
    #Read from binary files to create a vector of traveltimes for the stations in stanames
    #   at the index location ind, which is the 1D index
    #   stanames is a list of station names only
    from numpy import array
    ttvec=array([])
    for sta in stanames:
        fid=open('bin.'+sta+'.traveltime')
        ttvec=append(ttvec, read_binary_float(fid,ind) )
        fid.close()
    return ttvec

def fix_boundary_search(qx,nx):
#When performing a grid search on a subgrid, make sure you don't go off the edges
#  qx         search vectors, these will be modified then returned
#  nx         max index [li.nx]
    for ix in range(len(qx)):
        if qx[ix]<0:
            qx[ix]=0
        if qx[ix]>=nx:
            qx[ix]=nx-1
    newqx=uniq(qx)
    return newqx

def uniq(input):
#Remove duplicate items from a list. Preserves order.
  output = []
  for x in input:
    if x not in output:
      output.append(x)
  return output

class Linear_index():
    #Holds a 1D list of 3D indices and a 3D list of 1D indices
    # where iz varies fastest, then iy, then ix
    #  The speed of this can certainly be improved
    def __init__(self,nx,ny,nz):
        from numpy import empty
        self.nx=nx; self.ny=ny; self.nz=nz;
        self.i1D=[]
        self.i3D=empty((nx,ny,nz)) #python has weird index conventions
        ic=0
        for ix in range(nx): #Fuck it, just be explicit
            for iy in range(ny):
                for iz in range(nz):
                    self.i1D.append( (ix,iy,iz) )
                    self.i3D[ix,iy,iz]=ic
                    ic=ic+1
        self.i3D=self.i3D.astype(int)
    def get_1D(self,ix,iy,iz):
        return self.i3D[ix,iy,iz]
    def get_3D(self,iv):
        return self.i1D[iv]

class _Linear_index():
    #A class which can quickly convert 3D indices to a single vector index
    def __init__(self,nx,ny,nz):
        #We are using the FMM fast order which is x, y, z (lon,lat,depth)
        self.nx=nx; self.ny=ny; self.nz=nz

    def get_1D(self,ix,iy,iz):
        #Convert 3 indices to a single vector index
        #CAREFUL WITH PYTHON INDEXING
        iv = ix*self.ny*self.nz + iy*self.nz + iz #Check on this
        return iv

    def get_3D(self,iv):
        #Convert the vector index into 3 matrix indices
        ix= iv/(self.ny*self.nz)
        iy=(iv%(self.ny*self.nz))/self.nz
        iz=iv-ix*self.ny*self.nz-iy*self.nz
        return ix,iy,iz

class Parameters(): #THIS CAN PROBABLY JUST BE DELETED. I DON'T THINK I EVER USED IT!
#The basic parameters of the velocity and interface grids and whatnot
    def __init__(self):
        #Velocity grid
        self.vel_minlon=     -118.007
        self.vel_dlon=       0.0109
        self.vel_nlon=       273
        self.vel_minlat=     32.4613
        self.vel_dlat=       0.0091
        self.vel_nlat=       252
        self.vel_minr=       6325.0
        self.vel_dr=         1.0
        self.vel_nr=         51
        #derived
        self.vel_dlon_rad=  self.vel_dlon*pi/180.0 #These are the actual values in vgrids.in
        self.vel_dlat_rad=  self.vel_dlat*pi/180.0

        #Propagation grid !!!Later, I'll derive all these from the values for the velocity grid above
        self.prop_minlon=        -117.95
        self.prop_minlat=        32.5
        self.prop_minz=          0.0
        self.prop_nlon=          250
        self.prop_nlat=          230
        self.prop_nr=            46
        #derived
        self.prop_dlon=          self.vel_dlon
        self.prop_dlat=          self.vel_dlat
        self.prop_dr=            self.vel_dr
        #Build vectors of the grid axes. These loops are more accurate than linspace somehow
        self.prop_lonvec=[]
        self.prop_latvec=[]
        self.prop_rvec=[]
        for ii in range(self.prop_nlon):
            self.prop_lonvec.append( self.vel_minlon+self.vel_dlon*ii )
        for ii in range(self.prop_nlat):
            self.prop_latvec.append( self.vel_minlat+self.vel_dlat*ii )
        for ii in range(self.prop_nr):
            self.prop_rvec.append(   self.vel_minr+self.vel_dr*ii     )
#self.prop_rvec=linspace(self.vel_minr,self.vel_minr+(self.vel_dr*self.vel_nr),self.vel_nr)

class Station():
    """
    A container class to hold station metadata.
    """
    def __init__(self, name, lat, lon, elev):
        self.name, self.lat, self.lon, self.elev = name, lat, lon, elev

    def __str__(self):
        ret = 'Station Object\n--------------\n'
        ret += 'name:\t\t%s\n' % self.name
        ret += 'lat: \t\t%f\n' % self.lat
        ret += 'lon: \t\t%f\n' % self.lon
        ret += 'elev: \t\t%f\n' % self.elev

    def show(self): #show the contents of the class
        print '%s %5.4f %5.4f %5.2f'% (self.name,self.lon,self.lat,self.elev)

class StationList(list):
    """
    A container class for a list of Station objects.

    This class replaces the deprecated Stalist class.
    """
    def __init__(self, inp, is_db=True):
        if is_db: self._init_db(inp)
        else: self._init_scedc(inp)

    def __str__(self):
        ret = 'StationList Object\n------------------\n'
        for sta in self:
            ret += 'sta:\t\t%s\n' % sta.name
            ret += 'lat:\t\t%f\n' % sta.lat
            ret += 'lon:\t\t%f\n' % sta.lon
            ret += 'elev:\t\t%f\n' % sta.elev
        return ret

    def _init_db(self, db):
        """
        Initialize station list using a CSS3.0 database as input.
        """
        with closing(dbopen(db, 'r')) as db:
            tbl_site = db.schema_tables['site']
#The following line will be taken out
            tbl_site = tbl_site.subset('lon >= -117.80 && lat >= 32.5 && lon <= -115.4456 && lat <= 34.5475')
            tbl_site = tbl_site.sort('sta', unique=True)
            for record in tbl_site.iter_record():
                sta, lat, lon, elev = record.getv('sta', 'lat', 'lon', 'elev')
                self.append(Station(sta, lat, lon, elev))

    def _init_scedc(self, infile):
        """
        Initialize station list using SCEDC format flat file as input.
        """
        infile = open(infile, 'r')
        for line in infile:
            line = line.strip().split() #stripping may be uneccessary
            self.append(Station(line[0],
                                 float(line[1]),
                                 float(line[2]),
                                 float(line[3])))

class _Event():
    """
    A container class for earthquake event metadata.
    """
    #def __init__(self, time, lat, lon, depth, mag, magtype=None, evid=None):
    def __init__(self, *args, **kwargs):
        """
        Initialize Event object using one of two possible inputs.

        This first (standard) input method is to simply supply a 
        CSS3.0 database path and an evid.

        The second method requires an auxiliary read function which will parse
        an input flat file (probably SCEDC format) and pass the following
        arguments in the following order.
        EG. event = Event(time,
                          latitude,
                          longitude,
                          depth,
                          magnitude,
                          phase_list,
                          magtype=magnitude_type,
                          evid=event_id)
        The keyword arguments (magtype and evid) are optional.
        """
        if len(args) == 2: self._init_from_db(*args)
        else: self._init_from_aux(*args, **kwargs)

    def __str__(self):
        ret = 'Event Object\n------------\n'
        ret += 'evid:\t\t%d\n' % self.evid
        ret += 'time:\t\t%f\n' % self.time
        ret += 'lat:\t\t%f\n' % self.lat
        ret += 'lon:\t\t%f\n' % self.lon
        ret += 'depth:\t\t%f\n' % self.depth
        ret += 'mag:\t\t%f\n' % self.mag
        ret += 'magtype:\t%s\n' % self.magtype
        ret += 'year:\t\t%d\n' % self.year
        ret += 'month:\t\t%d\n' % self.month
        ret += 'day:\t\t%d\n' % self.day
        ret += 'hour:\t\t%d\n' % self.hour
        ret += 'minute:\t\t%d\n' % self.minute
        ret += 'second:\t\t%d\n' % self.second
        ret += 'arrivals:\n'
        if len(self.arrivals) == 0:
            ret += '\t\tNone\n'
        else:
            for i in range(len(self.arrivals)):
                for line in  ('%s' % self.arrivals[i]).split('\n'):
                    ret += '\t\t%s\n' % line
        return ret

    def _init_from_db(self, db, evid):
        """
        Initialize Event object using a CSS3.0 database as input.
        """
        if evid == None: raise(Exception('No \'evid\' supplied. Could not '
            'initialize Event object from CSS3.0 database.'))
        with closing(dbopen(db, 'r')) as db:
            view = db.schema_tables['event']
            view = view.join('origin')
            view = view.subset('evid == %s' % evid)
            view = view.subset('orid == prefor')
            #If for some reason this subset is empty, just take the first
            #solution as preferred. EG. prefor field is unitialized.
            if view.record_count == 0:
                view = db.schema_tables['origin']
                view = db.schema_tables['event']
                view = view.join('origin')
                view = view.subset('evid == %s' % evid)
            view = view.join('netmag', outer=True)
            view.record = 0
            evid, time, lat, lon, depth, mag, magtype =  view.getv('evid',
                'time', 'lat', 'lon', 'depth', 'magnitude', 'magtype')
            self.evid       = evid
            self.time       = time
            self.lat        = lat
            self.lon        = lon
            self.depth      = depth
            self.mag        = mag
            self.magtype    = magtype
            self.year       = int(epoch2str(time, '%Y'))
            self.month      = int(epoch2str(time, '%m'))
            self.day        = int(epoch2str(time, '%d'))
            self.hour       = int(epoch2str(time, '%H'))
            self.minute     = int(epoch2str(time, '%M'))
            self.second     = float(epoch2str(time, '%S.%s'))
            view = view.join('assoc')
            view = view.join('arrival')
            arrivals = [ record.getv('sta',
                                     'arrival.time',
                                     'phase') \
                                     + (None, ) \
                                     for record in view.iter_record()
                       ]
            self.arrivals = [ Phase(sta, time, phase, qual)
                            for sta, time, phase, qual in arrivals
                            ]

    def _init_from_aux(self, time, lat, lon, dpeth, mag, magtype=None, evid=None):
        """
        Initialize Event object using an auxiliary read function.

        Auxiliary read function must parse flat file and pass parameters 
        to Event contsructor.
        """
        self.time = time
        self.lat        = lat
        self.lon        = lon
        self.depth      = depth
        self.mag        = mag
        self.magtype    = magtype
        self.year = int(epoch2str(time, '%Y'))
        self.month = int(epoch2str(time, '%m'))
        self.day = int(epoch2str(time, '%d'))
        self.hour = int(epoch2str(time, '%H'))
        self.minute = int(epoch2str(time, '%M'))
        self.second = float(epoch2str(time, '%S.%s'))
        self.arrivals = phase_list

class Phase():
    """
    A container class for phase data.
    """
    def __init__(self, sta, time, phase, qual):
        self.sta = sta
        self.time = time
        self.phase = phase
        self.qual = qual

    def __str__(self):
        ret = 'Arrival Object\n--------------\n'
        ret += 'sta:\t\t%s\n' % self.sta
        ret += 'time:\t\t%f\n' % self.time
        ret += 'phase:\t\t%s\n' % self.phase
        ret += 'qual:\t\t%s\n'  % self.qual
        return ret

class Traveltime_header_file():
#This reads the header file passed as an argument by fnam, saves all of the header info to a class
    def __init__(self,fnam):
       fid = open(fnam,'r')
       tmp=fid.readline().strip().split()
       self.nz=num(tmp[0]);self.nlat=num(tmp[1]);self.nlon=num(tmp[2]) #Number of grid points
       tmp=fid.readline().strip().split()
       self.dz=num(tmp[0]);self.dlat=num(tmp[1]);self.dlon=num(tmp[2]) #Spacing of grid points
       tmp=fid.readline().strip().split()
       self.oz=num(tmp[0]);self.olat=num(tmp[1]);self.olon=num(tmp[2]) #Origin of grid points
       tmp=fid.readline().strip().split()
       self.arrsets=num(tmp[0])                         #number of sets of arrival times
       tmp=fid.readline().strip().split()
       self.srcpath=num(tmp[0])                         #source and path for arrival time ???
       fid.close()

def process_events(dbin, dbout, events, station_list, params):
    from run_3D_loc import run_3d_location_algorithm
    for event in events:
        event = _Event(dbin, event)
        solution = run_3d_location_algorithm(event, station_list, params)

def create_event_list(inp, fmt):
    """
    Create and return a list of Event objects.

    Arguments:
    inp - A (potentially subsetted) View of the Event table from the
    CSS3.0 database schema. Alternatively the path to a SCEDC format
    flat file.

    fmt - The format of the 'inp' argument; 'CSS3.0' or 'SCEDC'.

    Return Values:
    A list of Event objects.
    """
    from eqloc3d_classes import Event
    import time as pytime
    import calendar
    event_list = []
    if fmt == 'CSS3.0':
        for record1 in inp.iter_record():
            evid, evname, prefor, auth, commid, lddate = record1.getv('evid',
                                                                     'evname',
                                                                     'prefor',
                                                                     'auth',
                                                                     'commid',
                                                                     'lddate')
            event = Event(evid,
                          prefor,
                          evname=evname,
                          auth=auth,
                          commid=commid,
                          lddate=lddate)
            view2 = inp.subset('evid == %d' % evid)
            view2 = view2.join('origin')
            view2 = view2.separate('origin')
            for record2 in view2.iter_record():
                lat = record2.getv('lat')[0]
                lon = record2.getv('lon')[0]
                depth = record2.getv('depth')[0]
                time = record2.getv('time')[0]
                orid = record2.getv('orid')[0]
                evid = record2.getv('evid')[0]
                jdate = record2.getv('jdate')[0]
                nass = record2.getv('nass')[0]
                ndef = record2.getv('ndef')[0]
                ndp = record2.getv('ndp')[0]
                grn = record2.getv('grn')[0]
                srn = record2.getv('srn')[0]
                etype = record2.getv('etype')[0]
                review = record2.getv('review')[0]
                depdp = record2.getv('depdp')[0]
                dtype = record2.getv('dtype')[0]
                mb = record2.getv('mb')[0]
                mbid = record2.getv('mbid')[0]
                ms = record2.getv('ms')[0]
                msid = record2.getv('msid')[0]
                ml = record2.getv('ml')[0]
                mlid = record2.getv('mlid')[0]
                algorithm = record2.getv('algorithm')[0]
                auth = record2.getv('auth')[0]
                commid = record2.getv('commid')[0]
                lddate = record2.getv('lddate')[0]
                view3 = view2.subset('orid == %d' % orid)
                view3 = view3.join('assoc')
                view3 = view3.join('arrival')
                arrival_data = [record3.getv('sta',
                                             'arrival.time',
                                             'iphase')\
                                             + (None, )\
                                             for record3 in view3.iter_record()]
                arrivals = [Phase(sta, time, phase, qual)
                            for sta, time, phase, qual in arrival_data]
                event.add_origin(lat,
                                 lon,
                                 depth,
                                 time,
                                 orid,
                                 evid,
                                 auth,
                                 arrivals,
                                 jdate=jdate,
                                 nass=nass,
                                 ndef=ndef,
                                 ndp=ndp,
                                 grn=grn,
                                 srn=srn,
                                 etype=etype,
                                 review=review,
                                 depdp=depdp,
                                 dtype=dtype,
                                 mb=mb,
                                 mbid=mbid,
                                 ms=ms,
                                 msid=msid,
                                 ml=ml,
                                 mlid=mlid,
                                 algorithm=algorithm,
                                 commid=commid,
                                 lddate=lddate)
            event.set_preferred_origin(event.prefor)
            event_list += [event]
    elif fmt == 'SCEDC':
#evid should be set using an idserver!!!
        evid_ctr = 100000
        infile = open(inp, 'r')
        event = None
        for line in infile:
            line = line.strip().split()
            #In SCEDC format, event lines begin with '#'
            if line[0] == '#':
                if event != None:
                    event_list += [event]
                year = int(line[1])
                month = int(line[2])
                day = int(line[3])
                hour = int(line[4])
                minute = int(line[5])
                second = float(line[6])
                lat = float(line[7])
                lon = float(line[8])
                depth = float(line[9])
                mag = float(line[10])
                orid = int(line[14])
                time = calendar.timegm(pytime.struct_time(list([year,
                                              month,
                                              day,
                                              hour,
                                              minute,
                                              second,
                                              0, 0, 0])))
                event = Event(evid_ctr,
                              orid,
                              auth='SCEDC',
                              lddate=pytime.time())
                event.add_origin(lat, lon, depth, time, orid, evid_ctr, 'SCEDC')
                event.set_preferred_origin(orid)
                evid_ctr += 1
            else:
                sta = line[0]
                arrtime = float(line[1]) + time
                qual = float(line[2])
                iphase = line[3]
                event.preferred_origin.arrivals += [Phase(sta,
                                                          arrtime,
                                                          iphase,
                                                          qual)]

    else:
        raise Exception('Input format %s not recognized.' % fmt)
    return event_list
