import subprocess
import time
import struct
from numpy import *
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

def run_fmm():
#Run the traveltime solver with all of the input files currently in the directory
#  Note: Input files should already be set using write_sources_in, etc...
    subprocess.call('./fm3d',shell=True)

def write_sources_in(src_dep,src_lat,src_lon):
#Write  the sources.in file, which has this form:
    # 1                                number of sources
    # 0                                source is local/teleseismic (0/1)
    # 5.00  33.0   -116.0      position depth(km),lat(deg),long(deg)
    # 1                                number of paths from this source
    # 1                                number of sections on the path
    # 0 1           define the path sections
    # 1            define the velocity type along the path
    numsrc=1
    teleflag=0
    numpaths=1
    numsections=1
    veltype=1
    fid = open('sources.in','w')
    fid.write(" %i\n %i\n %.4f %.4f %.4f\n %i\n %i\n %i %i\n %i\n" % (numsrc,teleflag,src_dep,src_lat, src_lon, numpaths, numsections, 0,1,veltype))
    fid.close()

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
    print 'poop'

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

def gen_sta_tt_maps(stalist,if_write_binary=True):
#Generate a travel time map for each station in the station list
    start_time=time.time()
    print 'Starting travel time calculation at ',start_time,'\n'
    for sta in stalist:
        print 'Generating travel times for station', sta.name, '\n'
        write_sources_in(.10,sta.lat,sta.lon) # !!!! For now, depth is set to 0
        run_fmm()
        #Create output file name
        outfnam=sta.name+'.traveltime'
        subprocess.call('mv arrtimes.dat '+outfnam,shell=True)
        if if_write_binary:
            tt_ascii_to_binary(outfnam)
    elapsed_time=time.time()-start_time
    print 'Finished travel time calculations at ',time.time(),'\n'
    print 'Total calculation time: %8.4f seconds.\n' % (elapsed_time)

def find_containing_cube(px,py,pz,xvec,yvec,zvec):
#Find the 8 endpoints for the cell which contains point px,py
#  We take advantage of the regular grid
#  Assumes the point is inside the volume defined by xvec,yvec,zvec
#  Returns an array of size 8,3 where the rows contain x,y,z coordinates of the cubes endpoints
#  Also returns indexes of endpoints
    #Find the nearest node point and indexes <--"indices" sounds stupid to me
    xind,xnode=find_nearest(px,xvec)
    yind,ynode=find_nearest(py,yvec)
    zind,znode=find_nearest(pz,zvec)
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

def find_nearest(px,xvec):
#Find the nearest x in xvec 
#  returns index
    best_ind=0
    shortest=100000000.0;
    for ii in range(len(xvec)):
        if abs(xvec[ii]-px)<shortest:
            shortest=abs(xvec[ii]-px)
            best_ind=ii
    return best_ind,xvec[best_ind]

def tt_ascii_to_binary(fnam):
#Convert an ascii traveltime file to binary format
# Also puts the header in a separate ascii file
# just put "bin" and "hdr" in front of the filename
    binfnam='bin.'+fnam
    hdrfnam='hdr.'+fnam
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
    narr=num(tmp[0])                         #number of sets of arrival times
    tmp=fid.readline().strip().split()
    null1=num(tmp[0]);null2=num(tmp[1]);null3=num(tmp[2])                         #source and path for arrival time ???
    #Now read the traveltimes into a list
    lines=fid.readlines()
    data=[]
    for line in lines:
        data.append( float(line) )
    fid.close()
    # Now output in binary format
    fid=open(binfnam,'w')
    out_array=array('d',data)#It is apparently faster to turn this into an array
    out_array.tofile(fid)
    fid.close()

    #Now Write the header
    fid = open(hdrfnam,'w')
    fid.write(" %i %i %i\n %.4f %.6f %.6f\n %.5f %.5f %.5f\n %i\n %i %i %i" % (nz,nlat,nlon, dz,dlat,dlon, oz,olat,olon, narr, null1,null2,null3))
    fid.close()
    print 'Finished writing arrival times to '+binfnam+' and '+hdrfnam
    #return data #Don't forget to remove this!!

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

def grid_search_traveltimes(arrsta,qx,qy,qz,arrvec,li):
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

def fix_boundary_search(qx,nx):
#When performing a grid search on a subgrid, make sure you don't go off the edges
#  qx         search vectors, these will be modified then returned
#  nx         max index [li.nx]
    for ix in range(len(qx)):
        if qx[ix]<0:
            qx[ix]=0
        if qx[ix]>nx:
            qx[ix]=nx
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

class Parameters():
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

class Stalist(list):
#A class containing a list of stations and their locations
    def __init__(self,fnam):
        #Read the station file
        fid = open(fnam,'r')
        a = fid.readlines()
        #Parse into a structure containing name,lon,lat,elevation
        for line in a:
            tmp = line.strip().split()
            self.append(self.Sta(tmp[0],float(tmp[1]),float(tmp[2]),float(tmp[3])))
        fid.close()

    class Sta():
    #a class for each individual station. There is almost certainly a more elegant way to do this.
        def __init__(self,name,lat,lon,elev):
            self.name=name; self.lat=lat; self.lon=lon; self.elev=elev

class Phalist(list):
#A class containing the event information with arrival times for each event
    def __init__(self,fnam):
        #Read the file
        print 'Reading ' +fnam
        fid = open(fnam,'r')
        a = fid.readlines()
        #Parse into a structure containing events
        #self=list()
        ic = -1
        for line in a:
            tmp=line.strip().split()
            if tmp[0] is '#':
                ic=ic+1
                self.append( {'id':int(tmp[14]),'year':int(tmp[1]),'month':int(tmp[2]),'day':int(tmp[3]),'hour':int(tmp[4]),'min':int(tmp[5]),'sec':float(tmp[6]),'lon':float(tmp[8]),'lat':float(tmp[7]),'depth':float(tmp[9]),'mag':float(tmp[10]),'arrivals':list()} )
            else:
                self[ic]['arrivals'].append( {'staname':tmp[0],'ttime':float(tmp[1]),'qual':float(tmp[2]),'phase':tmp[3]} )
        print 'Finished reading '+fnam

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
