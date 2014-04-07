#Test least squares on real data
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from misc_tools import *
from mtools import Phalist,Stalist

earth_rad=6371

#Load events
fnam='pha.phase_220to9_july13.dat'
pha=Phalist(fnam,'scedc')

#Load stations
fnam='sta.phase_220to9_july13.dat'
sta=Stalist(fnam)

#Read the header for the last traveltime file
fnam = '%s.traveltime' % sta[0].name
hdr=Traveltime_header_file(fnam)
nlat=hdr.nlat
nlon=hdr.nlon
nr=hdr.nz
nx, ny, nz = nlon, nlat, nr
li = Linear_index(nlon, nlat, nz)

#Build vectors of geographic coordinates
qlon = arange(hdr.olon, hdr.dlon * hdr.nlon + hdr.olon, hdr.dlon)
qlat = arange(hdr.olat, hdr.dlat * hdr.nlat + hdr.olat, hdr.dlat)
qdep = earth_rad-arange(hdr.oz, hdr.dz * hdr.nz + hdr.oz, hdr.dz)
delta_x=qlon[1]-qlon[0]
delta_y=qlat[1]-qlat[0]
delta_z=qdep[1]-qdep[0]

for ev in pha:
    #Find the closest indices to the known ev location
    ix=nonzero(qlon==find_nearest(qlon,ev.lon))[0][0]
    iy=nonzero(qlat==find_nearest(qlat,ev.lat))[0][0]
    iz=nonzero(qdep==find_nearest(qdep,ev.depth))[0][0]
    iz=iz
    #Get station names for arrivals and a traveltime vector
    absvec=[]
    arrvec=[]
    arrsta=[]        #a list of station names
    for arrival in ev.arrivals:
        if arrival.phase is 'P':
            arrvec.append(arrival.ttime)
            absvec.append(arrival.epoch)
            arrsta.append(arrival.staname)
        if not os.path.isfile(arrival.staname+'traveltime'):
            continue
    absvec=asarray(absvec)
    arrvec=asarray(arrvec)

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

    #Calculate updated locations
    loc_change=c*[delta_x,delta_y,delta_z]
    newloc=[newlon,newlat,newz]=[ qlon[ix],qlat[iy],qdep[iz] ]+loc_change
    #Backwards
    bloc_change=bc*[delta_x,delta_y,delta_z]
    bnewloc=[bnewlon,bnewlat,bnewz]=[ qlon[ix],qlat[iy],qdep[iz] ]+bloc_change

    #Forward derivative
    print '---------------\nForward Derivative:'
    print 'Location changed by: ',loc_change,resid
    print 'Initial Location: ',[qlon[ix],qlat[iy],qdep[iz]]
    print 'New Location: ',newloc

    #Backward derivative
    print '---------------\nBackward Derivative:'
    print 'Location changed by: ',bloc_change,resid
    print 'Initial Location: ',[qlon[ix],qlat[iy],qdep[iz]]
    print 'New Location: ',bnewloc
    print '---------------------------'
