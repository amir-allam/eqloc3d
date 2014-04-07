"""
Test the 3-D location algorithm.
"""
from misc_tools import *
from mtools import *
from numpy import arange,asarray
from antelope.stock import pfread, pfin
params=pfin('eqloc3d.pf')
loc_params = params['location_parameters']
#These should go in parameter file.
#nr = int(loc_params['nr'])
#nlat = int(loc_params['nlat'])
#nlon = int(loc_params['nlon'])
#nx, ny, nz = nlon, nlat, nr
earth_rad=6371

#Load events
fnam='pha.phase_210to9_july2011.dat'
pha=Phalist(fnam,'scedc')

#Load stations
fnam='sta.phase_210to9_july2011.dat'
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
qdep = arange(hdr.oz, hdr.dz * hdr.nz + hdr.oz, hdr.dz)

for ev in pha:
    #Grid search for best location
    #    arrvec=array([]) #a vector of travel times
    #    absvec=array([]) #absolute arrival times
    start_time=time.time()
#    arrvec = asarray([arrival.ttime for arrival in ev.arrivals]) #Use numpy
#    absvec = asarray([arrival.epoch for arrival in ev.arrivals])
#    arrsta = [arrival.staname for arrival in ev.arrivals] #List of strings
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
    print 'Number of phases used: ',len(arrvec)

    #Search coarsely
    #dstep should go in parameter file.
    dstep = int(loc_params['dstep1'])
    dx, dy, dz = nlon/dstep, nlat/dstep, nr/dstep;
    qx, qy, qz = range(1, nlon, dx), range(1, nlat, dy), range(1, nr, dz)
    minx, miny, minz = grid_search_traveltimes_rms(arrsta, qx, qy, qz,
                                                    arrvec, li)

    #Finer search
    #buff should go in paramter file.
    buff = int(loc_params['buff1'])
    qx = range(minx - buff, minx + buff)
    qy = range(miny - buff, miny + buff)
    qz = range(minz - buff, minz + buff);
    qx = fix_boundary_search(qx, li.nx)
    qy = fix_boundary_search(qy, li.ny)
    qz = fix_boundary_search(qz, li.nz)
    minx, miny, minz = grid_search_traveltimes_rms(arrsta, qx, qy, qz,
                                                    arrvec,li)
    #Find the best-fit source location in geographic coordinates
    lon1, lat1, z1 = qlon[minx], qlat[miny], qdep[minz]

    #Search coarsely
    #dstep should go in parameter file.
    dstep = int(loc_params['dstep2'])
    dx, dy, dz = nlon / dstep, nlat / dstep, nr / dstep
    qx, qy, qz = range(1, nlon, dx), range(1, nlat, dy), range(1, nr, dz);
    minx, miny, minz, orgmin = grid_search_traveltimes_origin(arrsta, qx, qy,
                                                                qz, absvec, li)
    #Finer search
    buff = int(loc_params['buff2'])
    qx = range(minx - buff, minx + buff)
    qy = range(miny - buff, miny + buff)
    qz = range(minz - buff, minz + buff);
    qx = fix_boundary_search(qx, li.nx)
    qy = fix_boundary_search(qy, li.ny)
    qz = fix_boundary_search(qz, li.nz)
    minx, miny, minz, orgmin = grid_search_traveltimes_origin(arrsta, qx, qy,
                                                                qz, absvec, li)
    #Find the best-fit source location in geographic coordinates
    lon2, lat2, z2 = qlon[minx], qlat[miny], qdep[minz]
    elapsed_time=time.time()-start_time
    print 'RMS: ',lon1,lat1,earth_rad-z1,'|| Origin: ',lon2,lat2,earth_rad-z2,' || ',ev.lon,ev.lat,ev.depth,'|| Time change:',ev.epoch-orgmin
    print 'Elapsed time (s): ', elapsed_time
