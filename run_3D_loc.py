def run_3d_location_algorithm(event, station_list, params):
    """
    Perform 3-D location algorithm.
    """
    from misc_tools import Linear_index,\
                           Traveltime_header_file,\
                           grid_search_traveltimes_rms,\
                           fix_boundary_search,\
                           grid_search_traveltimes_origin
    from numpy import arange
    loc_params = params['location_parameters']

    #Return None if there are no arrival for events.
    #This should probably be recorded in log.
    print params
    if len(event.arrivals) == 0: return None

    #These should go in parameter file.
    nr = loc_params['nr']
    nlat = loc_params['nlat']
    nlon = loc_params['nlon']
    nx, ny, nz = nlon, nlat, nr
    li = Linear_index(nlon, nlat, nr)

    #Read the header for the last traveltime file
#    fnam=evpha[0]['arrivals'][0]['staname']+'.traveltime'
    fnam = '%s.traveltime' % event.arrivals[0].sta
    hdr=Traveltime_header_file(fnam)
    #Build vectors of geographic coordinates
    qlon = arange(hdr.olon, hdr.dlon * hdr.nlon + hdr.olon, hdr.dlon)
    qlat = arange(hdr.olat, hdr.dlat * hdr.nlat + hdr.olat, hdr.dlat)
    qdep = arange(hdr.oz, hdr.dz * hdr.nz + hdr.oz, hdr.dz)

    #Grid search for best location
#    arrvec=array([]) #a vector of travel times
#    absvec=array([]) #absolute arrival times
#    arrsta=[]        #a list of station names
    absvec = [arrival.time for arrival in event.arrivals]
    arrsta = [arrival.sta for arrival in event.arrivals]
#Deprecate the use of travel-times, the loop below should not be used.
#    for arrival in event.arrivals: #Loop over arrivals and make vectors
#        if not os.path.isfile(arrival['staname']+'traveltime'):
#            continue
#        arrvec=append(arrvec, float(arrival['ttime']) ) #Build vector of observed ttimes

    #Search coarsely
    #dstep should go in parameter file.
    dstep = loc_params['dstep1']
    dx, dy, dz = nlon/dstep, nlat/dstep, nr/dstep;
    qx, qy, qz = range(1, nlon, dx), range(1, nlat, dy), range(1, nr, dz)
    minx, miny, minz = grid_search_traveltimes_rms(arrsta, qx, qy, qz,
                                                   arrvec, li)

    #Finer search
    #buff should go in paramter file.
    buff = loc_params['buff1']
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
    dstep = loc_params['dstep2']
    dx, dy, dz = nlon / dstep, nlat / dstep, nr / dstep
    qx, qy, qz = range(1, nlon, dx), range(1, nlat, dy), range(1, nr, dz);
    minx, miny, minz, orgmin = grid_search_traveltimes_origin(arrsta, qx, qy,
                                                              qz, absvec, li)
    #Finer search
    #buff should go in parameter file.
    buff = loc_params['buff2']
    qx = range(minx - buff, minx + buff)
    qy = range(miny - buff, miny + buff)
    qz = range(minz - buff, minz + buff);
    qx = fix_boundary_search(qx, li.nx)
    qy = fix_boundary_search(qy, li.ny)
    qz = fix_boundary_search(qz, li.nz)
    minx, miny, minz, orgmin = grid_search_traveltimes_origin(arrsta, qx, qy,
                                                              qz, absvec, li)
#    print orgmin
    #Find the best-fit source location in geographic coordinates
    lon2, lat2, z2 = qlon[minx], qlat[miny], qdep[minz]

    #Find the best-fit source location in geographic coordinates
#    print 'RMS: ',lon1,lat1,z1,'|| Origin: ',lon2,lat2,z2,' || ',ev['lon'],ev['lat'],ev['depth']

    return 0
