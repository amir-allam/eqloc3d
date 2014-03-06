def run_3d_location_algorithm(stalist, evpha):
    # This script runs the whole earthquake location algorithm
    #  Just an outline for now

    #Initialize things? Set parameters? Read stations?
    from misc_tools import *
    from numpy import *
    import os.path

    nr= 46
    nlat=230
    nlon=250
    nx=nlon; ny=nlat; nz=nr
    li=Linear_index(nlon,nlat,nr)
    #Read station list
    stafnam='sta.phase_220to9_july13.dat'
    stalist=Stalist(stafnam)

    #Keep track of time
    start_time=time.time()
    print 'Starting earthquake location at ',time.ctime(start_time),'\n'

    #Read arrival times
    phafile='pha.phase_220to9_july13.dat'
    #phafile='tst_pha.dat'
    evpha=Phalist(phafile)

    #Read the header for the last traveltime file
    fnam=evpha[0]['arrivals'][0]['staname']+'.traveltime'
    hdr=Traveltime_header_file(fnam)
    #Build vectors of geographic coordinates
    qlon=arange(hdr.olon,hdr.dlon*hdr.nlon+hdr.olon,hdr.dlon)
    qlat=arange(hdr.olat,hdr.dlat*hdr.nlat+hdr.olat,hdr.dlat)
    qdep=arange(hdr.oz,hdr.dz*hdr.nz+hdr.oz,hdr.dz)

    #Grid search for best location
    for ev in evpha: #Loop over each event. We'll relocate everything on the list
        arrvec=array([]) #a vector of travel times
        absvec=array([]) # absolute arrival times
        arrsta=[] #a list of station names
        for arrival in ev['arrivals']: #Loop over arrivals and make vectors
    #        if not os.path.isfile(arrival['staname']+'traveltime'):
    #            continue
            arrvec=append(arrvec, float(arrival['ttime']) ) #Build vector of observed ttimes
            absvec=append(absvec, arrival['epoch'] ) #Absolute, not relative, times
            arrsta.append( arrival['staname'])
        if not arrsta: #Skip this event if no phases are present
            #print 'No arrivals for event ' +str(ev['id'])
            continue

        #Search coarsely
        dstep= 10;dx=nlon/dstep;dy=nlat/dstep;dz=nr/dstep;
        qx = range(1,nlon,dx); qy=range(1,nlat,dy); qz=range(1,nr,dz);
        minx,miny,minz=grid_search_traveltimes_rms(arrsta,qx,qy,qz,arrvec,li)

        #Finer search
        buff=15
        qx = range(minx-buff,minx+buff); qy=range(miny-buff,miny+buff); qz=range(minz-buff,minz+buff);
        qx=fix_boundary_search(qx,li.nx)
        qy=fix_boundary_search(qy,li.ny)
        qz=fix_boundary_search(qz,li.nz)
        minx,miny,minz=grid_search_traveltimes_rms(arrsta,qx,qy,qz,arrvec,li)
        #Find the best-fit source location in geographic coordinates
        lon1=qlon[minx]; lat1=qlat[miny]; z1=qdep[minz]

        #Search coarsely
        dstep= 10;dx=nlon/dstep;dy=nlat/dstep;dz=nr/dstep;
        qx = range(1,nlon,dx); qy=range(1,nlat,dy); qz=range(1,nr,dz);
        minx,miny,minz,orgmin=grid_search_traveltimes_origin(arrsta,qx,qy,qz,absvec,li)

        #Finer search
        buff=15
        qx = range(minx-buff,minx+buff); qy=range(miny-buff,miny+buff); qz=range(minz-buff,minz+buff);
        qx=fix_boundary_search(qx,li.nx)
        qy=fix_boundary_search(qy,li.ny)
        qz=fix_boundary_search(qz,li.nz)
        minx,miny,minz,orgmin=grid_search_traveltimes_origin(arrsta,qx,qy,qz,absvec,li)
        print orgmin
        #Find the best-fit source location in geographic coordinates
        lon2=qlon[minx]; lat2=qlat[miny]; z2=qdep[minz]

        #Find the best-fit source location in geographic coordinates
        print 'RMS: ',lon1,lat1,z1,'|| Origin: ',lon2,lat2,z2,' || ',ev['lon'],ev['lat'],ev['depth']
    elapsed_time=time.time()-start_time
    print 'Finished locating earthquake at ',time.ctime(time.time()),'\n'
    print 'Total location time: %8.4f seconds.\n' % (elapsed_time)
    return 0
