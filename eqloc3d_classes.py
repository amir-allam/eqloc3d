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
        ret += 'lat: \t\t%s\n' % self.lat
        ret += 'lon: \t\t%s\n' % self.lon
        ret += 'elev: \t\t%s\n' % self.elev

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
            ret += 'lat:\t\t%s\n' % sta.lat
            ret += 'lon:\t\t%s\n' % sta.lon
            ret += 'elev:\t\t%s\n' % sta.elev
        return ret

    def _init_db(self, db):
        """
        Initialize station list using a CSS3.0 database as input.
        """
        with closing(dbopen(db, 'r')) as db:
            tbl_site = db.schema_tables['site']
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

class Event():
    """
    A container class for earthquake event metadata.
    """
    #def __init__(self, time, lat, lon, depth, mag, magtype=None, evid=None):
    def __init__(self,
                 prefor,
                 evid=None,
                 evname=None,
                 auth=None,
                 commid=None,
                 lddate=None,
                 origins=None):
        """
        Initialize Event object using one of two possible inputs.
        """
        import time as pytime
        self.evid = evid
        self.evname = evname
        self.prefor = prefor
        self.auth = auth
        self.commid = commid
        self.lddate = lddate
        self.preferred_origin = None
        if origins == None: self.origins = []
        else: self.origins = origins

    def __str__(self):
        ret = 'Event Object\n------------\n'
        ret += 'evid:\t\t%s\n' % self.evid
        ret += 'evname:\t\t%s\n' % self.evname
        ret += 'prefor:\t\t%s\n' % self.prefor
        ret += 'auth:\t\t%s\n' % self.auth
        ret += 'commid:\t\t%s\n' % self.commid
        ret += 'lddate:\t\t%s\n' % self.lddate
        ret += 'origins:\n'
        if len(self.origins) == 0:
            ret += '\t\tNone\n'
        else:
            for i in range(len(self.origins)):
                for line in  ('%s' % self.origins[i]).split('\n'):
                    ret += '\t\t%s\n' % line
        return ret

    def set_preferred_origin(self, prefor):
        """
        Set self.preferred_origin to equal origin with orid == prefor.
        """
        for i in range(len(self.origins)):
            if self.origins[i].orid == prefor:
                self.preferred_origin = self.origins[i]
                return 0

    def add_origin(self,
                   lat,
                   lon,
                   depth,
                   time,
                   auth,
                   arrivals=[],
                   orid=None,
                   evid=None,
                   jdate=None,
                   nass=None,
                   ndef=None,
                   ndp=None,
                   grn=None,
                   srn=None,
                   etype=None,
                   review=None,
                   depdp=None,
                   dtype=None,
                   mb=None,
                   mbid=None,
                   ms=None,
                   msid=None,
                   ml=None,
                   mlid=None,
                   algorithm=None,
                   commid=None,
                   lddate=None):
        """
        Add an Origin object to the list of origins associated with this event.
        """
        self.origins += [Origin(lat,
                                lon,
                                depth,
                                time,
                                auth,
                                orid=orid,
                                evid=evid,
                                arrivals=arrivals,
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
                                lddate=lddate)]
    def write(self, out, fmt):
        """
        Write out newly authored data pertaining to event.

        Arguments:
        out - A datascope db pointer to an open CSS3.0 database for output.
        Alternatively the path to an output SCEDC format flat file.

        fmt - The format of the 'out' argument; 'CSS3.0' or 'SCEDC'.

        Return Values:
        0 - Success
        -1 - Failure
        """
        if fmt == 'CSS3.0':
            return self._write_CSS()
        elif fmt == 'SCEDC':
            #Do it the SCEDC way
            pass
        else:
            raise Exception('Output format %s not recognized' % fmt)
    def _write_css(self, db):
        """
        Write newly authored data to a CSS3.0 database.

        Arguments:
        db - Output database.
        """
        #There are three cases that must be handled when writing out to
        #a CSS3.0 schema database.
        #Case 1: Reading from and writing to SAME database.
        #Case 2: Reading from and writing to SEPARATE databases.
        #Case 3: Reading from SCEDC format file and writing to CSS3.0
        #database.
        for origin in self.origins:
            if origin.auth == 'eqloc3d':
                tbl_origin = out.schema_tables['origin']
                tbl_origin.record = tbl_origin.addnull()
                tbl_origin.putv('lat', origin.lat,
                                'lon', origin.lon,
                                'depth', origin.depth,
                                'time', origin.time,
                                'orid', origin.orid,
                                'evid', origin.evid,
                                'auth', origin.auth,
                                'jdate', origin.jdate,
                                'nass', origin.nass,
                                'ndef', origin.ndef,
                                'ndp', origin.ndp,
                                'grn', origin.grn,
                                'srn', origin.srn,
                                'etype', origin.etype,
                                'review', origin.review,
                                'depdp', origin.depdp,
                                'dtype', origin.dtype,
                                'mb', origin.mb,
                                'mbid', origin.mbid,
                                'ms', origin.ms,
                                'msid', origin.msid,
                                'ml', origin.ml,
                                'mlid', origin.mlid,
                                'algorithm', origin.algorithm,
                                'commid', origin.commid,
                                'lddate', origin.lddate)
                tbl_assoc = out.schema_tables['assoc']
                for arrival in origin.arrivals:
                    if arrival.arid == None:
                        pass
                        #arrival is from SCEDC source add row to arrival table
        return 0

class Origin():
    """
    A container class for origin data.
    """
    def __init__(self,
                 lat,
                 lon,
                 depth,
                 time,
                 auth,
                 arrivals=[],
                 orid=None,
                 evid=None,
                 jdate=None,
                 nass=None,
                 ndef=None,
                 ndp=None,
                 grn=None,
                 srn=None,
                 etype=None,
                 review=None,
                 depdp=None,
                 dtype=None,
                 mb=None,
                 mbid=None,
                 ms=None,
                 msid=None,
                 ml=None,
                 mlid=None,
                 algorithm=None,
                 commid=None,
                 lddate=None):
        self.lat = lat
        self.lon = lon
        self.depth = depth
        self.time = time
        self.orid = orid
        self.evid = evid
        self.auth = auth
        self.arrivals = arrivals
        self.jdate = jdate
        self.nass = nass
        self.ndef = ndef
        self.ndp = ndp
        self.grn = grn
        self.srn = srn
        self.etype = etype
        self.review = review
        self.depdp = depdp
        self.dtype = dtype
        self.mb = mb
        self.mbid = mbid
        self.ms = ms
        self.msid = msid
        self.ml = ml
        self.mlid = mlid
        self.algorithm = algorithm
        self.commid = commid
        self.lddate = lddate

    def __str__(self):
        """
        Return string representation of Origin object.
        """
        ret = 'Origin Object\n-------------\n'
        ret += 'lat:\t\t%s\n' % self.lat
        ret += 'lon:\t\t%s\n' % self.lon
        ret += 'depth:\t\t%s\n' % self.depth
        ret += 'time:\t\t%s\n' % self.time
        ret += 'orid:\t\t%s\n' % self.orid
        ret += 'evid:\t\t%s\n' % self.evid
        ret += 'auth:\t\t%s\n' % self.auth
        ret += 'jdate:\t\t%s\n' % self.jdate
        ret += 'nass:\t\t%s\n' % self.nass
        ret += 'ndef:\t\t%s\n' % self.ndef
        ret += 'ndp:\t\t%s\n' % self.ndp
        ret += 'grn:\t\t%s\n' % self.grn
        ret += 'srn:\t\t%s\n' % self.srn
        ret += 'etype:\t\t%s\n' % self.etype
        ret += 'review:\t\t%s\n' % self.review
        ret += 'depdp:\t\t%s\n' % self.depdp
        ret += 'dtype:\t\t%s\n' % self.dtype
        ret += 'mb:\t\t%s\n' % self.mb
        ret += 'mbid:\t\t%s\n' % self.mbid
        ret += 'ms:\t\t%s\n' % self.ms
        ret += 'msid:\t\t%s\n' % self.msid
        ret += 'ml:\t\t%s\n' % self.ml
        ret += 'mlid:\t\t%s\n' % self.mlid
        ret += 'algorithm:\t\t%s\n' % self.algorithm
        ret += 'commid:\t\t%s\n' % self.commid
        ret += 'lddate:\t\t%s\n' % self.lddate
        ret += 'arrivals:\n'
        for i in range(len(self.arrivals)):
            ret += '\t\t%s' % self.arrivals[i]
        return ret

class Phase():
    """
    A container class for phase data.
    """
    def __init__(self, sta, time, phase, qual=None, arid=None):
        self.sta = sta
        self.time = time
        self.phase = phase
        self.qual = qual
        self.arid = arid

    def __str__(self):
        ret = 'Arrival Object\n--------------\n'
        ret += 'sta:\t\t%s\n' % self.sta
        ret += 'time:\t\t%s\n' % self.time
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
