import os
import sys
sys.path.append('%s/data/python' % os.environ['ANTELOPE'])
#from antelope.datascope import closing, dbopen
from antelope.stock import str2epoch, epoch2str
from numpy import *

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
    """This class is deprecated by StationList class"""
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

    def show(self): #Show the contents of the station list
        for sta in self:
            sta.show()

    class Sta():
    #a class for each individual station. There is almost certainly a more elegant way to do this.
        def __init__(self,name,lat,lon,elev):
            self.name=name; self.lat=lat; self.lon=lon; self.elev=elev
        def show(self): #show the contents of the class
            print '%s %5.4f %5.4f %5.2f'% (self.name,self.lon,self.lat,self.elev)

class Phalist(list):
    """This class is deprecated by PhaseList class."""
#A class containing the event information with arrival times for each event
    def __init__(self,fnam,ftype='ANTDB'):
        # fnam -    the file to be read
        # ftype -   the file type; right now either 'SCEDC' or 'ANTDB'
        if ftype is 'ANTDB' or ftype is 'antdb':
            self.read_antdb(fnam)
        elif ftype is 'SCEDC' or ftype is 'scedc':
            self.read_scedc(fnam)

    def read_antdb(self,fnam): #Read Antelope Database format (default)
        ####DOES NOT WORK YET######
        def __init__(self, db, evid): #Mal, this shouldn't take evid as an argument
            with closing(dbopen(db, 'r')) as db:
                view = db.schema_tables['origin']
                view = view.subset('evid == %s' % evid)
                view = view.subset('orid == prefor')
                view = view.join('netmag', outer=True)
                evid, time, lat, lon, depth, mag =  view.getv('evid', 'time', 'lat',
                    'lon', 'depth', 'magnitude')
                year = int(epoch2str(time, '%Y'))
                month = int(epoch2str(time, '%m'))
                day = int(epoch2str(time, '%d'))
                hour = int(epoch2str(time, '%H'))
                minute = int(epoch2str(time, '%M'))
                second = float(epoch2str(time, '%S.%s'))
                self.append({'id': evid, 'year': year, 'month': month, 'day': day,
                    'hour': hour, 'min': minute, 'sec': second, 'lat': lat,
                    'lon': lon, 'depth': depth, 'mag': mag, 'arrivals': []})
                view = view.join('assoc')
                view = view.join('arrival')
                for record in view.iter_record():
                    sta, arr_time, timeres, phase = record.getv('sta',
                        'arrival.time', 'timeres', 'phase')
                    ttime = arr_time - time
                    self[-1]['arrivals'].append({'staname': sta, 'ttime': ttime,
                        'qual': timeres, 'phase': phase}) 

    def read_scedc(self,fnam): #Read SCEC Datacenter format
        #Read the file
        print 'Reading ' +fnam
        fid = open(fnam,'r')
        a = fid.readlines()
        ic = -1
        for line in a: #Go line-by-line
            tmp=line.strip().split()
            if tmp[0] is '#': #In SCEDC format, an event line begins with #
                ic=ic+1
                self.append(self.Event( int(tmp[14]),int(tmp[1]),int(tmp[2]),int(tmp[3]),int(tmp[4]),int(tmp[5]),float(tmp[6]),float(tmp[8]),float(tmp[7]),float(tmp[9]),float(tmp[10]) ))
            else: #If it didn't begin with #, this line is a phase arrival for the last event
                self[ic].arrivals.append( self.Arrival(tmp[0],float(tmp[1]),float(tmp[2]),tmp[3],self[ic].epoch ))
        print 'Finished reading '+fnam

    class Event():
    #A class containing event metadata (if it exists)
        def __init__(self,evid,yr,mon,day,hr,mins,sec,lon,lat,dep,mag):
            self.id=evid; #event id
            self.year=yr; self.month=mon; self.day=day;  #event time
            self.hour=hr; self.min=mins;  self.sec=sec;  #more event time
            self.lon=lon; self.lat=lat;   self.depth=dep #event location
            self.mag=mag; #event magnitude
            #Calculate epoch time
            str_date=str(mon)+'/'+str(day)+'/'+str(yr)+' '+str(hr)+':'+str(mins)+':'+str(sec)
            print str_date
            self.epoch=str2epoch(str_date) #time in seconds after 1/1/1970
            self.arrivals=list() #This will be a list of phase arrivals
        def show(self): #show the contents of the class in stdout
            print '%u %u %u %u %u %u %5.4f %5.4f %5.4f %5.2f %5.2f %14.4f'% (self.id,self.year,self.month,self.day,self.hour,self.min,self.sec,self.lon,self.lat,self.depth,self.mag,self.epoch)
            for pha in self.arrivals:
                pha.show()

    class Arrival():
    #A class containing arrivals associated with an event; these will be entries in a
    #   list under Event.arrivals
        def __init__(self,staname,ttime,qual,phase,ev_epoch):
            self.staname=staname; self.ttime=ttime; self.qual=qual; self.phase=phase
            self.epoch=ev_epoch+ttime;
        def show(self): #show the contents of the class in stdout
            print '%s %5.4f %5.1f %s %14.4f'% (self.staname,self.ttime,self.qual,self.phase,self.epoch)
