import sys
import os
sys.path.append('%s/data/python' % os.environ['ANTELOPE'])
import logging
import time
import subprocess
from misc_tools import *
from mtools import Stalist

tt_calculator = '../fortran/fm3d'

def _parse_args():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('db', type=str, help='Input database.')
    parser.add_argument('-b', '--binary', action='store_true',
        help='Write binary output file.')
    return parser.parse_args()

def _configure_logger():
    logger = logging.getLogger(sys.argv[0])
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter(fmt='%(asctime)s::%(levelname)s::%(message)s',
        datefmt='%Y%j %H:%M:%S')
    fh = logging.FileHandler('%s.log' % sys.argv[0])
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    sh = logging.StreamHandler()
    sh.setLevel(logging.INFO)
    sh.setFormatter(formatter)
    logger.addHandler(sh)

def run_fmm():
    #Run the traveltime solver with all of the input files currently in the directory
    #  Note: Input files should already be set using write_sources_in, etc...
        subprocess.call('./fm3d',shell=True)

def gen_sta_tt_maps(stalist,if_write_binary=True):
#Generate a travel time map for each station in the station list
    start_time=time.time()
    print 'Starting travel time calculation at ',start_time,'\n'
    for sta in stalist:
        print 'Generating travel times for station', sta.name, '\n'
        #Elevation can be set to a large negative number to glue the source to the surface
        _write_sources_file(sta.elev * -1, sta.lat, sta.lon) #Elevation is in km and negative
        #_write_sources_file(0.0,sta.lat,sta.lon) #!!!! SET SOURCE TO 0 DEPTH
        run_fmm()
        #Create output file name
        outfnam=sta.name+'.traveltime'
        subprocess.call('mv arrtimes.dat '+outfnam,shell=True)
        if if_write_binary:
            _tt_ascii_2_binary(outfnam)
    elapsed_time=time.time()-start_time
    print 'Finished travel time calculations at ',time.time(),'\n'
    print 'Total calculation time: %8.4f seconds.\n' % (elapsed_time)

def _generate_tt_maps(db, write_binary=True):
    logger = logging.getLogger(sys.argv[0])
    logger.debug('Begin travel-time map generateion.')
    with closing(dbopen(db, 'r')) as db:
        tbl_site = db.schema_tables['site']
        for record in tbl_site.iter_record():
            sta, lat, lon, elev = record.getv('sta', 'lat', 'lon', 'elev')
            logger.debug('Begin travel-time map generation for station %s'
                % sta)
            _write_sources_file(0.10, lat, lon)
            os.system(tt_calculator)
            logger.debug('End travel-time map generation for station %s'
                % sta)
#            outfile = '%s.traveltime' % sta
#            os.system('mv arrtimes.dat %s' % outfile)
#            if write_binary:
#                logger.debug('Begin writing binary file for station % s' % sta)
#                _tt_ascii_2_binary(outfile)
#                logger.debug('End writing binary file for station % s' % sta)

def _write_sources_file(depth, lat, lon):
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
    outfile = open('sources.in','w')
    outfile.write(" %i\n %i\n %.4f %.4f %.4f\n %i\n %i\n %i %i\n %i\n"
        % (numsrc, teleflag, depth, lat, lon, numpaths, numsections, 0, 2,
        veltype))
    outfile.close()
    #Path sections are described in more detail in the FMM README
    # For the first arrival from a source propagating at the surface [ 0 2 ]
    #  is the path section

def _tt_ascii_2_binary(fnam):
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

def gen_tt_map_malcolm(db):
    from antelope.datascope import closing, dbopen
    station_list = StationList(db, is_db=True)
    gen_sta_tt_maps(station_list)
