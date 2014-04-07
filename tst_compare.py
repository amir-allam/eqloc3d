#Compare the traveltime calculations of FMM and TomoDD
from numpy import asarray
from misc_tools import *
from mtools import *
import sys,os
import scipy.io
import numpy as np
import math
import utm

#Load events
fnam='pha.phase_220to9_july13.dat'
pha=Phalist(fnam,'scedc')

#Load stations
stafnam='sta.phase_220to9_july13.dat'
stalist=StationList(stafnam,False)

#Load dummy phase list
fnam='new_pha.dat'
outpha=Phalist(fnam,'scedc')

#fid=open('new_pha.dat','w')
iev=-1
for ev in pha:
    #Write out a new phase file for TomoDD including all possible stations
    #fid.write('# %i %i %i %i %i %8.4f %8.4f %8.4f %8.4f %8.4f 0.0 0.0 0.0 %i\n' % 
    #        (ev.year,ev.month,ev.day,ev.hour,ev.min,ev.sec,ev.lat,ev.lon,ev.depth,
    #            ev.mag,ev.id) )
    #for arr in ev.arrivals:
    #    fid.write('%s %8.4f %5.1f %s\n'% (arr.staname,arr.ttime,arr.qual,arr.phase) )
    #for sta in stalist: #Write out all stations with dummy values for calc_tt
    #    fid.write('%s 10.0 1.0 P\n'%sta.name)
    #
    #Write sources.in and run FMM for every event
    fid=open('sources.in','w')
    fid.write('1\n0\n %5.3f %8.4f %8.4f\n1\n1\n0 2\n1'%(ev.depth,ev.lat,ev.lon) )
    fid.close()
    outfile='./ttimes/%i.tt'%ev.id
    cmd='./fm3d>scratch'
    os.system(cmd)
    #Read the raw output, putting ttimes back into the pha struct
    fid=open('scratch','r')
    poop=fid.readlines()
    ic=-1 #Phase counter for output
    iev+=1#Event counter for output
    for line in poop:
        tmp=line.strip().split()
        if len(tmp)==0:
            continue
        elif tmp[0] == 'traced':
            ic+=1
            tt=float(tmp[10])
            outpha[iev].arrivals[ic].ttime=tt
    outpha[iev].show()
    fid.close()
#fid.close()
