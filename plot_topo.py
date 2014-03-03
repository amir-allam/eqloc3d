#Plot a netcdf topo file in python
from numpy import *
from scipy.io import netcdf_file as netcdf
from matplotlib import pyplot as plt
from misc_tools import *

fn='/Users/aallam/gmt/100m_poop.grd'
flon,flat=load_faults()

topo=netcdf(fn,'r')
subsamp=5
x=topo.variables['x'][::subsamp] #vector
y=topo.variables['y'][::subsamp] #vector
z=topo.variables['z'][::subsamp,::subsamp] #matrix
print 'poo'
plt.figure(figsize=(6,6))
print '1'
plt.plot(flon,flat,'k')
print '2'
plt.axis([-118,-115,32,35])
print '3'
plt.pcolor(x,y,z)
print '4'
