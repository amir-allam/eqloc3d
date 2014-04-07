#Plot a netcdf topo file in python
from numpy import *
from scipy.io import netcdf_file as netcdf
from matplotlib import pyplot as plt
from misc_tools import *
import scipy.interpolate

fn='/Users/aallam/gmt/100m_poop.grd'
flon,flat=load_faults()

topo=netcdf(fn,'r')
subsamp=5
x=topo.variables['x'][::subsamp] #vector
y=topo.variables['y'][::subsamp] #vector
z=topo.variables['z'][::subsamp,::subsamp] #matrix
itopo=Interface()
itopo.read_topo_netcdf(fn)
gnn=scipy.interpolate.RectBivariateSpline(itopo.x,itopo.y,itopo.z.transpose())
qq=gnn.__call__(x,y)
print 'poo'
plt.figure(figsize=(6,6))
plt.plot(flon,flat,'k')
plt.axis([-118,-115,32,35])
plt.pcolor(x,y,z)
print 1

plt.figure(figsize=(6,6))
plt.axis([-118,-115,32,35])
plt.pcolor(x,y,qq.transpose() )
