#Test the translation from 1D to 3D indices
from misc_tools import *

#data = read_arrtimes(fnam='BAC.traveltime')
#data=read_arrtimes()
ny,nx,nz=data.shape

lli=Linear_index(nx,ny,nz)

vdata=[]
for ix in range(nx):
    for iy in range(ny):
        for iz in range(nz):
            vdata.append( data[iy,ix,iz] )

st=stalist[38]
st.show()

#Find the nearest node to the station of interest
tmp=misc_tools.find_nearest(qlon,st.lon)
ilon=nonzero(qlon==tmp)[0][0]
tmp=misc_tools.find_nearest(qlat,st.lat)
ilat=nonzero(qlat==tmp)[0][0]
mm=misc_tools.find_nearest(qdep,earth_rad+st.elev/1000)
idep=nonzero(qdep==mm)[0][0]
print st.name,qlon[ilon],qlat[ilat],qdep[idep]

#Read the arrival time there
ind =li.get_1D(ilon,ilat,idep)
print vdata[ind]
#print read_binary_float(fid,ind)
