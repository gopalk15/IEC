import h5py
import matplotlib.pyplot as plt
import numpy as np
from config import Domain
from lib import get_params,plot_ion_density

cathode_radius = get_params('Grid','Cathode_Radius') # Cathode [m]
anode_radius = get_params('Grid','Anode_Radius')  # Anode [m]

fusor = Domain()
nx,ny = nodes = fusor.get_nodes()
ind_x,ind_y = cathode = fusor.build_grid(cathode_radius)
anode = fusor.build_grid(anode_radius)
dx = fusor.dx
dy = fusor.dy


''' Plot Chamber '''
cube = np.ones([nx,ny])*-2
idx = np.round(nx/2).astype(int)
idy = np.round(ny/2).astype(int)
cube[idx + 35 ,idy + 35 ] = 20
cube[idx + 15 ,idy + 25] = -20
PHI_BV2 = cube.T
xsh = PHI_BV2.shape[1]
ysh = PHI_BV2.shape[0]
xv = np.linspace(-fusor.chamber_radius, fusor.chamber_radius, xsh) 
yv = np.linspace(-fusor.chamber_height/2,fusor.chamber_height/2, ysh)
X, Y = np.meshgrid(xv, yv)
plt.figure(1)
plt.pcolormesh(X, Y, PHI_BV2, cmap='seismic',shading='auto')
plt.axis('equal')
plt.xlabel('Radial Direction [m]')
plt.ylabel('Height [m]')
plt.title('Computational Domain')

#open hdf5 file for reading 
with h5py.File('TestData2.h5','r') as hdf:
    base_items = list(hdf.items())
    print(f"Items in the base directory: {base_items}")
    G1 = hdf.get('DataSets/potential/set2')
    G1_items = list(G1.items())
    print(f"Items in Data Group: {G1_items}")
    position = np.array(G1.get('ParticlePosition'))
    den = np.array(G1.get('Density'))
    phi = np.array(G1.get('PHI'))
    efx = np.array(G1.get('electricFieldx'))
    efy = np.array(G1.get('electricFieldy'))
    G2 = hdf.get('DataSets/potential/set2/fusion')

    G2_items = list(G2.items())
    print(f"Items in Data Group: {G2_items}")

    fuse_pos = np.array(G2.get('Position'))
    fuse_rate = np.array(G2.get('RateData'))
    fuse_count = np.array(G2.get('FusionCount'))
    
    



x_data = position[:,0]

y_data = position[:,1]




plt.scatter(x_data,y_data)


plt.figure()
plot_ion_density(den,fusor.chamber_radius,fusor.chamber_height)

plt.figure()
plt.imshow(phi)

# plt.plot(xv*1000, den.T[round(den.shape[1]/2), :] / 1e15)
# plt.xlabel('Radial Direction [mm]')
# plt.ylabel('Ion Number Density [1^15 m^-3]')
plt.show()

print(f" count: {fuse_count}")
print(f"Rate Data: {fuse_rate}")
