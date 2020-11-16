import sys
sys.path.insert(1,'C:\\Users\\ellen\\git\\IEC')
from config import Domain,Particles 
import numpy as np
from lib import get_params,plot_chamber
import matplotlib.pyplot as plt


#Variable Inputs (for later)
potential = -100000   # [V] cathode potential 
cathode_radius = get_params('Grid','Cathode_Radius') # Cathode [m]
anode_radius = get_params('Grid','Anode_Radius')  # Anode [m]
radius = get_params('Source','sourceRadius') - 0.2
pressure = 7      # Pa (Not used yet)
Te = np.abs(potential)    #electron temperature in eV

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
plt.pcolormesh(X, Y, PHI_BV2, cmap='seismic',shading='auto')
plt.axis('equal')
plt.xlabel('Radial Direction [m]')
plt.ylabel('Height [m]')
plt.title('Computational Domain')








particles = Particles(nodes)

spray_radius,angle = particles.get_spray_values()

particles.generate(radius)

x_part = particles.pos[:,0]
y_part = particles.pos[:,1]

plt.scatter(x_part,y_part)

plt.axis('equal')
plt.show()

# # # Deleate repeate XY Pairs
# cat2 = np.zeros([2, len(ind_x)])
# cat2[0, :] = ind_x
# cat2[1, :] = ind_y
# cat2 = np.unique(cat2, axis=1)
# ind_x = cat2[0, :].astype(int)
# ind_y = Cathode_Id2[1, :].astype(int)

# ''' Visulize compuational domian '''
# PHI_BV = np.ones([nx, ny])*-2
# idx = np.round(nx/2).astype(int)
# idy = np.round(ny/2).astype(int)
# PHI_BV[cat2[0] + idx, cat2[1] + idy] = 20
# PHI_BV[xa_index + idx, ya_index + idy] = -20
# PHI_BV2 = PHI_BV.T
# xsh = PHI_BV2.shape[1]
# ysh = PHI_BV2.shape[0]
# xv = np.linspace(-fusor.chamber_radius, fusor.chamber_radius, xsh) 
# yv = np.linspace(0, fusor.chamber_height, ysh)
# X, Y = np.meshgrid(xv, yv)
# plt.pcolor(X, Y, PHI_BV2, cmap='seismic')
# plt.axis('equal')
# plt.xlabel('Radial Direction [m]')
# plt.ylabel('Height [m]')
# plt.title('Computational Domain')
# plt.show()