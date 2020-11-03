from config import Domain,Particles,ESPIC,Fusion
import numpy as np
from lib import get_params,plot_chamber,XtoL
import matplotlib.pyplot as plt

#Variable Inputs (for later)
potential = -100000   # [V] cathode potential 
cathode_radius = get_params('Grid','Cathode_Radius') # Cathode [m]
anode_radius = get_params('Grid','Anode_Radius')  # Anode [m]
pressure = 7      # Pa (Not used yet)
Te = np.abs(potential)    #electron temperature in eV


''' Initialise '''

fusor = Domain()
nodes = fusor.get_nodes()
cathode = fusor.build_grid(cathode_radius)
anode = fusor.build_grid(anode_radius)

loop = ESPIC(cathode,anode,nodes,600,potential)
particles = Particles(nodes)
fusion = Fusion()

PHI_G = np.zeros(nodes)
loop.build_potential_matrix(PHI_G)

























''' Visulize compuational domian '''
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








# ind_x = (cathode[0]/fusor.dx).astype(int) 
# ind_y = (cathode[1]/fusor.dy).astype(int) 

# an_x = (anode[0]/fusor.dx).astype(int)
# an_y = (anode[1]/fusor.dy).astype(int)

# print(ind_x)
# print(len(ind_x))
# # # Deleate repeate XY Pairs
# Cathode_Id1 = np.zeros([2, len(ind_x)])
# Cathode_Id1[0, :] = ind_x
# Cathode_Id1[1, :] = ind_y
# Cathode_Id2 = np.unique(Cathode_Id1, axis=1)
# ind_x = Cathode_Id2[0, :].astype(int)
# ind_y = Cathode_Id2[1, :].astype(int)

# print(len(ind_y))
# Anode_Id1 = np.zeros([2, len(an_x)])
# Anode_Id1[0, :] = an_x
# Anode_Id1[1, :] = an_y
# Anode_Id2 = np.unique(Cathode_Id1, axis=1)
# an_x = Anode_Id2[0, :].astype(int)
# an_y = Anode_Id2[1, :].astype(int)


# cube = np.ones([nx,ny])*-2
# idx = np.round(nx/2).astype(int)
# idy = np.round(ny/2).astype(int)
# cube[ind_x + 35 ,ind_y + 35 ] = 20
# cube[an_x + 15 ,an_y + 25] = -20
# PHI_BV2 = cube.T
# xsh = PHI_BV2.shape[1]
# ysh = PHI_BV2.shape[0]
# xv = np.linspace(-fusor.chamber_radius, fusor.chamber_radius, xsh) 
# yv = np.linspace(0, fusor.chamber_height, ysh)
# X, Y = np.meshgrid(xv, yv)
# plt.pcolormesh(X, Y, PHI_BV2, cmap='seismic',shading='auto')
# plt.axis('equal')
# plt.xlabel('Radial Direction [m]')
# plt.ylabel('Height [m]')
# plt.title('Computational Domain')
# plt.show()









# ''' Visulize compuational domian '''
# PHI_BV = np.ones([nx, ny])*-2
# idx = np.round(nx/2).astype(int)
# idy = np.round(ny/2).astype(int)
# PHI_BV[i_cat_x + idx, i_cat_y + idy] = 20
# PHI_BV[i_an_x + idx, i_an_y + idy] = -20
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


