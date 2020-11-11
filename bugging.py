import numpy as np
import matplotlib.pyplot as plt

#setup constants
EPS0 = 8.854e-12       #permittivity of free space
QE = 1.602e-19         #elementary charge
kb = 1                 #boltzmann constant
AMU = 1.661e-27        #atomic mass unit
m_ion = 2*AMU          #ion mass (deuterium)

# Physical Grid Dimensions
R1 = 0.05            # Cathode [m]
GEO_C = .95          # Geometric transparency cathode (manual for now)
R2 = .085            # Anode [m]
GEO_A = .95          # Geometric transparency anode   (manual for now)
#Physical Chamber Dimensions
R_Chamber = .1       # [m]
H_Chamber = 0.5     # [m]
r_wire = .80*1e-3 / 2 # Radius of 20 guage wire in m 
#Sourse Dimension
R_Sourse = .08 

#input settings
n0 = 4.6e13                 #electron background density in #/m^3
phi_Cathode = -100000        #cathode potential
phi0 = 0                    #reference potential 
Te = np.abs(phi_Cathode)    #electron temperature in eV
Ti = 0.1                    #ion velocity in eV (not used yet)
vth = np.sqrt(2*QE*Ti/m_ion)   #thermal velocity with Ti in eV
Operating_Pressure = 7      # Pa (Not used yet)
#calculate plasma parameters
#lD = np.sqrt(EPS0*Te/(n0*QE))      #Debye length (not used yet)
vth = np.sqrt(2*QE*Ti/m_ion)        #Thermal velocity with Ti in eV (not used yet)
lD = 0.004   # MANUAL SETTING FOR TESTING
dx = lD               #x cell size
dy = lD               #y cell size
print(lD)
#set simulation domain in the x dimension  
nx = np.floor((2*R_Chamber)/dx).astype(int) #number of nodes in x direction, with a physical 
                                            #domain +/- R_Chamber
#set simulation domain in the y dimension  
ny = np.floor((H_Chamber)/dy).astype(int)   #number of nodes in y direction, with a physical 
                                            #domain +/- H_Chamber / 2
ts = 800                                    #number of time steps

#calculate maximum expected velocity and timestep
E_av = (np.abs(phi_Cathode) - 0) / (R2 - R1)
a_av = E_av*QE / m_ion
v_max = np.sqrt(2*a_av*R2)
dt = 1e-10     #time step size, at vmax move 0.10dx

# Domain length in the x direction
Lx = (nx-1)*dx
# Domain length in the y direction 
Ly = (ny-1)*dy
print(nx)
print(dt)
print(vth)

# SET UP THE POINTS REPRESENTING THE GRIDS
def get_discrete_Value(y, dy):
    # Takes in an array y and a grid spacing dy and returns an array 
    # with the same physical dimensions as y but with a spacing of dy.
    actual = y / dy
    index = np.round(actual)
    return index*dy

# CATHODE SETUP
Thetav = np.arange(0, 2*np.pi, np.pi/1000)
x_cathode = R1*np.cos(Thetav)
y_cathode = R1*np.sin(Thetav)
X_cathode = get_discrete_Value(x_cathode, dx)
Y_cathode = get_discrete_Value(y_cathode, dy)

INDEX_X_Cathode = (X_cathode/dx).astype(int)
INDEX_Y_Cathode = (Y_cathode/dx).astype(int)

# ANODE SETUP
x_anode = R2*np.cos(Thetav)
y_anode = R2*np.sin(Thetav)
X_anode = get_discrete_Value(x_anode, dx)
Y_anode = get_discrete_Value(y_anode, dy)

INDEX_X_Anode = (X_anode/dx).astype(int)
INDEX_Y_Anode = (Y_anode/dy).astype(int)

# LEFT WALL SETUP
y_left_wall = np.arange(-H_Chamber/2, H_Chamber/2, dy)
x_left_wall = np.ones(y_left_wall.shape[0]) * (-R_Chamber)
INDEX_X_left_wall = (x_left_wall/dx).astype(int)
INDEX_Y_left_wall = (y_left_wall/dy).astype(int)
# RIGHT WALL SETUP
y_right_wall = np.arange(-H_Chamber/2, H_Chamber/2, dy)
x_right_wall = np.ones(y_right_wall.shape[0]) * (R_Chamber - dx)
INDEX_X_right_wall = (x_right_wall/dx).astype(int)
INDEX_Y_right_wall = (y_right_wall/dy).astype(int)
# BOTTOM WALL SETUP
x_bottom_wall = np.arange(-R_Chamber, R_Chamber, dx)
y_bottom_wall = np.ones(x_bottom_wall.shape[0]) * (-H_Chamber/2)
INDEX_X_bottom_wall = (x_bottom_wall/dx).astype(int)
INDEX_Y_bottom_wall = (y_bottom_wall/dy).astype(int)
# TOP WALL SETUP
x_top_wall = np.arange(-R_Chamber, R_Chamber, dx)
y_top_wall = np.ones(x_bottom_wall.shape[0]) * (H_Chamber/2 - dy)
INDEX_X_top_wall = (x_top_wall/dx).astype(int)
INDEX_Y_top_wall = (y_top_wall/dy).astype(int)
print(len(INDEX_X_Cathode))
print(INDEX_X_Cathode)

# Deleate repeate XY Pairs
Cathode_Id1 = np.zeros([2, INDEX_X_Cathode.shape[0]])
Cathode_Id1[0, :] = INDEX_X_Cathode
Cathode_Id1[1, :] = INDEX_Y_Cathode
print(Cathode_Id1)
Cathode_Id2 = np.unique(Cathode_Id1, axis=1)
INDEX_X_Cathode = Cathode_Id2[0, :].astype(int)
INDEX_Y_Cathode = Cathode_Id2[1, :].astype(int)
print(Cathode_Id2)

# Calculate specific weight and prepair to insert particles
np_insert = 400                      #insert 2 particles per anode cell.
flux = 4.6e20                          #Flux of entering ions [ions per second]
npt = flux*dt
spwt = npt/np_insert                    #specific weight, real particles per macroparticle
mp_q = 1                                #macroparticle charge
max_part = 200000                       #buffer size
#allocate particle array
part_x = np.zeros([max_part,2]) #particle positions
part_v = np.zeros([max_part,2]) #particle velocities
sourse_storage = np.zeros([max_part,2]) #particle sourses 

# plt.scatter(x_right_wall, y_right_wall)
# plt.scatter(x_left_wall, y_left_wall)
# plt.scatter(x_top_wall, y_top_wall)
# plt.scatter(x_bottom_wall, y_bottom_wall)
# plt.scatter(X_anode, Y_anode,  marker='*')
# plt.plot(x_anode, y_anode, 'b')
# plt.scatter(X_cathode, Y_cathode,  marker='*')
# plt.plot(x_cathode, y_cathode, 'r')
# plt.axis('equal')

# PHI_BV = np.ones([nx, ny])*-2
# idx = np.round(nx/2).astype(int)
# idy = np.round(ny/2).astype(int)
# print(len(INDEX_X_Anode))
# print(idy)
# PHI_BV[INDEX_X_Cathode + idx, INDEX_Y_Cathode + idy] = 20
# PHI_BV[INDEX_X_Anode + idx, INDEX_Y_Anode + idy] = -20
# PHI_BV2 = PHI_BV.T
# xsh = PHI_BV2.shape[1]
# ysh = PHI_BV2.shape[0]
# xv = np.linspace(-R_Chamber, R_Chamber, xsh) 
# yv = np.linspace(0, H_Chamber, ysh)
# X, Y = np.meshgrid(xv, yv)
# plt.pcolor(X, Y, PHI_BV2, cmap='seismic')
# plt.axis('equal')
# plt.xlabel('Radial Direction [m]')
# plt.ylabel('Height [m]')
# plt.title('Computational Domain')

# plt.show()