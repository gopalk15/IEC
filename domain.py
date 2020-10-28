import numpy as np
import matplotlib.pyplot as plt
from parameters import *

Ti = 0.1                    #ion velocity in eV (not used yet)
#calculate plasma parameters
#lD = np.sqrt(EPS0*Te/(n0*QE))      #Debye length (not used yet)
vth = np.sqrt(2*QE*Ti/m_ion)        #Thermal velocity with Ti in eV (not used yet)
lD = 0.004   # MANUAL SETTING FOR TESTING
dx = lD               #x cell size
dy = lD               #y cell size
print(lD)



#set simulation domain in the x dimension
nx = int(np.floor((2*R_Chamber)/dx)) #number of nodes in x direction, with a physical
                                            #domain +/- R_Chamber
#set simulation domain in the y dimension
ny = int(np.floor((H_Chamber)/dy))  #number of nodes in y direction, with a physical
                                            #domain +/- H_Chamber / 2
ts = 800                                    #number of time steps





# Domain length in the x direction
Lx = (nx-1)*dx
# Domain length in the y direction
Ly = (ny-1)*dy
print(nx)
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

# Deleate repeate XY Pairs
Cathode_Id1 = np.zeros([2, INDEX_X_Cathode.shape[0]])
Cathode_Id1[0, :] = INDEX_X_Cathode
Cathode_Id1[1, :] = INDEX_Y_Cathode
Cathode_Id2 = np.unique(Cathode_Id1, axis=1)
INDEX_X_Cathode = Cathode_Id2[0, :].astype(int)
INDEX_Y_Cathode = Cathode_Id2[1, :].astype(int)



plt.scatter(x_right_wall, y_right_wall)
plt.scatter(x_left_wall, y_left_wall)
plt.scatter(x_top_wall, y_top_wall)
plt.scatter(x_bottom_wall, y_bottom_wall)
plt.scatter(X_anode, Y_anode,  marker='*')
plt.plot(x_anode, y_anode, 'b')
plt.scatter(X_cathode, Y_cathode,  marker='*')
plt.plot(x_cathode, y_cathode, 'r')
plt.axis('equal')
plt.show()
