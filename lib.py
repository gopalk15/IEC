import configparser
import numpy as np
import math
import matplotlib.pyplot as plt 


def get_params(section,variable,type=None):
    ''' Returns Paramenters Stored in .cfg file '''
    config = configparser.ConfigParser()
    config.read("simulation_parameters.cfg")
    config.sections()
    value = config[section][variable]

    if type=='int':
        return int(value)
    else:
        return float(value)
    

def get_discrete_values(y, dy):
    # Takes in an array y and a grid spacing dy and returns an array 
    # with the same physical dimensions as y but with a spacing of dy.
    actual = y / dy
    index = np.round(actual)
    return index*dy

def phi_matrix(cathode,anode,cathode_potential,anode_potenial,nx,ny):
    cathode_index = np.unique(cathode, axis=1)
    index_offset_x = np.round(nx/2).astype(int)
    index_offset_y = np.round(ny/2).astype(int)

    PHI = np.zeros([nx,ny])
    PHI[cat2[0] + idx, cat2[1] + idy] = cathode_potential
    PHI[xa_index + idx, ya_index + idy] = anode_potential


def plot_chamber(cathode,anode,chamber):
    plt.plot(cathode[0],cathode[1])
    plt.plot(anode[0],anode[1])
    plt.scatter(chamber[0][0],chamber[0][1])
    plt.scatter(chamber[1][0],chamber[1][1])
    plt.scatter(chamber[2][0],chamber[2][1])
    plt.scatter(chamber[3][0],chamber[3][1])
    plt.axes().set_aspect('equal', 'datalim')
    plt.show()


def XtoL(position,dx,dy):
    ''' Takes 2D dimentional postions of  '''
    lc = [position[0]/dx,position[1]/dy]
    return lc 

def nodal_position(lc,dx,dy):
    ''' Returns dimentional postions of nodes in meters'''
    pos = [lc[0]*dx,lc[1]*dy]
    return pos

def sampleIsotropicVel(vth):
    #pick a random angle
    theta = 2*np.pi*random()
    
    #pick a random direction for n[2]
    R = -1.0+2*random()
    a = np.sqrt(1-R*R)
    n = (np.cos(theta)*a, np.sin(theta)*a, R)
    
    #pick maxwellian velocities
    vm = np.zeros(2)
    vm[0:2] = np.sqrt(2)*vth*(2*(random()+random()+random()-1.5))
    
    vel = (n[0]*vm[0], n[1]*vm[1]) 
    return vel