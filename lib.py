import configparser
import numpy as np
import math
from random import random
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

def plot_chamber(cathode,anode,chamber):
    plt.plot(cathode[0],cathode[1])
    plt.plot(anode[0],anode[1])
    plt.scatter(chamber[0][0],chamber[0][1])
    plt.scatter(chamber[1][0],chamber[1][1])
    plt.scatter(chamber[2][0],chamber[2][1])
    plt.scatter(chamber[3][0],chamber[3][1])
    plt.axes().set_aspect('equal', 'datalim')
    plt.show()


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

def fusion_cross_section(vx, vy,massIon):
    # Takes in velocity componants in m/s and returns a cross section in barns
    E = .5*massIon*(vx**2 + vy**2)
    E = E*6.242e15  # convert J to KeV 
    A1 = 46.097
    A2 = 372
    A3 = 4.36e-4
    A4 = 1.22
    A5 = 0
    AA1 = 47.88
    AA2 = 482
    AA3 = 3.08e-4
    AA4 = 1.177
    AA5 = 0
    term1 = A5 + A2/((A4 - A3*E)**2 + 1)
    term2 = E*(np.exp(A1/np.sqrt(E)) - 1)
    term3 = AA5 + AA2/((AA4 - AA3*E)**2 + 1)
    term4 = E*(np.exp(AA1/np.sqrt(E)) - 1)
    cs = ((term1/term2) + (term3/term4))*1e-28
    return cs

def plot_ion_density(density,chamber_radius,chamber_height):
    PHI_B2 = density.T
    xsh = PHI_B2.shape[1]
    ysh = PHI_B2.shape[0]
    xv = np.linspace(-chamber_radius, chamber_radius, xsh) 
    yv = np.linspace(0, chamber_height, ysh)
    X, Y = np.meshgrid(xv, yv)
    plt.contourf(X, Y, PHI_B2, levels=60, cmap='inferno')
    #plt.pcolor(X, Y, PHI_B2, cmap='inferno')
    plt.colorbar()
    plt.axis('equal')
    print(PHI_B2.shape)
    plt.xlabel('Radial Direction [m]')
    plt.ylabel('Height [m]')
    plt.title('Ion Number Density')

def gather(Matrix,particles_position,index,chamber_height,chamber_radius,dh):

    fi = (particles_position[index, 0] + chamber_radius-dh)/dh    #real i index of particle's cell
    i = np.floor(fi).astype(int)          #integral part
    hx = fi - i                              #the remainder
    fj = (particles_position[index,1] + (chamber_height/2)-dh)/dh #real i index of particle's cell
    j = np.floor(fj).astype(int)              #integral part
    hy = fj - j 

    #  Tk = Te*11604.45 #Kelvin
    # chamber_pressure = pressure*133.322 #Pascal

    rho = Matrix[i,j]*(1-hx)*(1-hy)
    rho = rho + Matrix[i+1,j]*hx*(1-hy)
    rho = rho + Matrix[i,j+1]*(1-hx)*hy
    rho = rho + Matrix[i+1,j+1]*hx*hy
    return rho
