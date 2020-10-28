import configparser
import numpy as np
import matplotlib.pyplot as plt 


def get_params(section,variable):
    ''' Returns Paramenters Stored in .cfg file '''
    config = configparser.ConfigParser()
    config.read("simulation_parameters.cfg")
    config.sections()
    return float(config[section][variable])
    

def get_discrete_values(y, dy):
    # Takes in an array y and a grid spacing dy and returns an array 
    # with the same physical dimensions as y but with a spacing of dy.
    actual = y / dy
    index = np.round(actual)
    return index*dy

def plot_domain(cathode,anode,chamber):
    plt.plot(cathode[0],cathode[1])
    plt.plot(anode[0],anode[1])
    plt.scatter(chamber[0][0],chamber[0][1])
    plt.scatter(chamber[1][0],chamber[1][1])
    plt.scatter(chamber[2][0],chamber[2][1])
    plt.scatter(chamber[3][0],chamber[3][1])
    plt.show()