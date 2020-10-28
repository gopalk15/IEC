from config import Domain
import numpy as np
from lib import get_params,plot_domain
import matplotlib.pyplot as plt

#Variable Inputs (for later)
potential = -100000   # [V] cathode potential 
cathode_radius = get_params('Grid','Cathode_Radius') # Cathode [m]
anode_radius = get_params('Grid','Anode_Radius')  # Anode [m]
pressure = 7      # Pa (Not used yet)

Te = np.abs(potential)    #electron temperature in eV

fusor = Domain()

cathode = fusor.build_grid(cathode_radius)
anode = fusor.build_grid(anode_radius)
chamber = fusor.build_chamber()

plot_domain(cathode,anode,chamber)