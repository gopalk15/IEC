import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
from extract_data import Data
import dataLib as lib



sim_type = ['ratio','potential','radius']

for sim in sim_type:
    print(f"Processing Data ({sim})...")
    df = pd.DataFrame(columns=['Potential (keV)','Cathode Radius (m)','Anode Radius (m)','Fusion Rate (m^2/sec)'])
    if sim == 'ratio':
        ''' Ratio test: Anode Grid kept constant while cathode is varied '''
        potential = np.ones(10)*-100 #[keV]
        anode_radius = np.ones(10)*0.085 # [m]
        ratio = np.array([0.9,10/17,0.8,0.5,0.1,1/3,3/4,3/7,0.2,0.95])
        cathode_radius = np.multiply(anode_radius,ratio)
    elif sim == 'radius':
        ''' Radius Test: Cathode grid is kept constant while the anode is varied '''
        potential = np.ones(10)*-100 #[keV]
        cathode_radius =np.ones(10)*0.05 #[m]
        ratio = np.array([1.7,1.2,2,2.5,3,3.5,2.2,2.8,1.5,3.2]) 
        anode_radius = np.multiply(cathode_radius,ratio)
    elif sim == 'potential':
        ''' Potential Test: Radius of grids held constant while potential is  varied'''
        potential = np.array([-80,-90,-100,-110,-120,-130,-140,-150,-160,-170])
        anode_radius = np.ones(10)*0.085 # [m]
        cathode_radius = np.ones(10)*0.05 #[m]

    for k in range(10):

        # Create data object 
        data = Data(sim,(k+1))

        # Rettrive Data
        DEN = data.get_density_data()
        velocity = data.get_particle_velocity()
        position = data.get_particle_position()

        #Process Data
        speed = lib.get_speed(velocity)
        sigma = lib.get_cross_section(speed)
        rho = lib.gather_density(position,DEN)

        #Create Distribution 
        bin_array = np.linspace(0,10,100)*1e07
        normalized_particles,bin_edges = np.histogram(speed,bins=bin_array,density=True)
        x_scatter = bin_edges[:-1]+ 0.5*(bin_edges[1:] - bin_edges[:-1])

        #Find Fusion rate
        sig_v = lib.get_vel_average(x_scatter,normalized_particles,sigma)
        fusion_rate = 0.5*(rho**2)*sig_v

        #Add to dataframe
        df = df.append({'Potential (keV)': potential[k],'Cathode Radius (m)': cathode_radius[k],
                'Anode Radius (m)': anode_radius[k],'Fusion Rate (m^2/sec)': fusion_rate}, ignore_index=True)
    df.to_csv(f"{sim}.csv",header=True,index=False)
    del df
    print(f"Completed Processing {sim}.")




# print(sigma)

# bin_array = np.linspace(0,10,100)*1e07
# plt.figure(2)
# n,bin_edges = np.histogram(speed,bins=bin_array,density=True)
# # n,bins,patches=plt.hist(speed,bins=bin_array,density=True, edgecolor='black')
# x_scatter = bin_edges[:-1]+ 0.5*(bin_edges[1:] - bin_edges[:-1])
# plt.scatter(x_scatter, n, marker='x', c='red', s=40, alpha=1)
























