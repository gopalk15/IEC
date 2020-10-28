from lib import get_params,get_discrete_values 
import numpy as np

class Parameters: 
    #Inputs
    time_steps = get_params('Input','TimeStep')  
    #Constants
    EPS0 = 8.854e-12  # Permittivity of Free Space
    QE = 1.602e-19 # Elementary Charge
    kb = 1  # Boltzmann Constant
    AtomicMassUnit = 1.661e-27 
    massIon = 3.334e-27  # Deuterium Ion Mass 

    # Physical Grid Dimensions
    cathode_gt = get_params('Grid','Cathode_gt')   # Geometric transparency cathode (manual for now)
    anode_gt = get_params('Grid','Anode_gt')       # Geometric transparency anode   (manual for now)

    #Physical Chamber Dimensions
    chamber_radius = get_params('Chamber','Radius')      # [m]
    chamber_height = get_params('Chamber', 'Height')     # [m]
    wire_radius = get_params('Chamber','wireRadius')/2 # Radius of 20 guage wire in m

    #Source Dimension
    R_Sourse = get_params('Source','R')

    #input settings
    n0 = 4.6e13                 #electron background density in #/m^3
    phi0 = 0                    #reference potential 
    Ti = 0.1                    #ion velocity in eV (not used yet)
    vth = np.sqrt(2*QE*Ti/massIon)   #thermal velocity with Ti in eV
    
    #calculate plasma parameters
    #lD = np.sqrt(EPS0*Te/(n0*QE))      #Debye length (not used yet)
    vth = np.sqrt(2*QE*Ti/massIon)        #Thermal velocity with Ti in eV (not used yet)
    debye_length = 0.004   # [m] MANUAL SETTING FOR TESTING
  


class Domain(Parameters):
    def __init__(self):
        self.dx = self.debye_length # [m] cell-size in x & y direction 
        self.dy = self.debye_length # [m] cell-size in x & y direction 

    def build_chamber(self):
        x_axis = self.chamber_radius/2
        y_axis = self.chamber_height/2

        wall_y = np.arange(-y_axis,y_axis,self.dy)
        wall_x = np.ones(len(wall_y))*x_axis
        
        ceiling_x = np.arange(-x_axis,x_axis,self.dx)
        ceiling_y = wall_x = np.ones(len(ceiling_x))*y_axis

        left = (-wall_x,wall_y)
        top = (ceiling_x,ceiling_y)
        right = (wall_x,wall_y)
        bottom = (ceiling_x,-ceiling_y)

        chamber = (left,top,right,bottom)

        return chamber

    def build_grid(self,radius):
        theta = np.arange(0,2*np.pi,np.pi/1000)
        x_array = radius*np.cos(theta)
        y_array = radius*np.sin(theta)
        X = get_discrete_values(x_array, self.dx) #distrised cathode array in x-direction
        Y = get_discrete_values(y_array, self.dy) #distrised cathode array in y-direction
        grid = (X,Y)
        return grid
       
    def get_nodes(self):
        nx = int(np.floor((2*self.chamber_radius)/self.dx))
        ny = int(np.floor(self.chamber_height/self.dy))
        nodes = (nx,ny)
        return nodes
    



class Fusion:
    def __init__(self):
        pass 
    
    def get_potential(self):
        pass 

    def sample_source(self):
        pass 
    
    def get_cross_section(self):
        pass
