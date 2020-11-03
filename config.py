from lib import get_params,get_discrete_values 
import numpy as np
from lib import XtoL

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
    source_radius = get_params('Source','R')

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
        self.top_counter = 0 
        self.bottom_counter = 0
        self.right_counter = 0 
        self.left_counter = 0 
        self.anode_counter = 0
        self.cathode_counter = 0 
    
    def build_chamber(self):
        x_axis = self.chamber_radius/2
        y_axis = self.chamber_height

        wall_y = np.arange(0,y_axis,self.dy)
        wall_x = np.ones(len(wall_y))*x_axis
    
        ceiling_x = np.arange(-x_axis,x_axis,self.dx)
        ceiling_y = np.ones(len(ceiling_x))*y_axis

        left = (-wall_x,wall_y)
        top = (ceiling_x,ceiling_y)
        right = (wall_x,wall_y)
        bottom = (ceiling_x,ceiling_y*0)
 
        chamber = (left,top,right,bottom)


        return chamber

    def build_grid(self,radius):
        theta = np.arange(0,2*np.pi,np.pi/1000)
        x_array = radius*np.cos(theta) 
        y_array = radius*np.sin(theta)  # offset in y-direction
        X = get_discrete_values(x_array, self.dx) #distrised cathode array in x-direction
        Y = get_discrete_values(y_array, self.dy) #distrised cathode array in y-direction
        x_index = list(map(int,X/self.dx))
        y_index = list(map(int,Y/self.dy))
        
        return [x_index,y_index]
    
    def get_nodes(self):
        ''' Returns the number of nodes in the x and y direction'''
        nx = int(np.floor((2*self.chamber_radius)/self.dx))
        ny = int(np.floor(self.chamber_height/self.dy))
        nodes = [nx,ny]
        return nodes
    
    

class Particles(Parameters):
    def __init__(self,nodes):
        self.cell = Parameters.debye_length
        self.count = 0 #Number of particles
        self.nx,self.ny = nodes
    
    def insert_particles(self,particles_per_cell):
        np = (self.ny - 1)*particles_per_cell #insert n particles per cell


class ESPIC(Parameters):

    def __init__(self,cathode,anode,nodes,iterations,cathode_potential):
        self.iters = iterations
        self.cathode = cathode 
        self.anode = anode 
        self.col_counter = 0
        self.nx,self.ny = nodes
        self.cathode_potential = cathode_potential
        
        self.cathode_index = np.unique(self.cathode, axis=1)
        self.index_offset_x = np.round(self.nx/2).astype(int)
        self.index_offset_y = np.round(self.ny/2).astype(int)
    
    def build_potential_matrix(self,PHI):
        PHI[self.cathode_index[0] + self.index_offset_x, self.cathode_index[1] + self.index_offset_y] = self.cathode_potential
        PHI[self.anode[0] + self.index_offset_x, self.anode[1] + self.index_offset_y] = self.phi0

    def get_potential(self,PHI,density):
        for _ in range(inters):
            old_PHI = PHI
            for i in range(1, self.nx-2):
                for j in range(1, self.ny-2):
                    ni = den[i,j] # number of ions
                    rho = self.QE*(ni - self.n0*np.exp((old_PHI[i,j] - self.phi0)/(self.kb*self.Te))) / self.EPS0
                    chrg = -rho*self.dx**2
                    PHI[i,j] = (chrg - PHI[i+1, j] - PHI[i-1, j] - PHI[i, j+1] - PHI[i, j-1])/(-4)
            self.build_potential_matrix(PHI)

        return PHI


    def sample_source(self):
        pass 
    
    def get_cross_section(self):
        pass
    

    
class Fusion(Parameters):
    
    def __init__(self):
        self.counter = 0 
        self.x_position = np.array([])
        self.y_position = np.array([])
        self.time = np.array([])
    
