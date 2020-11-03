from lib import get_params,get_discrete_values 
import numpy as np
from lib import XtoL
from random import random


class Parameters: 
    #Inputs
    time_steps = get_params('Input','TimeStep',type='int')
    max_particles = get_params('Input','MaximumParticles',type='int')
    #Constants
    EPS0 = 8.854e-12  # Permittivity of Free Space
    QE = 1.602e-19 # Elementary Charge
    kb = 1  # Boltzmann Constant
    AtomicMassUnit = 1.661e-27 
    massIon = 3.334e-27  # Deuterium Ion Mass 
    Flux = 4.6e20 # Flux of ions entering [ions per second]
    marcoCharge = 1 #macroparticle charge

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
    n0 = 4.6e13                 #electron background self.DENsity in #/m^3
    phi0 = 0                    #reference potential 
    Ti = 0.1                    #ion velocity in eV (not used yet)
    vth = np.sqrt(2*QE*Ti/massIon)   #thermal velocity with Ti in eV
    
    #calculate plasma parameters
    #lD = np.sqrt(EPS0*Te/(n0*QE))      #Debye length (not used yet)
    vth = np.sqrt(2*QE*Ti/massIon)        #Thermal velocity with Ti in eV (not used yet)
    debye_length = 0.004   # [m] MANUAL SETTING FOR TESTING
    dt = 1e-10 # time step size, at vmax move 0.1dx
  


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
        self.pos = np.zeros([Parameters.max_particles,2])
        self.vel = np.zeros([Parameters.max_particles,2])
        self.insert = 400 #insert 2 ions per anode cell
    
    def get_specificWeight(self):
        specific_weight = (Parameters.Flux*Parameters.dt)/self.insert
        return specific_weight
    
    def generate_particles(self):

        for _ in range(self.insert):
            pass 
        

        


 
class ESPIC(Parameters):

    def __init__(self,cathode,anode,nodes,iterations,cathode_potential):
        self.iters = iterations
        self.cathode = cathode 
        self.anode = anode 
        self.counter = 0
        self.nx,self.ny = nodes
        self.cathode_potential = cathode_potential
        
        self.cathode_index = np.unique(self.cathode, axis=1)
        self.index_offset_x = np.round(self.nx/2).astype(int)
        self.index_offset_y = np.round(self.ny/2).astype(int)

        self.EFX = np.zeros(nodes) # Electric Field Matrix, x-component
        self.EFY = np.zeros(nodes) # Electric Field Matrix, y-component
        self.DEN = np.zeros(nodes) # Number density Matrix
    
    def build_potential_matrix(self,PHI):
        PHI[self.cathode_index[0] + self.index_offset_x, self.cathode_index[1] + self.index_offset_y] = self.cathode_potential
        PHI[self.anode[0] + self.index_offset_x, self.anode[1] + self.index_offset_y] = self.phi0

    def get_potential(self,PHI,Te,dh):
        for _ in range(self.iters):
            old_PHI = PHI
            for i in range(1, self.nx-2):
                for j in range(1, self.ny-2):
                    ni = self.DEN[i,j] # number of ions
                    rho = self.QE*(ni - self.n0*np.exp((old_PHI[i,j] - self.phi0)/(self.kb*Te))) / self.EPS0
                    chrg = -rho*dh**2
                    PHI[i,j] = (chrg - PHI[i+1, j] - PHI[i-1, j] - PHI[i, j+1] - PHI[i, j-1])/(-4)
            self.build_potential_matrix(PHI)

        return PHI

    
    def get_electricField(self,PHI,dx,dy):
        self.EFX[1:self.nx-2,:] = PHI[0:self.nx-3,:] - PHI[2:self.nx-1,:]   
        self.EFY[:,1:self.ny-2] = PHI[:,0:self.ny-3] - PHI[:,2:self.ny-1]   
        self.EFX[0,:] = 2*(PHI[0,:] - PHI[1,:])             
        self.EFX[self.nx-1,:] = 2*(PHI[self.nx-2,:] - PHI[self.nx-1,:])     
        self.EFY[:,0] = 2*(PHI[:,0] - PHI[:,1])             
        self.EFY[:,self.ny-1] = 2*(PHI[:,self.ny-2] - PHI[:,self.ny-1])     
        self.EFX = self.EFX / (dy**2)
    
    def get_density(self,CHARGE,specific_weight,dh):
        self.DEN = specific_weight*Parameters.marcoCharge*CHARGE / (dh**2)
        self.DEN[0,:] = 2*self.DEN[0,:]              # Double self.DENsity since only half volume contributing
        self.DEN[self.nx-1,:] = 2*self.DEN[self.nx-1,:]    
        self.DEN[:,0] = 2*self.DEN[:,0]
        self.DEN[:,self.ny-1] = 2*self.DEN[:,self.ny-1]
        self.DEN = self.DEN + 1e3                    # Add self.DENsity floor to help the solver

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
    
