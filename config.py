from lib import get_params,get_discrete_values,sampleIsotropicVel,fusion_cross_section 
import numpy as np
from random import uniform
import warnings

class Parameters: 
    #Inputs
    time_steps = get_params('Input','TimeStep',type='int')
    max_particles = get_params('Input','MaximumParticles',type='int')
    #Constants
    EPS0 = 8.854e-12  # Permittivity of Free Space
    QE = 1.602e-19 # Elementary Charge
    kb = 1  # [J/K]Boltzmann Constant
    AtomicMassUnit = 1.661e-27 
    massIon = 3.334e-27  # Deuterium Ion Mass 
    
    marcoCharge = 1 #macroparticle charge
    deturonDiameter = 2.1413e-15


    # Physical Grid Dimensions
    cathode_gt = get_params('Grid','Cathode_gt')   # Geometric transparency cathode (manual for now)
    anode_gt = get_params('Grid','Anode_gt')       # Geometric transparency anode   (manual for now)

    #Physical Chamber Dimensions
    chamber_radius = get_params('Chamber','Radius')      # [m]
    chamber_height = get_params('Chamber', 'Height')     # [m]
    wire_radius = get_params('Chamber','wireRadius')/2 # Radius of 20 guage wire in m
    chamber_pressure = get_params('Chamber','Pressure') # [Torr]

    #Source Dimension
    source_radius = get_params('Source','sourceRadius')

    #input settings
    n0 = 4.6e13                 #electron background self.DENsity in #/m^3
    phi0 = 0                    #reference potential 
    Ti = 0.1                    #ion velocity in eV (not used yet)
    vth = np.sqrt(2*QE*Ti/massIon)   #thermal velocity with Ti in eV
    
    #calculate plasma parameters
    #debye_length = np.sqrt(EPS0*Te/(n0*QE))      #Debye length (not used yet)
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
        self.lost = 0 # Number of Particles lost
        self.nx,self.ny = nodes
        self.pos = np.zeros([Parameters.max_particles,2])
        self.vel = np.zeros([Parameters.max_particles,2])
        self.insert = 600 #insert 2 ions per anode cell
    
    def get_specificWeight(self,cathode_potential):
        Flux = 1.6e21*np.abs(cathode_potential)/100000 # Flux of ions entering [ions per second]
        specific_weight = (Flux*self.dt)/self.insert
        return specific_weight
    
    def get_spray_values(self):
        ''' Returns values for spray arch'''
        angle = np.tan((2/3)*self.chamber_radius/self.chamber_height)
        radius = np.sqrt((self.chamber_radius**2) + (self.chamber_height**2)) - 0.1

        return (radius,angle)

    def generate(self,radius):

       # largest angle == pi/7
        
        # xv = np.zeros([self.insert])
        # yv = np.zeros([self.insert])
        # for i in range(self.insert):
        #     theta = np.random.rand(1)*2*np.pi   # Generate random polar angle
        #     x = self.source_radius*np.cos(theta)                # Get x position
        #     y = self.source_radius*np.sin(theta)                # Get y position
        #     xv[i] = x
        #     yv[i] = y

        # self.pos[self.count:self.count+self.insert,0] = xv
        # self.pos[self.count:self.count+self.insert,0] = yv

        # pt1 = np.random.rand(self.insert)
        # pt2 = np.random.rand(self.insert)
        # pt3 = np.random.rand(self.insert)
        # pt11 = np.random.rand(self.insert)
        # pt12 = np.random.rand(self.insert)
        # pt13 = np.random.rand(self.insert)
        # self.vel[self.count:self.count+self.insert,0] = (-1.5 + pt1 + pt2 + pt3)*self.vth     # x velocity
        # self.vel[self.count:self.count+self.insert,1] = (-1.5 + pt11 + pt12 + pt13)*self.vth  # y velcoity

        for i in range(self.insert):
            theta = np.random.rand(1)*2*np.pi
            self.pos[i + self.count,0] = radius*np.cos(theta)
            self.pos[i + self.count,1] = radius*np.sin(theta) 
            self.vel[i + self.count,0],self.vel[i+self.count,1] = sampleIsotropicVel(self.vth)

        self.count += self.insert

    
    def get_index(self,index,dh):
        fi = (self.pos[index, 0] + self.chamber_radius-dh)/dh    #real i index of particle's cell
        i = np.floor(fi)          #integral part
        hx = fi - i                              #the remainder
        fj = (self.pos[index,1] + (self.chamber_height/2)-dh)/dh #real i index of particle's cell
        j = np.floor(fj)              #integral part
        hy = fj - j 

        return np.array([i,j,hx,hy])

    def get_fusion_prob(self,index,density,specific_weight):
        vx = self.vel[index,0]
        vy = self.vel[index,1]
        del_vx = vx*self.dt*1e4 
        del_vy = vy*self.dt*1e4
        path_length = np.sqrt((del_vx**2)+(del_vy**2))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sigma = fusion_cross_section(vx,vy,self.massIon)
        probability = density*sigma*path_length*specific_weight
        return probability
    
    def kill(self,index):
        self.pos[[index,self.count-1]] = self.pos[[self.count-1,index]]
        self.vel[[index,self.count-1]] = self.vel[[self.count -1 ,index]]
        #Rest last particle to 0
        self.pos[[self.count -1]] = np.array([0,0])
        self.vel[[self.count-1]] = np.array([0,0])
        self.count += -1
        self.lost += 1

 
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
        self.EFX = self.EFX / (dx**2)
        self.EFY = self.EFY/ (dy**2)
    
    def get_density(self,CHARGE,specific_weight,dh):
        self.DEN = specific_weight*Parameters.marcoCharge*CHARGE / (dh**2)
        self.DEN[0,:] = 2*self.DEN[0,:]              # Double self.DENsity since only half volume contributing
        self.DEN[self.nx-1,:] = 2*self.DEN[self.nx-1,:]    
        self.DEN[:,0] = 2*self.DEN[:,0]
        self.DEN[:,self.ny-1] = 2*self.DEN[:,self.ny-1]
        self.DEN = self.DEN + 1e3                    # Add self.DENsity floor to help the solver

    
    def gather_electricField(self,i,j,hx,hy):

        E = np.array([self.EFX[i,j],self.EFY[i,j]])*(1-hx)*(1-hy)
        E = E + np.array([self.EFX[i+1,j], self.EFY[i+1,j]])*hx*(1-hy)            #(i+1,j)
        E = E + np.array([self.EFX[i,j+1], self.EFY[i,j+1]])*(1-hx)*hy            #(i,j+1)
        E = E + np.array([self.EFX[i+1,j+1], self.EFY[i+1,j+1]])*hx*hy

        return E

    def gather_density(self,i,j,hx,hy):
        rho = self.DEN[i,j]*(1-hx)*(1-hy)
        rho = rho + self.DEN[i+1,j]*hx*(1-hy)
        rho = rho + self.DEN[i,j+1]*(1-hx)*hy
        rho = rho + self.DEN[i+1,j+1]*hx*hy
        return rho
    
class Fusion(Parameters):
    
    def __init__(self):
        self.events = 0 
        self.position = []
        self.rate_data = []
        self.counter = 0
    
    def occured(self,fuse_position):
        self.events += 1
        self.counter += 1
        self.position.append(fuse_position)
        

    
