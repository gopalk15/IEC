import numpy as np

''' __Helper Functions___'''

def get_cross_section(speed): 
    ''' Takes Velocity array'''
    massIon = 3.334e-27  # Deuterium Ion Mass 
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

    cross = []
    for v in speed:
        E = (1/2)*massIon*(v**2)*6.242e15 #covert J to KeV
        term1 = A5 + A2/((A4 - A3*E)**2 + 1)
        term2 = E*(np.exp(A1/np.sqrt(E)) - 1)
        term3 = AA5 + AA2/((AA4 - AA3*E)**2 + 1)
        term4 = E*(np.exp(AA1/np.sqrt(E)) - 1)
        cs = ((term1/term2) + (term3/term4))*1e-28
        cross.append(cs)
    
    return np.array(cross)

def calc_density(position):
    chamber_radius = 0.2 # [m]
    area = np.pi*(chamber_radius**2)

    density = len(position)/area

    return density


def gather_density(position,DEN):
    ''' Needs postion array and density matrix to be defined'''
    #Chamber sizes
    chamber_height = 0.5 # [m]
    chamber_radius = 0.2 # [m]
    dh = 0.004 # [m]

    density = 0
    for index in range(len(position)):
        fi = (position[index, 0] + chamber_radius-dh)/dh    #real i index of particle's cell
        i = np.floor(fi).astype(int)          #integral part
        hx = fi - i                              #the remainder
        fj = (position[index,1] + (chamber_height/2)-dh)/dh #real i index of particle's cell
        j = np.floor(fj).astype(int)            #integral part
        hy = fj - j 


        rho = DEN[i,j]*(1-hx)*(1-hy)
        rho = rho + DEN[i+1,j]*hx*(1-hy)
        rho = rho + DEN[i,j+1]*(1-hx)*hy
        rho = rho + DEN[i+1,j+1]*hx*hy

        density += rho

    return density.astype(float)

def get_speed(velocity):
    return np.hypot(velocity[:,0],velocity[:,1])

def get_mean_KE(speeds):
    '''Returns the  total kinetic energy of all particles in scaled units '''
    speed_sum_sq = 0
    for speed in speeds:
        square = speed**2
        speed_sum_sq += square

    KE = 0.5*massIon*speed_sum_sq
    
    return KE/(len(speeds) + 1)

def get_vel_average(x_array,y_array,sigma):
    ''' Intergartes where {x_array,y_array,sigma} are n sized arrays '''
    average_velocity = 0
    for i in range(len(x_array)-1):
        intergrand = sigma[i]*y_array[i]*(x_array[i+1]-x_array[i])
        average_velocity += intergrand
    
    return average_velocity
