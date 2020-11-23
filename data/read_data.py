import h5py 
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform


#Global Parameters 

dt = 1e-10 #[s]
chamber_height = 0.5 # [m]
chamber_radius = 0.2 # [m]
dh = 0.004 # [m]
num_insert = 600 # number of particles inserted per timestep
massIon = 3.334e-27  # Deuterium Ion Mass 
rscale = 1/chamber_radius
R = 8.314 #J/Kmol (Gas Constant)
M = 4.0282404 #g/mol (molar mass of deuturim)

#Take the mean speed of the root-mean-square velocity of deuterium at 10 keV
sbar = np.sqrt((3*R*(10e03*1.160e04))/M) * (rscale/dt)




sim_file = 'potential'
data_set = 5


class Histogram:
    """A class to draw a Matplotlib histogram as a collection of Patches."""

    def __init__(self, data, xmax, nbars, density=False):
        """Initialize the histogram from the data and requested bins."""
        self.nbars = nbars
        self.density = density
        self.bins = np.linspace(0, xmax, nbars)
        self.hist, bins = np.histogram(data, self.bins, density=density)

        # Drawing the histogram with Matplotlib patches owes a lot to
        # https://matplotlib.org/3.1.1/gallery/animation/animated_histogram.html
        # Get the corners of the rectangles for the histogram.
        self.left = np.array(bins[:-1])
        self.right = np.array(bins[1:])
        self.bottom = np.zeros(len(self.left))
        self.top = self.bottom + self.hist
        nrects = len(self.left)
        self.nverts = nrects * 5
        self.verts = np.zeros((self.nverts, 2))
        self.verts[0::5, 0] = self.left
        self.verts[0::5, 1] = self.bottom
        self.verts[1::5, 0] = self.left
        self.verts[1::5, 1] = self.top
        self.verts[2::5, 0] = self.right
        self.verts[2::5, 1] = self.top
        self.verts[3::5, 0] = self.right
        self.verts[3::5, 1] = self.bottom

    def draw(self, ax):
        """Draw the histogram by adding appropriate patches to Axes ax."""
        codes = np.ones(self.nverts, int) * path.Path.LINETO
        codes[0::5] = path.Path.MOVETO
        codes[4::5] = path.Path.CLOSEPOLY
        barpath = path.Path(self.verts, codes)
        self.patch = patches.PathPatch(barpath, fc='tab:green', ec='k',
                                  lw=0.5, alpha=0.5)
        ax.add_patch(self.patch)

    def update(self, data):
        """Update the rectangle vertices using a new histogram from data."""
        self.hist, bins = np.histogram(data, self.bins, density=self.density)
        self.top = self.bottom + self.hist
        self.verts[1::5, 1] = self.top
        self.verts[2::5, 1] = self.top





def get_index(index,dh):
    fi = (self.pos[index, 0] + self.chamber_radius-dh)/dh    #real i index of particle's cell
    i = np.floor(fi)          #integral part
    hx = fi - i                              #the remainder
    fj = (self.pos[index,1] + (self.chamber_height/2)-dh)/dh #real i index of particle's cell
    j = np.floor(fj)              #integral part
    hy = fj - j 

    return np.array([i,j,hx,hy])

def get_cross_section(velocity): 
    ''' Takes Velocity array'''
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
    for v in velocity:
        E = (1/2)*massIon*(v**2)*6.242e15 #covert J to KeV
        term1 = A5 + A2/((A4 - A3*E)**2 + 1)
        term2 = E*(np.exp(A1/np.sqrt(E)) - 1)
        term3 = AA5 + AA2/((AA4 - AA3*E)**2 + 1)
        term4 = E*(np.exp(AA1/np.sqrt(E)) - 1)
        cs = ((term1/term2) + (term3/term4))*1e-28
        cross.append(cs)
    
    return np.array(cross)


def get_density(index):
    ''' Needs postion array and density matrix to me defined'''
    fi = (pos[index, 0] + chamber_radius-dh)/dh    #real i index of particle's cell
    i = np.floor(fi).astype(int)          #integral part
    hx = fi - i                              #the remainder
    fj = (pos[index,1] + (chamber_height/2)-dh)/dh #real i index of particle's cell
    j = np.floor(fj).astype(int)            #integral part
    hy = fj - j 


    rho = DEN[i,j]*(1-hx)*(1-hy)
    rho = rho + DEN[i+1,j]*hx*(1-hy)
    rho = rho + DEN[i,j+1]*(1-hx)*hy
    rho = rho + DEN[i+1,j+1]*hx*hy
    return rho

def get_specificWeight(cathode_potential):
    Flux = 1.6e21*np.abs(cathode_potential)/100000 # Flux of ions entering [ions per second]
    specific_weight = (Flux*dt)/num_insert
    return specific_weight






# Getting Relavant Data

DEN = get_density_data(sim_file,data_set)
speed = get_particle_velocity(sim_file,data_set)
pos = get_particle_position(sim_file,data_set)
cathode_potential = -100e03



spwt = get_specificWeight(cathode_potential)
temp = np.abs(cathode_potential)
time = np.arange(0,len(speed))*dt

speed_hist = Histogram(speed, 2 * sbar, 50, density=True)

mean_KE = get_mean_KE(speed)
a = massIon / 2 / mean_KE
print(mean_KE)

DPI = 100
width, height = 1000, 500
fig = plt.figure(figsize=(width/DPI, height/DPI), dpi=DPI)
speed_ax = fig.add_subplot(122)
speed_ax.set_xlabel('Speed $v\,/m\,s^{-1}$')
speed_ax.set_ylabel('$f(v)$')

# TODO don't hardcode the upper limit for the histogram speed axis.
ticks = np.linspace(0, 600, 7, dtype=int)
speed_ax.set_xticks(ticks * rscale/dt)
speed_ax.set_xticklabels([str(tick) for tick in ticks])
speed_ax.set_yticks([])

fig.tight_layout()


# Maximum value of the 2D Maxwell-Boltzmann speed distribution.
fmax = np.sqrt(massIon / mean_KE / np.e)
speed_ax.set_ylim(0, fmax)

'''
Use a high-resolution grid of speed points so that the exact distribution
looks smooth.
'''
sgrid_hi = np.linspace(0, speed_hist.bins[-1], 200)
f = 2 * a * sgrid_hi * np.exp(-a * sgrid_hi**2)
mb_line, = speed_ax.plot(sgrid_hi, f, c='0.7')
# Maximum value of the 2D Maxwell-Boltzmann speed distribution.
fmax = np.sqrt(massIon / mean_KE / np.e)
speed_ax.set_ylim(0, fmax)

# For the distribution derived by averaging, take the abcissa speed points from
# the centre of the histogram bars.
sgrid = (speed_hist.bins[1:] + speed_hist.bins[:-1]) / 2
mb_est_line, = speed_ax.plot([], [], c='r')
mb_est = np.zeros(len(sgrid))

# A text label indicating the time and step number for each animation frame.
xlabel, ylabel = sgrid[-1] / 2, 0.8 * fmax
label = speed_ax.text(xlabel, ylabel, '$t$ = {:.1f}s, step = {:d}'.format(0, 0))


plt.show()

# #Calculate Maxwellian Distribution
# dist = []
# for i in range(len(pos)):
#     n = get_density(i)
#     max_dist = (n*(massIon/(2*np.pi*temp))**(3/2))*np.exp(-massIon*(vel[i]**4)/(2*temp))
#     dist.append(max_dist)

# maxwell = np.array(dist)
# sigma = get_cross_section(vel)
# print(np.multiply(sigma,maxwell))
# speed = np.linspace(0,2.5e06,len(dist))





# plt.scatter(speed,dist)
# plt.show()