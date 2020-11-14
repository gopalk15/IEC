from config import Domain,Particles,ESPIC,Fusion
from lib import get_params
import logging
import numpy as np
import h5py
import concurrent.futures
import matplotlib.pyplot as plt
from time import perf_counter

''' Setup Logger '''
logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)

# formatter = logging.Formatter('%(name)s :: %(process)d :: %(message)s')
logging.basicConfig(level=logging.WARNING, format='%(name)s :: %(message)s')

#Variable Inputs (for later)
cathode_radius = get_params('Grid','Cathode_Radius') # Cathode [m]
anode_radius = get_params('Grid','Anode_Radius')  # Anode [m]



def main(params):
    '''
        params = (cathode_potential,data_set,grid_ratio:None,spray_with:None)
    '''
    cathode_potential = params[0]
    data_set = params[1]

    if params[2] is None:
        grid_ratio = cathode_radius/anode_radius #Anode radius stays fixed
    else:
        grid_ratio = params[2]

    logger.info("INIT")
    Te = np.abs(cathode_potential)    #electron temperature in eV
    fusor = Domain()
    nx,ny = nodes = fusor.get_nodes()
    cathode = fusor.build_grid(anode_radius*grid_ratio)
    anode = fusor.build_grid(anode_radius)
    dx = fusor.dx
    dy = fusor.dy

    loop = ESPIC(cathode,anode,nodes,10,cathode_potential)
    particles = Particles(nodes)
    spray_radius,spray_angle = particles.get_spray_values()
    if params[3] is None:
        spray_width = spray_angle
    
    else:
        spray_width = params[3]

    fusion = Fusion()
    PHI_G = np.zeros(nodes)
    loop.build_potential_matrix(PHI_G)
    timeStep = loop.time_steps

    spec_wt = particles.get_specificWeight()

    for step in range(0,timeStep):
        logger.info(f"""
                    Begining Time Step # {step} 
                    Particle Count: {particles.count}
                        """)
        
        loop.DEN = np.zeros(nodes) # Number Density Matrix
        loop.EFX = np.zeros(nodes) # Electric Field Matrix, x-component
        loop.EFY = np.zeros(nodes) # Electric Field Matrix, y-component
        CHG = np.zeros(nodes) # Charge Density Matrix
        loop.counter = 0
        particles.lost = 0
        fusion.counter = 0
        
        ''' 1. Compute Charge Density '''
        for p in range(particles.count):
            fi = (particles.pos[p, 0] + fusor.chamber_radius-dx)/dx    #real i index of particle's cell
            i = np.floor(fi).astype(int)           #integral part
            hx = fi - i                              #the remainder
            fj = (particles.pos[p,1] + (fusor.chamber_height/2)-dx)/dx #real i index of particle's cell
            j = np.floor(fj).astype(int)              #integral part
            hy = fj - j                              #the remainder
            
            #interpolate charge to nodes
            CHG[i, j] = CHG[i, j] + (1-hx)*(1-hy)
            CHG[i+1, j] = CHG[i+1, j] + hx*(1-hy)
            CHG[i, j+1] = CHG[i, j+1] + (1-hx)*hy
            CHG[i+1, j+1] = CHG[i+1, j+1] + hx*hy

        # Calculate Densisty 
        loop.get_density(CHG,spec_wt,dx)

            
        ''' 2. Compute Electric Potential '''
        PHI_G = loop.get_potential(PHI_G,Te,dx)
        logger.info("Computed Electric Potential")



        ''' 3. Compute Electric Field '''
        logger.info('Computing Electric Field')
        loop.get_electricField(PHI_G,dx,dy)
        
        
        ''' 4. Generate Particles '''
        logger.info('Injecting Particles')
        particles.generate(spray_radius,spray_angle)

    
        ''' 5. Move Particles '''
        logger.info("Moving Particles")
        m = 0 # Particle index

        while m < particles.count:
            i,j,hx,hy = particles.get_index(m,dx)
        
            i = int(i)
            j = int(j)

    
            E_field = loop.gather_electricField(i,j,hx,hy)

            F = fusor.QE*E_field
            a = F/fusor.massIon
            particles.vel[m,:] = particles.vel[m,:] + a*loop.dt
            particles.pos[m,:] = particles.pos[m,:] + particles.vel[m,:]*loop.dt

        
            rho = loop.gather_density(i,j,hx,hy)
            fusion_prob = particles.get_fusion_prob(m,rho,spec_wt)
            fusion_chance = np.random.rand(1)*fusion_prob
            logging.debug(f"Fusion Probability: {fusion_prob}")
            logging.debug(f"Fusion Chance: {fusion_chance}")
            radial_distance = np.sqrt(particles.pos[m,0]**2 + particles.pos[m,1]**2)

            # Top Wall 
            if particles.pos[m,1] > fusor.chamber_height/2:
                particles.kill(m)  
                m += -1 
                fusor.top_counter += 1
            
            #Bottom Wall
            elif particles.pos[m,1] < -fusor.chamber_height/2:
                particles.kill(m)
                m += -1 
                fusor.bottom_counter += 1
            
            #Left Wall
            elif particles.pos[m,0] < -fusor.chamber_radius:
                particles.kill(m)
                m += -1
                fusor.left_counter += 1

            #Right Wall 
            elif particles.pos[m,0] > fusor.chamber_radius:
                particles.kill(m)
                m += -1
                fusor.right_counter += 1
            
        
            #Anode 
            elif (radial_distance < anode_radius + fusor.wire_radius) and (radial_distance > anode_radius - fusor.wire_radius):
                prob = np.random.rand(1)
                if prob > fusor.anode_gt:
                    particles.kill(m)
                    m += -1 
                    fusor.anode_counter += 1
            
            #Cathode
            elif (radial_distance < cathode_radius + fusor.wire_radius) and (radial_distance > cathode_radius - fusor.wire_radius):
                prob = np.random.rand(1)
                if prob > fusor.cathode_gt:
                    particles.kill(m)
                    m += -1 
                    fusor.cathode_counter += 1
            
            #Check Fusion
            elif fusion_chance[0] >= fusion_prob:
                fuse_position = particles.pos[[m]]
                fusion.occured(fuse_position[0])
                particles.kill(m)
                m += -1 
            
            m += 1                                      # move onto next particle
        
        #Get Fusion data
        fusion.rate_data.append([fusion.counter,step*fusor.dt])

        logger.info('Finshed moving Particles')
        logger.info(f"""
                        Net Charge: {particles.insert - particles.lost}
                        Particles lost: {particles.lost}
                        Fusion Events / dt: {fusion.counter}
                        Total Fusion Events: {fusion.events}
                    """)

    with h5py.File(f'data\\potential{data_set}.h5','w') as hdf:
        G2 = hdf.create_group("DataSets/potential/")
        dataset1 = G2.create_dataset('ParticlePosition',data=particles.pos)
        dataset2 = G2.create_dataset('ParticleVelocity', data=particles.vel)
        dataset6 = G2.create_dataset('PHI',data=PHI_G) 

        G3 = hdf.create_group("DataSets/electricfield/")
        dataset4 = G3.create_dataset('electricFieldx',data=loop.EFX)
        dataset5 = G3.create_dataset('electricFieldy',data=loop.EFY)

        G4 = hdf.create_group("DataSets/density/")
        dataset3 = G4.create_dataset('Density',data=loop.DEN)

        
        G1 = hdf.create_group(f"DataSets/fusion")
        dataset7 = G1.create_dataset('Position',data=fusion.position)
        dataset9 = G1.create_dataset('RateData',data = fusion.rate_data)
        dataset10 = G1.create_dataset('FusionCount',data = fusion.events)
        

        groups = [G1, G2, G3, G4]

        for group in groups:
            group.attrs['gridPotential'] = f'{cathode_potential/1000} kV'
            group.attrs['sprayAngle'] = f'{spray_width} radians'
            group.attrs['gridRatio'] = f'{grid_ratio}'
    
    print(f"Completed Simulation :: {data_set}")
            


if __name__ == '__main__':

    parameters = ((-10000,1,None,None),(-50000,2,None,None),(-70000,3,None,None),(-80000,4,None,None),(-90000,5,None,None),(-200000,6,None,None),(-50000,7,None,None),(-700000,8,None,None))

    start = perf_counter()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(main,parameters)



    end = perf_counter() 

    logger.warning(f'Completed in {round(end-start)/60} minutes(s) ')
   




