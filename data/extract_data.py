import h5py 
import numpy as np 


class Data():

    def __init__(self,simulation_file,data_set):
        self.file = f"{simulation_file}{data_set}.h5"

    def clean_data(self,array):
        size = (len(array) - 1 )

        for i in range(size + 1):
            if np.all(array[size -i] == 0 ) and np.all(array[size - (i+1)] == 0):
                index = size - (i+1) 
        new_array = np.empty_like(array[:index])
        new_array = array[:index]

        return new_array
        

    def get_particle_position(self):
        with h5py.File(self.file) as hdf:
            group = hdf.get('DataSets/potential/')
            postional_array = np.array(group.get('ParticlePosition'))
        
        position = self.clean_data(postional_array)

        return position 

    def get_particle_velocity(self):
        with h5py.File(self.file) as hdf:
            group = hdf.get('DataSets/potential/')
            particle_velocity = np.array(group.get('ParticleVelocity'))
        clean_velocity = self.clean_data(particle_velocity)
        return clean_velocity
        # vx = clean_velocity[:,0]
        # vy = clean_velocity[:,1]
        # v = []
        # for x,y in zip(vx,vy):
        #     square = np.sqrt((x**2)+(y**2))
        #     v.append(square)
        
        # return v

    def get_density_data(self):
        with h5py.File(self.file) as hdf:
            group = hdf.get('DataSets/density/')
            den = np.array(group.get('Density')) 
        
        return den


    def get_electric_field(self,axis=None):
        with h5py.File(self.file) as hdf:
            group = hdf.get('DataSets/electricfield/')
            efx = np.array(group.get('electricFieldx'))
            efy = np.array(group.get('eletricFieldy'))

        
        if axis == 'X':
            return efx 
        
        elif axis == 'Y':
            return efy
        
        elif axis is None:
            print("Please an axis variable")

    


