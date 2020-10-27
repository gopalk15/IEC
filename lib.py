import configparser

class Config:

    def get_params(section,variable):
        ''' Returns Paramenters Stored in .cfg file '''
        config = configparser.ConfigParser()
        config.read("simulation_parameters.cfg")
        config.sections()
        return float(config[section][variable])
    
    #Constants
    EPS0 = 8.854e-12  # Permittivity of Free Space
    QE = 1.602e-19 # Elementary Charge
    kb = 1  # Boltzmann Constant
    AtomicMassUnit = 1.661e-27 
    massIon = 3.334e-27  # Deuterium Ion Mass 

    #Parameters
    R_Chamber = get_params('Chamber','Radius')      # [m]
    H_Chamber = get_params('Chamber', 'Height')     # [m]
    r_wire = get_params('Chamber','wireRadius')/2 # Radius of 20 guage wire in m

    #Source Dimension
    R_Sourse = get_params('Source','R')
    

class Domain(Config):
    def __init__(self):
        pass

class Fusion:
    def __init__(self):
        pass 
    
    def get_potential(self):
        pass 

    def sample_source(self):
        pass 
    
    def get_cross_section(self):
        pass

