import configparser

config = configparser.ConfigParser()
config.read("Config.ini")
config.sections()

def file(section,variable):
    return float(config[section][variable])



#Constants
EPS0 = file('Constant','Epsilon')       #permittivity of free space
QE = file('Constant','ElementaryCharge')         #elementary charge
kb = file('Constant','BoltzmanConstant')                #boltzmann constan t
m_ion = file('Constant','IonMass')          #ion mass (deuterium)


#Physical Chamber Dimensions
R_Chamber = file('Chamber','Radius')      # [m]
H_Chamber = file('Chamber', 'Height')     # [m]
r_wire = file('Chamber','wireRadius')/2 # Radius of 20 guage wire in m

#Source Dimension
R_Sourse = file('Source','R')
