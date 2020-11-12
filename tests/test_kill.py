import sys 
sys.path.insert(1,'C:\\Users\\gopal\\git\\IEC')

import pytest 
from lib import kill_particle

def test_kill(particles):
    radius,angle = particles.get_spray_values()
    particles.generate(radius,angle)

    a = particles.pos[[0]]
    b = particles.pos[[3]]
    print(particles.pos[:5,:])

    kill_particle(particles.pos,0,4)
    
    print(particles.pos[:5,:])

    



    





