import sys 
sys.path.insert(1,'C:\\Users\\gopal\\git\\IEC\\fusor')

import pytest
from lib import sampleIsotropicVel as sample_vel

@pytest.mark.usefixtures("particles")
class TestGenerateParticles:
    def test_method(self,particles):
        print(particles)
        assert particles
        
    def test_sample_vel(self):
        result = sample_vel(300)
        print(len(result))

    def test_generate_particles(self,particles):
        radius,angle = particles.get_spray_values()
        assert particles.count == 0
        particles.generate(radius,angle)
        # assert particles.pos[300,0] != 0.0
        assert particles.count == 400
        particles.generate(radius,angle)
        assert particles.count == 800

       





