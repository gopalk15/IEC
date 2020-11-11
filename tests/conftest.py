import sys
sys.path.insert(1,'C:\\Users\\gopal\\git\\IEC\\fusor')

import pytest 
import config

''' Setup Nodes '''
fusor = config.Domain()
nx,ny = nodes = fusor.get_nodes()


@pytest.fixture(scope="session")
def particles():
    return config.Particles(nodes)
