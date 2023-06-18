# disutils: language=c++
from dimod.binary.binary_quadratic_model import BinaryQuadraticModel
def hc_sample(bqm: BinaryQuadraticModel):
    print("Hello from cython!")
    return 1