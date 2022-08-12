
from ore_algebra import *
from sage.rings.integer_ring import * 

class InitialValueProblem:
    def __init__(self,ode,inits,path,precision,name):
        """Set-up the holomorphic continuation of a vector valued function whose
        entries are all annihilated by a single ordinary differential equation."""
        self.ode = ode # ordinary differential equation
        self.inits = inits # matrix of initial values, each column stores initial values of a complex valued function
        self.path = path # list of points describing a rectiliniear path in the complex plane for holomorphic continuation, inits are given at the first point
        self.precision = IntegerRing()(precision) # the number of digits of requested precision
        self.name = name # for sorting different IVPs

    # this process should be parallelizable
    def compute_transition_matrix(self):
        self.transition_matrix = self.ode.numerical_transition_matrix(self.path, 10**(-self.precision), assume_analytic=True)
        print("ODE", self.name, "is complete. Max error: ", max(self.transition_matrix.apply_map(lambda x : x.diameter()).list()))

