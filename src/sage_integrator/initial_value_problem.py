from sage.all import *
from ore_algebra import *
# from sage.rings.integer_ring import * 

#TODO: make IVP take a field instead of precision? this way I don't have to create a field everytime ivp is called
class InitialValueProblem:
    def __init__(self,ode,inits,path,precision,field,label):
        """Set-up the holomorphic continuation of a vector valued function whose
        entries are all annihilated by a single ordinary differential equation."""
        self.ode = ode # ordinary differential equation
        self.inits = Matrix(inits) # matrix of initial values, each column stores initial values of a complex valued function
        self.path = path # list of points describing a rectiliniear path in the complex plane for holomorphic continuation, inits are given at the first point
        self.precision = IntegerRing()(precision) # the number of digits of requested precision
        self.field = field # complex ball field
        self.label = label # for sorting different IVPs

    # this process should be parallelizable
    def compute_transition_matrix(self):
        self.transition_matrix = self.ode.numerical_transition_matrix(self.path, 10**(-self.precision), assume_analytic=True)
        print("ODE", self.label, "is complete. Max error: ", max(self.transition_matrix.apply_map(lambda x : x.diameter()).list()))

    def holomorphic_continuation(self):
        """Return the value of the holomorphic vector valued function at the end of the path 
        by holomorphic continuation. Assumes the ODE is smooth at the end of the path point 
        (apparent singularities are also fine).
        
        If the ODE is singular at the end point, then the output has a different interpretation.
        Use asymptotic_holomorphic_continuation instead.
        """
        if not hasattr(self, 'transition_matrix'):
            self.compute_transition_matrix()
        inits=self.inits.change_ring(self.field)
        return self.transition_matrix.row(0)*inits
    
    def asymptotic_holomorphic_continuation(self):
        """
        Return the coordinates of the vector valued function at the end of the path in terms of the
        Frobenius basis for the solution space of the ODE. 
        
        The output is a matrix, each column stores the coordinate in the solution space 
        of the corresponding holomorphic coordinate function.
        """
        if not hasattr(self, 'transition_matrix'):
            self.compute_transition_matrix()
        inits=self.inits.change_ring(self.field)
        return self.transition_matrix*inits

def disassemble_cbf_matrix(cbf_matrix):
    """Take a CBF matrix and return the tuple consisting of
    the matrix of midpoints and the matrix of errors. This circumvents
    the unpickleable nature of cbf matricies."""
    return (cbf_matrix.apply_map(lambda x : x.mid()),cbf_matrix.apply_map(lambda x : x.diameter()))

def reassemble_cbf_matrix(cb_field,disassembled_cbf_matrix):
    """Reasseble a CBF matrix from the tuple of its midpoints and errors so that the base field is the given cb_field."""
    matrix = Matrix(cb_field,disassembled_cbf_matrix[0])
    errs = Matrix(disassembled_cbf_matrix[1])
    nrows = matrix.nrows()
    ncols = matrix.ncols()
    return Matrix(nrows,ncols, lambda i,j : matrix[i,j].add_error(errs[i,j]))

