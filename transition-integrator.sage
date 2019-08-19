pathToSuite="/usr/people/avinash/Gauss-Manin/PeriodSuite/";
load(pathToSuite+"ivpdir.sage")

import time
import pickle
import multiprocessing as mp
from ore_algebra import *

print("Beginning integration...")

### FORMAT OF META.SAGE FILE:
# steps
# precision
# d
# fermat_type
# reduce
load(ivpdir+"meta.sage")

DOP, t, D     = DifferentialOperators()
bit_precision = ceil(log(10^(precision+10))/log(2))+100
field         = ComplexBallField(bit_precision)

class ARBMatrixCerealWrap:
    """
    A wrapper class to enable serialization of complex arb matrix objects. The original arb matrix
    can be constructed via the `arb_matrix` method.
    """
    def __init__(self, arb_mat):
        self.nrows = arb_mat.nrows()
        self.ncols = arb_mat.ncols()
        self.arb_entries = [ [x.mid(), x.diameter()] for x in arb_mat.list()]

    def ball_field_elem(self, x):
        return field(x[0]).add_error(x[1])
    
    def arb_matrix(self):
        return matrix(field, self.nrows , self.ncols, map(self.ball_field_elem, self.arb_entries)  )
        
    def list(self):
        return self.arb_entries

"""
An ode label is of the form (step_number, equation_number). The (i,j)-th equation is the $j$-th ode
of the i-th step.

The function integrates this equation and determines the functional for how the initial conditions
change in terms of the original initial conditions.

The return type of this function is a tuple [ TMconv, Bool, label], where:

-- Bool is a variable that is always False, indicating the number of loops.
-- TMconv is the deconstructed form of a complex ball matrix. A tuple [A,B] where A is the matrix
   of midpoints and B is the matrix of errors.

   NOTE: The only reason for this format is that Sage's @parallel decorator is miserable. In our code
         I think we should just parallelize by hand.

-- label is the dictionary label for the ode.

Given an ode_label identifying which of the 21


====================================================
EXAMPLE FORMAT OF THE "IVP-*-*.sage" FILE:

ode=D^2 + 3*t/(t^2 - 4)*D + 3/4/(t^2 - 4)
init=[
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ],
[ 0, 0, 0, 0, 1/4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
]
path=[ 0, 1 ]
label=(1,2)
loop_position=-1

"""
def integrate_ode(ode_label):
    print("Integration started")    

    # Load the file with the ode_label.
    load(ode_label)
    # File defines values: ode, init, path, label, loop_position

    M = Matrix(init)   # The matrix of initial conditions.
    
    # TODO: At this point, we should decide on
    # -- integration path (via voronoi method),   
    tm = ode.numerical_transition_matrix(path, 10^(-precision), assume_analytic=true)

    max_err = max( x.diameter() for x in tm.list() )    
    print "\tODE", label, "is complete. Max error: ", max_err
    
    # Check if M has been computed approximately or exactly.
    # We basically always receive exact input.
    if not (M.base_ring() == Rationals() or M.base_ring() == Integers()):
        M = MatrixSpace(ComplexBallField(tm.parent().base().precision()), M.nrows(), M.ncols())(M)

    # Compute the transition row from the ODE and initial conditions.
    # due to a bug with the Arb-Sage interface, we convert to a portable object.
    TM = matrix(tm.row(0)*M)    
    TMcomp = ARBMatrixCerealWrap(TM)

    # Alternatively, the rows can be written to files to be dealt with later.
    return [TMcomp,False,label] 
    

###########################################################################################
# Main script start.

## Look for load files in the directory.
t0=time.time()
ivps=[]
for file in os.listdir(ivpdir):
    if file.startswith("IVP-") and file.endswith(".sage"):
        ivps.append(os.path.join(ivpdir,file))
        
ivps.sort()

## Transition matrices made compatible with base changes
#  The transition matrices are computed from MAGMA.
basis_change_files=[]
for file in os.listdir(ivpdir):
    if file.startswith("BaseChange-") and file.endswith(".sage"):
        basis_change_files.append(os.path.join(ivpdir,file))

# TODO: This goes badly wrong if there are 11 base-change files...
basis_change_files.sort()


if __name__ == '__main__':

    ## Transition matrices but without intermediate base changes
    tms={}    

    ## most time consuming part: where the solutions are actually tracked.
    pool = mp.Pool(mp.cpu_count() - 10)
    results = pool.map(integrate_ode, ivps)
    for solution in results:
        label      = solution[-1]
        tms[label] = solution[0]

    print "Integration completed in",time.time()-t0,"seconds."
    
## Function for patching together data from the parallelization.
print "rearranging the matrices"

def ith_compatible_matrix(i):
    # File contains `change_coordinates`
    load(basis_change_files[i])
    tm =  matrix( field, [tms[k].list() for k in tms.keys() if k[0] == i+1] )
    return change_coordinates*tm

total_transition_mat = prod( ith_compatible_matrix(i) for i in range(steps) )

## Write to file.
with open(ivpdir+"transition_mat",'w') as  outfile:
    pickle.dump( ARBMatrixCerealWrap(total_transition_mat), outfile )
