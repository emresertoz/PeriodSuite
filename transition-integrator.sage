pathToSuite="/usr/people/avinash/Gauss-Manin/PeriodSuite/";

# Magma source has been changed so that the Temp-directory is sent via the System call.
# This enables us to specify the load directory from a sage application, as well as prevent
# a concurrency issue.
""" Example statement: ivpdir="/usr/people/avinash/Gauss-Manin/PeriodSuite/ode_storage/test/" """

### TIMEOUT SUPPORT ###
#
# If `timeout` is set via command line, the process will self-destruct after the set time.

import sys
import signal

class Alarm(Exception):
    pass

def alarm_handler(signum, frame):
    sys.exit(1)

try:
    signal.signal(signal.SIGALRM, alarm_handler)
    signal.alarm(timeout)
except NameError:
    pass


load(pathToSuite + "voronoi_path.sage")

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

DOP, t, D       = DifferentialOperators()
digit_precision = precision
bit_precision   = ceil(log(10^(precision+10))/log(2))+100
field           = ComplexBallField(bit_precision)

load(pathToSuite + "arb_matrix_cereal_wrap.sage")

# class ARBMatrixCerealWrap:
#     """
#     A wrapper class to enable serialization of complex arb matrix objects. The original arb matrix
#     can be constructed via the `arb_matrix` method.
#     """
#     def __init__(self, arb_mat):
#         self.nrows = arb_mat.nrows()
#         self.ncols = arb_mat.ncols()
#         self.arb_entries = [ [x.mid(), x.diameter()] for x in arb_mat.list()]
#         self.base_ring = arb_matrix.base_ring()

#     def ball_field_elem(self, x):
#         return field(x[0]).add_error(x[1])
    
#     def arb_matrix(self):
#         return matrix(self.field, self.nrows , self.ncols, map(self.ball_field_elem, self.arb_entries)  )
        
#     def list(self):
#         return self.arb_entries

"""
An ode label is of the form (step_number, equation_number). The (i,j)-th equation is the $j$-th ode
of the i-th step.

The function integrates this equation and determines the functional for how the initial conditions
change in terms of the original initial conditions.

The return type of this function is a tuple [ TMconv, Bool, label], where:

-- Bool is a variable that is always False, indicating the number of loops.
-- TMconv is the deconstructed form of a complex ball matrix. A tuple [A,B] where A is the matrix
   of midpoints and B is the matrix of errors.

-- label is the dictionary label for the ode.

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
sngular_locus=[1]

"""
# For the sake of code clarity, we pass the value back as a proper function call.
# Alternatively, the rows can be written to files to be dealt with later instead.
#
def integrate_ode(ode_label):
    print("Integration started")    

    P.<t>=PolynomialRing(QQ)
    load(ode_label)    # File defines values: ode, init, path, label, loop_position, singular_locus
    
    initial_conditions = Matrix(init)
    path, singpts = voronoi_path(singular_locus)

    path = [0,1];
    
    transition_mat = complex_numerical_transition_matrix(ode, path, precision)

    # Harmonize base rings.
    if not is_exact_ring(initial_conditions.base_ring()):
        initial_conditions = initial_conditions.change_ring(transition_mat.base_ring())

    with open("dummy_file_name-"+ode_label[-10:-5], "w") as FF:
        FF.write(str(transition_mat.row(0)*initial_conditions))

        
    # Status update.
    max_err = max( x.diameter() for x in transition_mat.list() )    
    print "\tODE {:8} is complete. Max error: {}".format(label, max_err)
    
    # due to a bug with the Arb-Sage interface, convert to a portable object.
    transition_row = ARBMatrixCerealWrap(matrix( transition_mat.row(0)*initial_conditions ))
        
    return [transition_row,False,label] 
#####


""" Formatted return of the ODE solver. """
def complex_numerical_transition_matrix(ode, path, precision):
    tm = ode.numerical_transition_matrix(path, 10^(-precision), assume_analytic=true)
    return  tm.change_ring(ComplexBallField(tm.base_ring().precision()))

def is_exact_ring(ring):
    return (ring == Rationals() or ring == Integers() )

###########################################################################################
# Main script start.

## Look for load files in the directory.
#  The IVPs and the basis compatibility matrices are computed from MAGMA.
t0=time.time()
ivps=[]
basis_change_files=[]
for file in os.listdir(ivpdir):
    if file.startswith("IVP-") and file.endswith(".sage"):
        ivps.append(os.path.join(ivpdir,file))
    elif file.startswith("BaseChange-") and file.endswith(".sage"):
        basis_change_files.append(os.path.join(ivpdir,file))

ivps.sort()
basis_change_files.sort()

if len(ivps) == 0:
    print("No IVPs found in directory: {}.".format(ivpdir))
    sys.exit(1)


## The most time consuming part: where the solutions are actually tracked.
if __name__ == '__main__':    
    tms={}    
    #pool = mp.Pool(mp.cpu_count() - 10)
    #results = pool.map(integrate_ode, ivps)

    results = map(integrate_ode, ivps)
    for solution in results:
        label      = solution[-1]
        tms[label] = solution[0]

    print "Integration completed in",time.time()-t0,"seconds."
    

def ith_compatible_matrix(i):
    # File contains `change_coordinates`
    load(basis_change_files[i])
    tm =  matrix( field, [tms[k].list() for k in tms.keys() if k[0] == i+1] )
    return change_coordinates*tm

## Write to file.
print "Rearranging the matrices. Writing to file..."

with open(ivpdir+"transition_mat.sobj",'w') as outfile:
    total_transition_mat = prod( ith_compatible_matrix(i) for i in range(steps) )

    print ARBMatrixCerealWrap(total_transition_mat).arb_entries
    pickle.dump( ARBMatrixCerealWrap(total_transition_mat), outfile )

    #TODO: Also save the digit_precision somewhere sensible.

exit()
