import os, sys, getopt

############################################################
# Retrieve ivpdir and timeout from the command line.

# Default options
myargv = sys.argv[1:]

# Parse the input configuration.
opts, args = getopt.getopt(myargv, "", ["ivpdir=", "timeout="])

# Check to make sure nothing bad happened.
if not args == []:
    print("ERROR: options misinterpreted as arguments. Please check the input.")
    sys.exit(1)

for opt, arg in opts:
    if opt == "--timeout":
        timeout = eval(arg)
    
    elif opt == "--ivpdir":
        ivpdir = arg

    else:
        print("ERROR: Invalid option: {}".format(opt))
        sys.exit(1)


############################################################
# Timeout handling.
#
# If `timeout` is set via command line, the process will self-destruct after the set time.

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

        
############################################################
# Begin main script.

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

digit_precision = precision
bit_precision   = ceil(log(10^(precision+10))/log(2))+100
field           = ComplexBallField(bit_precision)

load(pathToSuite + "arb_matrix_cereal_wrap.sage")

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
singular_locus=[1]

"""

def common_voronoi_path(ivps):
    R = PolynomialRing(Rationals(), "t")
    t = R.gens()[0]

    locals_dic = {'t' : t}
    polys = []
    for fname in ivps:
        str_poly_list = parse_suite_file(fname)['singular_locus']
        polys += sage_eval(str_poly_list, locals=locals_dic)

    path, sing_pts = voronoi_path(polys)
    return path

# For the sake of code clarity, we pass the value back as a proper function call.
# Alternatively, the rows can be written to files to be dealt with later instead.
#
# TODO: Create Voronoi path first, then pass to function.
from parse_suite_data import *
def integrate_ode(ode_label, path):
    print("Integration started")    

    #Code below replaces load(ode_label)

    DOP, t, D  = DifferentialOperators()
    locals_dic = {'D':D, 't':t}
    id_dic = parse_suite_file(ode_label)

    # File specified by ode_label defines values:
    #      ode, init, path, label, loop_position, singular_locus
    #
    # Bind identifiers read from file.
    ode  = sage_eval(id_dic['ode'], locals=locals_dic)
    init = sage_eval(id_dic['init'], locals=locals_dic)
    label = sage_eval(id_dic['label'], locals=locals_dic)

    #singular_locus = sage_eval(id_dic['singular_locus'], locals=locals_dic)

    #####
    # Do the main integration.
    
    initial_conditions = Matrix(init)
    #path, singpts = voronoi_path(singular_locus)
    
    transition_mat = complex_numerical_transition_matrix(ode, path, precision)

    # Harmonize base rings.
    if not is_exact_ring(initial_conditions.base_ring()):
        initial_conditions = initial_conditions.change_ring(transition_mat.base_ring())
        
    # Status update.
    max_err = max(x.diameter() for x in transition_mat.list())
    print "\tODE {:8} is complete. Max error: {}".format(label, max_err)
    
    # due to a bug with the Arb-Sage interface, convert to a portable object.
    transition_row = ARBMatrixCerealWrap(matrix(transition_mat.row(0)*initial_conditions))
        
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

    common_V_path = common_voronoi_path(ivps)
    
    results = map((lambda x : integrate_ode(x, common_V_path)), ivps)
    for solution in results:
        label      = solution[-1]
        tms[label] = solution[0]

    print "Integration completed in",time.time()-t0,"seconds."
    

def ith_compatible_matrix(i):
    """Function to align the cohomology basis between steps in the homotopy. """

    # File contains `change_coordinates`
    load(basis_change_files[i])

    # Might need some magic to parse/bind identifiers correctly.
    #with open(basis_change_files[i]) as F:
    #    identifier = sage_eval(F.read())
    
    # Complex ball fields will interpret `field([a,b])` as an interval. We ned to parse
    # the element correctly first before coercing it into the complex ball field.

    # Filter the IVP solutions by label, and ensure they are in the correct order.
    # These form the rows of the transition matrix.
    entries = [tms[k].entries_as_arbs() for k in sorted(tms.keys()) if k[0] == i+1]
    tm =  matrix(entries)    
    return change_coordinates*tm

## Write to file.
print "Rearranging the matrices. Writing to file..."

with open(ivpdir+"transition_mat.sobj",'w') as outfile:
    total_transition_mat = prod(ith_compatible_matrix(i) for i in range(steps))
    pickle.dump(ARBMatrixCerealWrap(total_transition_mat), outfile)

    #TODO: Also save the digit_precision somewhere sensible.

exit()
