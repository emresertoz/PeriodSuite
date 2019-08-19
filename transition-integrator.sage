pathToSuite="/usr/people/avinash/Gauss-Manin/PeriodSuite/";
load(pathToSuite+"ivpdir.sage")


## Because whoever wrote the boiler plate for the ARB library did it wrong.


print("Beginning integration...")
load(ivpdir+"meta.sage")

### FORMAT OF META.SAGE FILE:
# steps
# precision
# d
# fermat_type
# reduce


import time
import pickle
import multiprocessing as mp
from ore_algebra import *

DOP, t, D = DifferentialOperators()

bit_precision=ceil(log(10^(precision+10))/log(2))+100
field=ComplexBallField(bit_precision)

class ARBMatrixCerealWrap:
    """
    A wrapper class to enable serialization of complex arb matrix objects. The original arb matrix
    can be constructed.
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

The return type of this function is a tuple [ TMconv, Bool, label]. Here

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

## This decorator is evil...
#@parallel
def integrate_ode(ode_label):
    print("Integration started")    

    # Load the file with the ode_label.
    load(ode_label)
    # File defines values:
    #
    # ode
    # init
    # path
    # label
    # loop_position

    # Form matrix of initial conditions.
    M = Matrix(init)

    
    # At this point, we should decide on
    # -- integration path (via voronoi method),
    
    tm = ode.numerical_transition_matrix(path, 10^(-precision), assume_analytic=true)

    max_err = max( x.diameter() for x in tm.list() )    
    print "\tODE", label, "is complete. Max error: ", max_err
    
    # Check if M has been computed approximately or exactly.
    # We basically always receive exact input.
    if not (M.base_ring() == Rationals() or M.base_ring() == Integers()):
        M = MatrixSpace(ComplexBallField(tm.parent().base().precision()), M.nrows(), M.ncols())(M)

    # Compute the correct transition row from the ODE and initial conditions.
    # due to a bug with the Arb-Sage interface, we have to convert to a portable object.
    TM = matrix(tm.row(0)*M)    
    TMcomp = ARBMatrixCerealWrap(TM)

    
    # write the row to a file
    #with open(ode_label + '-solution','w') as sol_file:
    #    #sol_file.write( str( ARBMatrixCerealWrap(TM).arb_entries )
    #    #pickle.dump(TMcomp, sol_file)
        
    return [TMcomp,False,label]



########################################################################################

## Conversion file for pickling complex_ball objects.

## Save

# Takes a complex ball matrix and converts it into a pair of matrices `[A,B]`, with `A` the midpoints
# of the ball and `B` the errors of the ball.
# def convert(cbf_matrix):
#     return [cbf_matrix.apply_map(lambda x : x.mid()), cbf_matrix.apply_map(lambda x : x.diameter())]


########################################################################################

## Format output and save.
#  Takes an unconverted complex ball matrix. Converts it and writes to a file.
# def output_to_file(compl_ball_mat,filename):
#     output_file = open(ivpdir+filename,'w')
#     maximal_error=max(compl_ball_mat.apply_map(lambda x : x.diameter()).list());
#     compl_ball_mat_mid=compl_ball_mat.apply_map(lambda x : x.mid());
#     print "Accumulated maximal error:", maximal_error
#     if maximal_error == 0:
#         attained_precision=precision
#     else:
#         attained_precision=-maximal_error.log(10).round()

#     # Write the precision and digits to the first two lines.
#     output_file.write(str(attained_precision)+"\n")
#     digits=ceil(attained_precision*1.2);
#     output_file.write(str(digits)+"\n")

#     # Write the rest of the matrix.
#     # Force the matrix to print entry by entry.
#     t0=time.time()
#     print "Writing the complex_ball_matrix to file."
#     numrows = compl_ball_mat_mid.nrows()
#     numcols = compl_ball_mat_mid.ncols()
#     for i in [1..numrows]:
#         output_file.write(str(compl_ball_mat_mid[i-1].list()))
#         if i < numrows: output_file.write("\n")
#     output_file.close()


###################################################################################################
    

## Transition matrices but without intermediate base changes
tms={}

t0=time.time()
ivps=[]
for file in os.listdir(ivpdir):
    if file.startswith("IVP-") and file.endswith(".sage"):
        ivps.append(os.path.join(ivpdir,file))
        
ivps.sort()

#######################
# Parallelization by hand.

if __name__ == '__main__':

    pool = mp.Pool(mp.cpu_count() - 10) # In case someone is using the cores.

    ## most time consuming part: where the solutions are actually tracked.
    results = pool.map(integrate_ode, ivps)
    
    for solution in results:
        label      = solution[-1]
        tms[label] = solution[0]

print "Integration completed in",time.time()-t0,"seconds."


t0=time.time()
keys=tms.keys()
compatible_tms=[1..steps]

#####
## Helper function for "construct matrix with error" Turns a pair of matrices into a complex ball matrix.
# def convert_to_matrix_with_error(matrix,error):
#     new_matrix=MatrixSpace(field,matrix.nrows(),matrix.ncols())(matrix)
#     for i in [1..matrix.nrows()]:
#         for j in [1..matrix.ncols()]:            
#             new_matrix[i-1,j-1]=new_matrix[i-1,j-1].add_error(error[i-1,j-1]) 
#     return new_matrix

# ## Essentially the inverse to "convert", but operating on the entire data structure at the same time.
# def construct_matrix(keys, dictionary):
#     # given an ordered list of keys, get the corresponding rows in dictionary 
#     # and form the mxn matrix
#     data = [dictionary[key] for key in keys]
#     mat=Matrix([dat[0] for dat in data])
#     err=Matrix([dat[1] for dat in data])
#     return convert_to_matrix_with_error(mat,err)

##### END FUN BLOCK

## Transition matrices made compatible with base changes
#  The transition matrices are computed from MAGMA.
basis_change_files=[]
for file in os.listdir(ivpdir):
    if file.startswith("BaseChange-") and file.endswith(".sage"):
        basis_change_files.append(os.path.join(ivpdir,file))
basis_change_files.sort()


## Function for patching together data from the parallelization.
print "rearranging the matrices"
for i in [1..steps]:
    #nrows = max([k[1] for k in keys if k[0] == i])
    #tm_with_error = construct_matrix([(i,j) for j in [1..nrows]],tms)
    tm = matrix( field, [tms[k].list() for k in keys if k[0] == i] )
    load(basis_change_files[i-1])
    compatible_tms[i-1]=change_coordinates*tm


## Some formatting of these complex ball matrices required.
with open(ivpdir+"transition_mat",'w') as  outfile:
    total_transition_mat = prod(compatible_tms)
    pickle.dump( ARBMatrixCerealWrap(total_transition_mat), outfile )


output_to_file( prod(compatible_tms) , "transition_mat-raw")
