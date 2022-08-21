import os, time
import argparse
import multiprocessing
from ore_algebra import *

# my functions
from pathToSuite import *
import input_output as io
from fermat_periods import FermatPeriods
from initial_value_problem import *
from lmhs import *

# # FIXME: uncomment this 
# # This script must be called by specifying the location of the initial value problems (IVPs = ODE + initial conditions etc.) like so:
# #               sage integrator.sage "path/to/suite/ode_storage/incinerator/"
# # We retrieve the IVP directory here and assign it to ivpdir. 
# parser = argparse.ArgumentParser()
# parser.add_argument('ivpdir')
# ivpdir=parser.parse_args().ivpdir

####### Read-in files #######
# meta.sage: stores global information, e.g., degree, dimension, precision, fermat_type, reduce 
load(ivpdir+"meta.sage")
degree = d
# Get all the file paths for the individual IVPs from ivpdir
ivp_paths=[]
for file in os.listdir(ivpdir):
    if file.startswith("IVP-") and file.endswith(".sage"):
        ivp_paths.append(os.path.join(ivpdir,file))
ivp_paths.sort()
# read in the coordinate changes for the period matrices in between families
base_change_files=[]
for file in os.listdir(ivpdir):
    if file.startswith("BaseChange-") and file.endswith(".sage"):
        base_change_files.append(os.path.join(ivpdir,file))
base_change_files.sort()
base_changes=[]
for file in base_change_files:
    load(os.path.join(ivpdir,file))
    # change_coordinates: rational coordinate change matrix, to switch from a cohom basis suitable for the family, to the grevlex basis at target hypersurface in each family.
    # note that all initial conditions are recorded in terms of the grevlex basis at initial hypersurface
    base_changes.append(change_coordinates)
# setup notation to interpret the ivp files
DOP, t, D = DifferentialOperators()
bit_precision=ceil((precision+10)*log(10)/log(2))+100
field=ComplexBallField(bit_precision)
# Read the ivp files and form the IVP objects
ivps=[]
for file in ivp_paths:
    load(file)
    ivps.append(InitialValueProblem(ode,init,path,precision,field,label))
####### Reading files completed #######

labels=sorted([ivp.label for ivp in ivps])
max_index = max(i for i,j in labels)


#Get CPU count
ncpus=multiprocessing.cpu_count()
print("Beginning integration...")
t0=time.time()
@parallel(ncpus=ncpus)
def holo_cont(ivp):
    hc=ivp.holomorphic_continuation()
    return [disassemble_cbf_matrix(hc),ivp.label]
# the solution to our ivps are rows of period transition matrices
rows_of_period_tms={}

smooth_ivps = [ivp for ivp in ivps if ivp.label[0] < max_index]
for solution in holo_cont(smooth_ivps):
    label=solution[-1][-1]
    rows_of_period_tms[label]=reassemble_cbf_matrix(field,solution[-1][-2])
print("Integration completed in",time.time()-t0,"seconds.")
print("Now putting it all together ...")
# reconstruct the period transition matrices from the rows
period_tms=[]
for ii in range(1,max_index): 
#WARNING: preparing periods of penultimate hypersurface
    row_nums = sorted([j for i,j in labels if i == ii])
    period_tms.append(Matrix([rows_of_period_tms[(ii,j)][0] for j in row_nums]))

# multiply everything together to express penultimate periods in terms of periods of initial hypersurface (=Fermat)
# reverse the order in which these appear in the lists
#LMHS special: throw out the last base change
base_changes = base_changes[:-1]
base_changes.reverse(); period_tms.reverse()
master_transition_matrix=prod(B*P for B,P in zip(base_changes,period_tms))
periods_of_fermat=FermatPeriods(d,fermat_type,precision).approximate_period_matrix.change_ring(field)
# periods: this contains the periods of the target hypersurface
periods=master_transition_matrix*periods_of_fermat

## Major new part for LMHS
print("Periods of the penultimate hypersurface computed. Now we turn to the final family of hypersurfaces.")

# writing the period matrix to file
# io.output_to_file(periods,ivpdir+"periods")

final_ivps = [ivp for ivp in ivps if ivp.label[0] == max_index]
lmhs = LMHS(final_ivps,periods)

print("Final batch of integration begun...")
t0=time.time()
@parallel(ncpus=ncpus)
def limit_holo_cont(ivp):
    mats=lmhs.integrate_ivp_with_loop(ivp)
    return [[disassemble_cbf_matrix(mat) for mat in mats],ivp.label]

rows1={} # rows of the period matrix of a nearby hypersurface
rows2={} # rows of the period matrix of a nearby hypersurface, after monodromy
limit_period_rows={} # coordinates of limiting periods
for solution in limit_holo_cont(final_ivps):
    label=solution[-1][-1]
    rows1[label] = reassemble_cbf_matrix(field,solution[-1][-2][0]) 
    limit_period_rows[label] = reassemble_cbf_matrix(field,solution[-1][-2][1]) 
    rows2[label] = reassemble_cbf_matrix(field,solution[-1][-2][2]) 
print("Integration completed in",time.time()-t0,"seconds.")

print("Computing monodromy")
row_nums = sorted([j for i,j in labels if i == max_index])
P1 = Matrix([rows1[(max_index,j)][0] for j in row_nums])
P2 = Matrix([rows2[(max_index,j)][0] for j in row_nums])
monod_float = P1.inverse()*P2 # P2 = P1*monod
monod = monod_float.apply_map(lambda x : round(x.mid().real()))
monod_err = RealField(10)(max((monod_float - monod.change_ring(field)).apply_map(lambda x : abs(x.mid())).list()))
print("Error in rounding the monodromy matrix to nearest integer matrix, check that this is small: ", monod_err)

logMon,unip,mult=logarithm_of_monodromy(monod)
# dimension of hypersurface is read from meta.sage, this bounds the index of nilpotency of logMon and gives the correct degrees for the weight filtration (Schmid)
W_matrix,W_dims = weight_filtration_of_nilpotent_matrix(logMon,dimension) 
change_to_W_basis = W_matrix.change_ring(field).inverse() 

# expansions = lmhs.expansions_as_log_puiseux_series()
expansions = lmhs.expansions_as_laurent_series(mult)
limit_periods = matrix([(matrix(expansions[j-1])*limit_period_rows[(max_index,j)]*change_to_W_basis)[0] for j in row_nums])

# get tthe dimensions of the pieces of the Hodge filtration
Fp_dims = hodge_filtration_dimensions(dimension,degree)

##
# Going to compute the Plücker coordinates of the Hodge flag in the limit as t->0
##

# get rid of duplicate dimensions, 0, and the full space
fps={fp for fp in Fp_dims if fp != 0 and fp != max(Fp_dims)}


######### preparing matrix for minor computation #########
L = limit_periods;
# remove complex ball coefficients that contain zero to get honest valuations
L = clean_cbs_matrix(L)
min_row_valuations = [min(l.valuation() for l in row) for row in L.rows()]
# getting all entries to have non-negative valuation by scaling rows (i.e. cohomology classes)
R = L.base_ring(); t = R.gens()[0];
L = Matrix([t**(-v)*row for v,row in zip(min_row_valuations,L.rows())])
L = L.change_ring(R.power_series_ring())

######## Lowering CBF precision for computing valuations of minors #################
lowPrec = 100 #FIXME: this bit precision might need to be adjusted based on size of the matrix L
RR = L.base_ring().change_ring(ComplexBallField(lowPrec)); 
L_low_prec = clean_cbs_matrix(L.change_ring(RR))
####################################################################################

## TODO:
# in general, I might not want all Fp's, e.g., with a limit K3 I might only want F2.
# then having the option to skip the larger flags will be hugely beneficial
##

t0=time.time()
Fp_matrices = []
for fp in fps:
    t1=time.time()
    # considering the vector space generated by the first fp rows
    Fp_low_prec = Matrix(L_low_prec[:fp])
    # find the maximal minor with smallest valuation
    print(f"Determening the limit of a flag of dimension {fp}.")
    print("Computing valuations of all minors. If matrix is large, this will take a lot of time...")
    smallest_val, index_of_smallest_val = min_val_of_maximal_minors(Fp_low_prec)
    #### FIXME
    # the process above is way too slow. I could just invert one submatrix randomly, 
    # then try to substitute zero to newFp and check if it has full rank,
    # if not, then I can compute Pluecker coordinates of newFp, which seems much faster
    ####
    # Now pass to high precision and invert the submatrix giving the smallest valuation
    Fp = Matrix(L[:fp])
    submat = submatrix_from_column_indices(Fp,index_of_smallest_val)
    print("Inverting a minor with high precision. This may take time.")
    t2 = time.time()
    inv=invert_cbs_matrix(submat)
    print(f"Inversion completed in {time.time()-t2} seconds.")
    newFp=clean_cbs_matrix(inv*Fp)
    if sum(newFp.apply_map(lambda l : l.valuation() < 0).list()) > 0:
        print("Something is wrong. Let me know with your example.") 
        # Possibly lowPrec needs to be increased
        # go back to step 1 and increase precision
    Fp_matrices.append(newFp.apply_map(lambda l : l[0])) # l[0] is evaluation at 0
    print(f"Finding the limit of the flag of dimension {fp} is done in {time.time()-t1} seconds.")
print(f"All limit flags are computed in {time.time()-t0} seconds.")

# TODO:
# It might be preferable to have a single matrix whose first fp rows (for all p) store the p'th flag.
# to do this one could replace the first fp rows of L with newFp
# one could also use the fact that a maximal minor is now the identity

# The result will be a single CBF matrix, the Fp_dims, W's and W_dims.
# Try it on the curve and the degenerate K3, see if you can get the Picard rank
