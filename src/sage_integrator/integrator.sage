import os, time
import argparse
import multiprocessing
from ore_algebra import *

# my functions
from pathToSuite import *
import input_output as io
from fermat_periods import FermatPeriods
from initial_value_problem import *

# This script must be called by specifying the location of the initial value problems (IVPs = ODE + initial conditions etc.) like so:
#               sage integrator.sage "path/to/suite/ode_storage/incinerator/"
# We retrieve the IVP directory here and assign it to ivpdir. 
parser = argparse.ArgumentParser()
parser.add_argument('ivpdir')
ivpdir=parser.parse_args().ivpdir

####### Read-in files #######
# meta.sage: stores global information, e.g., precision, fermat_type, reduce 
load(ivpdir+"meta.sage")
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
family_indices = set(sorted(i for i,j in labels))

# voronoi paths
# FIXME: get this out of here
# TODO: Ideally, if there are no (non-apparent) singularities along [0,1] we should take the straigh path
from voronoi_path import *
P=PolynomialRing(Rationals(),'x')
for ii in family_indices:
    ivp_batch = [ivp for ivp in ivps if ivp.label[0] == ii]
    polys = [ivp.ode.monic().denominator() for ivp in ivp_batch]
    polys = [P(p) for p in polys if p != 1]
    path = voronoi_path(polys)
    for ivp in ivp_batch:
        ivp.path = path

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
for solution in holo_cont(ivps):
    label=solution[-1][-1]
    rows_of_period_tms[label]=reassemble_cbf_matrix(field,solution[-1][-2])
print("Integration completed in",time.time()-t0,"seconds.")
print("Now putting it all together ...")
# reconstruct the period transition matrices from the rows
period_tms=[]
for ii in family_indices:
    row_nums = sorted([j for i,j in labels if i == ii])
    period_tms.append(Matrix([rows_of_period_tms[(ii,j)][0] for j in row_nums]))

# multiply everything together to express final periods in terms of periods of initial hypersurface (=Fermat)
# reverse the order in which these appear in the lists
base_changes.reverse(); period_tms.reverse()
master_transition_matrix=prod(B*P for B,P in zip(base_changes,period_tms))
periods_of_fermat=FermatPeriods(degree_of_hypersurface,fermat_type,precision).approximate_period_matrix.change_ring(field)
# periods: this contains the periods of the target hypersurface
periods=master_transition_matrix*periods_of_fermat
# writing the period matrix to file
io.output_to_file(periods,ivpdir+"periods")

