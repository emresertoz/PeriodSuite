pathToSuite="/path/to/suite/";
ivpdir=pathToSuite+"incinerator/";
ncpus=100

print("Beginning integration...")
load(ivpdir+"meta.sage")

import time
from ore_algebra import *
DOP, t, D = DifferentialOperators()

@parallel(ncpus=ncpus)
def integrate_ode(ode_label):
    load(ode_label)
    tm=ode.numerical_transition_matrix(path, 10^(-precision), assume_analytic=true)
    print "\tODE", label, "is complete. Max error: ", max(tm.apply_map(lambda x : x.diameter()).list())
    M=Matrix(init)
    if not (M.base_ring() == Rationals() or M.base_ring() == Integers()):
      M=MatrixSpace(tm.base_ring(),M.nrows(),M.ncols())(M)
    TM=tm.row(0)*M
    return [[x.mid() for x in TM],[x.diameter() for x in TM],label]

field=ComplexBallField(round(log(10^(precision+10))/log(2)))
def convert_to_matrix_with_error(matrix,error):
    new_matrix=MatrixSpace(field,matrix.nrows(),matrix.ncols())(matrix)
    for i in [1..matrix.nrows()]:
      for j in [1..matrix.ncols()]:
        new_matrix[i-1,j-1]=new_matrix[i-1,j-1].add_error(error[i-1,j-1]) 
    return new_matrix

## Transition matrices but without intermediate base changes
tms={}

t0=time.time()
ivps=[]
for file in os.listdir(ivpdir):
  if file.startswith("IVP-") and file.endswith(".sage"):
    ivps.append(os.path.join(ivpdir,file))
ivps.sort()
## most time consuming part
for solution in integrate_ode(ivps):
    tms[solution[-1][-1]]=solution[-1][0:2]
print "Integration completed in",time.time()-t0,"seconds."

## Transition matrices made compatible with base changes
base_change_files=[]
for file in os.listdir(ivpdir):
  if file.startswith("BaseChange-") and file.endswith(".sage"):
    base_change_files.append(os.path.join(ivpdir,file))
base_change_files.sort()

t0=time.time()
keys=tms.keys()
steps=len(base_change_files)
compatible_tms=[1..steps]

print "rearranging the matrices"
for i in [1..steps]:
    nrows=max([k[1] for k in keys if k[0] == i])
    ncols=len(tms[(i,1)][0])
    tm=Matrix([tms[(i,j)][0] for j in [1..nrows]])
    err=Matrix([tms[(i,j)][1] for j in [1..nrows]])
    tm_with_error=convert_to_matrix_with_error(tm,err)
    load(base_change_files[i-1])
    compatible_tms[i-1]=change_coordinates*tm_with_error
    
compatible_tms.reverse()
period_tm=prod(compatible_tms)

# reduce is in the meta file designating if the initial conditions have been reduced to the periods of the Fermat hypersurface or were precomputed already
if reduce:
  print "Computing periods of Fermat"
  t0=time.time()
  load(pathToSuite+"fermat_periods.sage")
# fermat_type and the degree d of fermat are loaded from the file "meta.sage"
  fermat_period_matrix=periods_of_fermat(fermat_type)
  print "Fermat periods computed in", time.time() -t0, "seconds."
  fpm_rows=fermat_period_matrix.nrows()
  fpm_cols=fermat_period_matrix.ncols()
  fpm=MatrixSpace(field,fpm_rows,fpm_cols)(fermat_period_matrix);
  periods=period_tm*fpm
else:
  periods=period_tm

# records the period matrix for magma to read
load(pathToSuite+"to_magma.sage")

