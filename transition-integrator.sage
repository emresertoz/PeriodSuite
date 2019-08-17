pathToSuite="/usr/people/avinash/Gauss-Manin/PeriodSuite/";
load(pathToSuite+"ivpdir.sage")
ncpus=100

print("Beginning integration...")
load(ivpdir+"meta.sage")

import time
from ore_algebra import *
DOP, t, D = DifferentialOperators()

bit_precision=ceil(log(10^(precision+10))/log(2))+100
field=ComplexBallField(bit_precision)

@parallel(ncpus=ncpus)
def integrate_ode(ode_label):
    load(ode_label)
    if loop_position > -1:
        return integrate_ode_with_loop(ode_label)
    tm=ode.numerical_transition_matrix(path, 10^(-precision), assume_analytic=true)
    print "\tODE", label, "is complete. Max error: ", max(tm.apply_map(lambda x : x.diameter()).list())
    M=Matrix(init)
    if not (M.base_ring() == Rationals() or M.base_ring() == Integers()):
        M=MatrixSpace(ComplexBallField(tm.parent().base().precision()),M.nrows(),M.ncols())(M)
    TM=tm.row(0)*M
    return [convert(TM),False,label]

def integrate_ode_with_loop(ode_label):
    load(ode_label)
    path1=path[0:loop_position]
    path2=path[loop_position-1:-1]
    ## not super efficient!
    #path3=[path[loop_position],path[len(path)-1]]
    #path3=[path[0],path[-1]]
    path3=path
    tm1=ode.numerical_transition_matrix(path1, 10^(-precision), assume_analytic=true)
    #print "\tODE", label, "is completed until loop. Max error: ", max(tm1.apply_map(lambda x : x.diameter()).list())
    tm2=ode.numerical_transition_matrix(path2, 10^(-precision), assume_analytic=true)
    #print "\tODE", label, "loop completed. Max error: ", max(tm2.apply_map(lambda x : x.diameter()).list())
    tm3=ode.numerical_transition_matrix(path3, 10^(-precision), assume_analytic=true)
    max_err=max([max(tm.apply_map(lambda x : x.diameter()).list()) for tm in [tm1,tm2,tm3]])
    print "\tODE", label, "loop and limit completed. Max error: ", max_err
    M=Matrix(init)
    if not (M.base_ring() == Rationals() or M.base_ring() == Integers()):
        M=MatrixSpace(ComplexBallField(tm1.parent().base().precision()),M.nrows(),M.ncols())(M)
    TM1=tm1.row(0)*M
    TM2=(tm2.row(0)*tm1)*M
    #TM3=tm3*tm1*M 
    TM3=tm3*M 
    return [convert(TM1),convert(TM2),convert(TM3),ode,True,label]
    
def convert(cbf_matrix):
    return [cbf_matrix.apply_map(lambda x : x.mid()),cbf_matrix.apply_map(lambda x : x.diameter())]

def convert_to_matrix_with_error(matrix,error):
    new_matrix=MatrixSpace(field,matrix.nrows(),matrix.ncols())(matrix)
    for i in [1..matrix.nrows()]:
        for j in [1..matrix.ncols()]:
            new_matrix[i-1,j-1]=new_matrix[i-1,j-1].add_error(error[i-1,j-1]) 
    return new_matrix

def construct_matrix(keys, dictionary):
    # given an ordered list of keys, get the corresponding rows in dictionary 
    # and form the mxn matrix
    data = [dictionary[key] for key in keys]
    mat=Matrix([dat[0] for dat in data])
    err=Matrix([dat[1] for dat in data])
    return convert_to_matrix_with_error(mat,err)

def compute_periods_of_fermat():
    print "Computing periods of Fermat"
    t0=time.time()
    load(pathToSuite+"fermat_periods.sage")
    # fermat_type and the degree d of fermat are loaded from the file "meta.sage"
    fermat_period_matrix=periods_of_fermat(fermat_type)
    print "Fermat periods computed in", time.time() -t0, "seconds."
    fpm_rows=fermat_period_matrix.nrows()
    fpm_cols=fermat_period_matrix.ncols()
    return MatrixSpace(field,fpm_rows,fpm_cols)(fermat_period_matrix);

def output_to_file(periods,filename):
    output_file = open(ivpdir+filename,'w')
    maximal_error=max(periods.apply_map(lambda x : x.diameter()).list());
    periods_mid=periods.apply_map(lambda x : x.mid());
    print "Accumulated maximal error:", maximal_error
    if maximal_error == 0:
        attained_precision=precision
    else:
        attained_precision=-maximal_error.log(10).round()
    output_file.write(str(attained_precision)+"\n")
    digits=ceil(attained_precision*1.2);
    output_file.write(str(digits)+"\n")
    t0=time.time()
    print "Writing the periods to file."
    numrows=periods_mid.nrows()
    numcols=periods_mid.ncols()
    for i in [1..numrows]:
        output_file.write(str(periods_mid[i-1].list()))
        if i < numrows: output_file.write("\n")
    output_file.close()



## Transition matrices but without intermediate base changes
tms={}
## the transition matrices around loops, if any
loops={}
## in case of loop, the final transition matrix and odes
limit_tm={}
limit_periods={}
final_odes={}

t0=time.time()
ivps=[]
for file in os.listdir(ivpdir):
    if file.startswith("IVP-") and file.endswith(".sage"):
        ivps.append(os.path.join(ivpdir,file))
ivps.sort()
## most time consuming part
for solution in integrate_ode(ivps):
    has_loop=solution[-1][-2]
    label=solution[-1][-1]
    if has_loop:
        tms[label]=solution[-1][0]
        loops[label]=solution[-1][1]
        limit_tm[label]=solution[-1][2]
        final_odes[label]=solution[-1][3]
    else:
        tms[label]=solution[-1][0]
print "Integration completed in",time.time()-t0,"seconds."
loop_keys=loops.keys()
loop_keys.sort()

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


## Function for patching together data from the parallelization.
print "rearranging the matrices"
for i in [1..steps]:
    nrows=max([k[1] for k in keys if k[0] == i])
    #ncols=len(tms[(i,1)][0])
    tm_with_error=construct_matrix([(i,j) for j in [1..nrows]],tms)
    load(base_change_files[i-1])
    compatible_tms[i-1]=change_coordinates*tm_with_error


with open(ivpdir+filename,'w') as  outfile:
    outfile.write(tm_with_error)
    outfile.write(compatible_tms)
    

