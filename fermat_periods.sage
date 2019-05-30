# the environment calling this file should have the following global variables
# pathToSuite : directory of current file
# precision : target precision
# d : degree of fermat hypersurface

storage_dir=pathToSuite+"fermat_data/";
bit_precision=ceil(log(10^precision)/log(2))+50
xi=N(exp(2*pi*I/d),bit_precision) #d-th root of unity

def gamma_values_to_compute(alphas):
    input_values={alphas[0][0]/d}
    for alpha in alphas:
        for a in alpha[0:len(alpha)-1]:
            input_values.add(a/d)
        input_values.add(sum(alpha[0:len(alpha)-1])/d)
    return input_values

@parallel(ncpus=ncpus)
def gamma_value(input_value):
    return N(gamma(input_value),bit_precision)

gamma_values={}
def compute_all_gamma_values(alphas):
    input_values=gamma_values_to_compute(alphas)
    #parallelization
    global gamma_values
    for value in gamma_value([i for i in input_values]):
        gamma_values[value[0][0][0]]=value[-1]
    
@parallel(ncpus=ncpus)
def minus_dth_root(val):
    return N(val^(-1/d),bit_precision)

scalars=[]
def compute_inverse_roots(fermat_type):
    vals=fermat_type[0:-1]
    vals.append(-fermat_type[-1])
    vals=seq(set(vals)) #remove duplicates
    root_values={}
    for r in minus_dth_root(vals):
        root_values[r[0][0][0]]=r[-1]
    global scalars
    scalars=[root_values[i] for i in fermat_type[0:-1]]
    scalars.append(root_values[-fermat_type[-1]])


def alpha_beta(alpha,beta):
    l=sum(alpha)/d
    scalar0=prod([(1-alpha[-1]/(j*d)) for j in [1..(l-1)]])
    scalar1=prod([scalars[i]^alpha[i] for i in [0..(len(alpha)-1)]])
    scalar2=prod([(1-xi^-a)/d for a in alpha[0:-1]])
    scalar3=xi^sum([alpha[i]*beta[i] for i in [0..(len(alpha)-1)]])
    scalar4=prod(gamma_values[a/d] for a in alpha[0:-1])/gamma_values[sum([a/d for a in alpha[0:-1]])]
    return scalar0*scalar1*scalar2*scalar3*scalar4
    
def periods_of_fermat(fermat_type):
    n=len(fermat_type)-2
    # load alphas and betas
    load(storage_dir+"fermat-"+str(n)+"-"+str(d)+".sage")
    compute_all_gamma_values(alphas)
    compute_inverse_roots(fermat_type)
    return Matrix([[alpha_beta(a,b) for b in betas] for a in alphas])
