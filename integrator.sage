#######################################################
#### Set this value to the directory of period-suite ##
#######################################################
pathToSuite="/Users/sertoez/git-projects/period-suite/"
#######################################################

import time
from ore_algebra import *
DOP, t, D = DifferentialOperators()
# The following may be overridden in input file, but need not be specified there
#
#
precision=100; # will be OVERRIDDEN by user specified precision, under normal usage
# if increasing this precision manually, do not forget to recompute initial conditions to higher precision.
#
#
#
precisionList=[]; # if left empty, every ODE will be integrated using the same precision
paths=[]; #unless overridden we will convert this to a list of straight paths

t0=time.time();
load(pathToSuite+"current.sage")
t1=time.time();
print "Loaded file in", t1-t0, "seconds";
steps=len(allODEs);
cohomologies=[];
allTransMats=[1..steps];

if len(paths)==0:
    paths=[[0,1] for j in [1..steps]]
    
if len(precisionList)==0:
    precisionList=[10^-precision for j in [1..steps]]


t0=time.time();


print "Integrating systems", 1, "through", steps;
for j in [1..steps]:
    noOfdeqs=len(allODEs[j-1]);
    deqs=[deq.numerator() for deq in allODEs[j-1][0:noOfdeqs]];
    print "\nIntegrating system number", j;
    transitionMatricies = []
    for i in [1..noOfdeqs]:
        tt0=time.time()
        tm=deqs[i-1].numerical_transition_matrix(paths[j-1], precisionList[j-1],assume_analytic=true)
        print "\tODE number", i, "took", time.time()-tt0, "seconds to integrate"
        transitionMatricies=transitionMatricies+[tm]
    allTransMats[j-1]=transitionMatricies;
    print "Maximal error in transition matricies:"
    print [max(tm.apply_map(lambda x : x.diameter()).list()) for tm in transitionMatricies]
    print "The largest error:"
    print max([max(tm.apply_map(lambda x : x.diameter()).list()) for tm in transitionMatricies])
    cohomology=[1..noOfdeqs];
    for i in [1..noOfdeqs]:
        if j==1:
            init=Matrix(inits[i-1]);
            a = init.nrows(); b = init.ncols();
            init =  MatrixSpace(ComplexBallField(2000), a, b)(init);
        else:
            init=Matrix(rawInits[j-2][i-1])*cohomologies[j-2];
        transitionMatrix=transitionMatricies[i-1];
        cohomology[i-1]=(transitionMatrix * init).row(0);
    cohomologies=cohomologies+[Matrix(cohomology)];
    
periodsWithError=Matrix(cohomologies[-1].rows()[0:noOfdeqs]);
maximalError=max(periodsWithError.apply_map(lambda x : x.diameter()).list());
periods=periodsWithError.apply_map(lambda x : x.mid());
print "\nAccumulated maximal error:", maximalError
t1=time.time();
print "\nIntegration took", t1-t0, "seconds in total.\n"
#print Matrix(periods).str()


# Write the periods to file
    
outputFile = open(pathToSuite+"lastPeriods",'w')

precision=-maximalError.log(10).round()
outputFile.write(str(precision)+"\n")

digits=ceil(precision/100)*100;
outputFile.write(str(digits)+"\n")

bits=ceil(log(10^digits)/log(2))+10
Rr=RealField(bits)
def print_number(num):
    re=Rr(num.real())
    im=Rr(num.imag())
    re_str=abs(re).str(digits=digits,no_sci=2)+"p"+str(digits)
    im_str=abs(im).str(digits=digits,no_sci=2)+"p"+str(digits)+"*I"
    if re.is_square():
        num_str=re_str
    else:
        num_str="-"+re_str
    if im.is_square():
        num_str=num_str+"+"+im_str
    else:
        num_str=num_str+"-"+im_str
    return num_str

numrows=periods.nrows()
numcols=periods.ncols()

prefix="Matrix(CC," +str(numrows) +","+str(numcols)+", [CC|"

numels=numrows*numcols
flat_periods=periods.list()
mat=print_number(flat_periods[0])

for i in [1..(numels-1)]:
    mat=mat+","+print_number(flat_periods[i])
    
outputFile.write(prefix+mat+"])")
outputFile.close()

