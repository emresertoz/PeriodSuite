# import sage.all
# import os, sys, getopt
from sage.functions.other import ceil



def output_to_file(periods,filepath):
    output_file = open(filepath,'w')
    maximal_error=max(periods.apply_map(lambda x : x.diameter()).list());
    periods_mid=periods.apply_map(lambda x : x.mid());
    print("Accumulated maximal error:", maximal_error)
    attained_precision=-maximal_error.log(10).round()
    output_file.write(str(attained_precision)+"\n")
    digits=ceil(attained_precision*1.2);
    output_file.write(str(digits)+"\n")
    print("Writing the periods to file.")
    numrows=periods_mid.nrows()
    numcols=periods_mid.ncols()
    for i in range(1,numrows+1):
        output_file.write(str(periods_mid[i-1].list()))
        if i < numrows: output_file.write("\n")
    output_file.close()


