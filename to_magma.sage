output_file = open(ivpdir+"periods",'w')

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
#prefix="Matrix(CC," +str(1) +","+str(numcols)+", "
for i in [1..numrows]:
  output_file.write(str(periods_mid[i-1].list()))
  if i < numrows: output_file.write("\n")

output_file.close()
print "Periods written to file in", time.time()-t0,"seconds."
