BEGIN { 
  FS=","
  OFS=""
  buffer=""
  abort=0
}

# if you hit a commented line
# stop processing the file
/^\#/ {abort=1; print $0; next}
abort==1 {print "#" $0; next}

# the first line is not modified
NR==1 {print $0; next}

# if a sqr bracket is opened at the end
# print the entire line
/\[$/ {print $0; next}

# if a bracket is closed at the beginning of line
# then print buffer 
/^\]/ {
  if (NF==1) str="]" 
  else str="],"
  if (buffer == "") print str
  else print buffer "\n" str
  buffer=""
  next
  }

# if you find a line ending ] then
# complete the text in buffer,
# print out all the intermediate terms with 
# commas in between
/\]$/ {
  if ((length($0)==1) && (buffer!="")) {
  print buffer "\n" $0; 
  buffer = ""
  next
  }
  if (NF==1) print buffer $1
  else {
    if (buffer != "") print buffer $1 ","
    else printf $1 ","
    for (i=2;i<=NF-1;i++) printf $i ","
    print $NF
  }
  buffer=""; next
}

# if you have a line ending with , then
# complete the text in buffer,
# print out all the intermediate terms with 
# commas in between
/,$/ {
  if (buffer != "") print buffer $1 ","
  else printf $1 ","
  for (i=2;i<=NF-2;i++) printf $i ","
  # the last non-empty entry needs a new line as well as a come
  if (NF>2) print $(NF-1) ","
  # ommit $NF which is empty
  buffer=""; next
}

# if the last symbol is a backslash
# we need to merge the last field
# with the next line (and so on) until
# a field seperator (,) is hit
/\\$/ { 
  if (NF==1) {
    buffer = buffer substr($0,1,length($0)-1)
  }
  else {
    if (buffer != "") print buffer $1 ","
    else printf $1 ","
    for (i=2;i<=NF-1;i++) printf $i ","
    buffer = substr($NF,1,length($NF)-1)
  }
      next
}

# if this last block is activated
# then all the other checks failed
# this happens only when we have field
# that need to be joined together
# no additional check is necessary
{ 
  if (NF==1) {
    buffer = buffer $0
  }
  else {
    if (buffer != "") printf buffer $1 ",\n"
    else printf $1 ","
    for (i=2;i<=NF-1;i++) printf $i ","
    buffer = $NF
  }
}

