######
# Utilities to parse data from suite files. Typical format is:
"""
ode=D^2 + 3*t/(t^2 - 4)*D + 3/4/(t^2 - 4)
init=[
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ],
[ 0, 0, 0, 0, 1/4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
]
path=[ 0, 1 ]
label=(1,2)
loop_position=-1
singular_locus=[1]
"""

def parse_suite_file(filename):
    """
    Parses the output of a period suite file. Returns a dictionary of the form
    {identifier:value}. Value will be a string, with the intention to evaluate
    later.
    """

    id_dic = {}
    
    with open(filename) as F:
        for line in F:
            s = line.rstrip().split('=')

            # Detect if an identifier occurs on the line
            if not len(s) == 1:
                some_id  = s[0].strip()
                id_dic[some_id] = s[1]
            else:
                id_dic[some_id] += s[0]
    return id_dic
