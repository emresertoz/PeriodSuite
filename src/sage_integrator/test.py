## This script is called by specifying the location of the initial value problems (IVPs = ODE + initial conditions etc.) like so:
##               sage integrator.sage --ivpdir="path/to/suite/ode_storage/incinerator/"

## Parse the input configuration.
#opts, _ = getopt.getopt(sys.argv[1:], "", ["ivpdir=", "timeout="])
#for opt, arg in opts:
#    if opt == "--timeout":
#        timeout = eval(arg)
#    elif opt == "--ivpdir":
#        ivpdir = arg
#    else:
#        print("ERROR: Invalid option: {}".format(opt))
#        sys.exit(1)
#############################################################

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('ivpdir')
print(parser.parse_args().ivpdir)
