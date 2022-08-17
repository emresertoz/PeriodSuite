from sage.all import *
from ore_algebra import *

class LMHS:
    def __init__(self,ivps,initial_period_matrix):
        self.ivps = ivps
        self.odes = [ivp.ode for ivp in ivps]
        self.initial_period_matrix = initial_period_matrix
        self.path = ivps[0].path # the paths of integration should all be the same
        self.field = ivps[0].field
        self.refine_path()
        self.mark_path_for_loop()

    def refine_path(self):
        #FIXME make the penultimate point closer to the final point
        # But I won't do it now, so this function does nothing
        print("We should refine the path in the future. For now, see
        if it is OK.")
        print(self.path)
        self.path = self.path

    def mark_path_for_loop(self):
        path = self.path
        a = path[-2]; b=path[-1]; # penultimate and ultimate points
        self.marked_path = path[:-2] + [[a],[b],[a]]

# this must be done in parallel
    def integrate_ivp_with_loop(self,ivp):
        end_point = self.path[-1]
        ode = ivp.ode
        # marked path will give a seqence of 3 transition matrices
        T=ode.numerical_transition_matrix(self.marked_path,ivp.precision,assume_analytic=True)
        T=[t.change_ring(self.field) for t in T]
        M=LocalMonodromy(ode,end_point).monodromy
        inits=self.inits.change_ring(self.field)
# take the output matrices: T0, T1, T2 and the local monodromy M
# hit each of these below with *inits* (then with the initial period matrix)
    # T1[0], T2, ((T3*T2.inverse())*M*T2)[0]
    # this gives you perids at penultimate, ultimate, and penultimate after monodromy
    # the one at ultimate needs all of the coordinates, so we keep all the matrix
# now reconstruct the period matrix at penultimate
        newT:=[T[0][0],T[1],(T[2]*T[1].inverse()*M*T[2])[0]]
        return [t*inits*self.initial_period_matrix for t in newT]
        
# find an intermediate point not far away from the end point and on the path.
# give it to the numerical_tms [a1,...,a{n-2},[a{n-1}],[a{n}],[a{n-1}]]



