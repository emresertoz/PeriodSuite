from sage.all import *
from ore_algebra import *
import logging
from monodromy import LocalMonodromy
logger = logging.getLogger(__name__)

class LMHS:
    def __init__(self,ivps,initial_periods):
        self.ivps = sorted(ivps, key=lambda ivp : ivp.label)
        self.odes = [ivp.ode for ivp in ivps]
        # self.initial_period_matrix = initial_period_matrix
        self.path = ivps[0].path # the paths of integration should all be the same
        self.initial_periods = initial_periods

        self.field = ivps[0].field
        self.refine_path()
        self.mark_path_for_loop()

    def refine_path(self):
        #FIXME make the penultimate point closer to the final point
        # But I won't do it now, so this function does nothing
        print("We should refine the path in the future. For now, see if it is OK.")
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
        T=ode.numerical_transition_matrix(self.marked_path,10**(-ivp.precision),assume_analytic=True)
        T=[t[1].change_ring(self.field) for t in T]
        M=LocalMonodromy(ode,end_point,self.field).monodromy
        inits=ivp.inits.change_ring(self.field)
# take the output matrices: T0, T1, T2 and the local monodromy M
# hit each of these below with *inits* (later we will multiply with the initial period matrix)
    # T1[0], T2, ((T3*T2.inverse())*M*T2)[0]
    # this gives you perids at penultimate, ultimate, and penultimate after monodromy
    # the one at ultimate needs all of the coordinates, so we keep all the matrix
        newT=[T[0][0],T[1],(T[2]*T[1].inverse()*M*T[1]*T[0].inverse())[0]]
        return [t*inits*self.initial_periods for t in newT]

    def get_compatible_expansions(self):
        """For the list of odes, increase the order of expansion at target point so that the maximal valuations match, and return these expansions."""
        odes = self.odes
        point = self.path[-1]
        vals = {}; orders = {}; new_orders = {};
        for ode in odes:
            vals[ode] = get_max_valuation_of_expansion(ode,point)
            orders[ode] = get_order_of_expansion(ode,point)
        maxVal = max(vals[ode] for ode in odes)
        for ode in odes:
            # adding "1" because I want to kill terms of maxVal+1 and I want some wiggle room
            new_orders[ode] = orders[ode]+maxVal-vals[ode]+1
        expansions = [ode.local_basis_expansions(point, order=new_orders[ode]) for ode in odes]
        # print([vals[ode] for ode in odes])
        # print([orders[ode] for ode in odes])
        # print([new_orders[ode] for ode in odes])
        # print([ceil(max([s[1].n for sol in sols for s in sol if s[0] != 0])) for sols in expansions])
        return expansions, maxVal

    def expansions_as_log_puiseux_series(self, coef_ring = None):
        """ Compute all expansions of odes at target point to sufficiently high valuation and then cast them to CC{{t}}[log(t)].  """
        expansions, maxVal = self.get_compatible_expansions()
        if coef_ring == None:
            coef_ring = self.field
        # FIXME!!
        # R = PuiseuxSeriesRing(coef_ring,'t',default_prec = maxVal+1)
        R = LaurentSeriesRing(coef_ring,'t',default_prec = maxVal+1)
        S = R['logt']
        t=R.gen(); L = S.gen();
        def convert(sol):
            # FIXME!!!!!!!!!!
            # return S(sum(s[0]*t**(s[1].n)*L**(s[1].k) for s in sol))
            return S(sum(s[0]*t**(2*s[1].n)*L**(s[1].k) for s in sol))
        return [[convert(sol) for sol in sols] for sols in expansions]

    def expansions_as_laurent_series(self,monodromy_mult, coef_ring = None):
        """ Set log terms to zero and substitute t -> t^monodromy_mult to kill denominators of fractional powers of t. Returns a Laurent series."""
        expansions, maxVal = self.get_compatible_expansions()
        if coef_ring == None:
            coef_ring = self.field
        R = LaurentSeriesRing(coef_ring,'t',default_prec = maxVal+1)
        t=R.gen();
        def convert(sol):
            return R(sum(s[0]*t**(monodromy_mult*s[1].n) for s in sol if s[1].k == 0))
        return [[convert(sol) for sol in sols] for sols in expansions]

#########################################
## Functions outside of the LMHS class ##
#########################################

def hirzebruch_series(d, prec = 10):
    r"""
    Return the generating series of the (primitive) Hodge numbers $h^{p,q}(d)-\delta_{p,q}$ of the middle cohomology of a smooth hypersurface of degree d, i.e., 
     $\sum_{p,q \ge 0} (h^{p,q}(d)-\delta_{p,q}) x^py^q$.
    """
    R = PowerSeriesRing(QQ,['x','y'], default_prec = prec)
    x,y = R.gens()
    # In Arapura Alg Geo over C, Thm 17.3.4, this is attributed to Hirzebruch
    return ((1+y)**(d-1) - (1+x)**(d-1))/((1+x)**d*y-(1+y)**d*x)

def hodge_numbers(n,d):
    r"""Return the primitive Hodge number $h^{p,n-p} - \delta_{p,n-p}$ of the middle cohomology of the projective hypersurface of dimension n and degree d."""
    H = hirzebruch_series(d, prec = n+10)
    cofs = H.coefficients()
    R = parent(H); x,y = R.gens()
    hp = []
    for p in range(0,n+1):
        q = n-p
        if x**p*y**q in cofs.keys():
            hp.append(cofs[x**p*y**q])
        else:
            hp.append(0)
    return hp

def hodge_filtration_dimensions(n,d):
    """See hodge_numbers."""
    hp = hodge_numbers(n,d)
    fp = []
    for ii in range(0,len(hp)):
        fp.append(sum(hp[:ii+1]))
    return fp


def get_max_valuation_of_expansion(ode,point,order=None):
    """ Compute the smallest integer bounding the valuation of the coefficients of the series expansions. 

    This computes the series expansion so is not very efficient.
    """
    sols = ode.local_basis_expansions(point, order = order)
    return ceil(max([s[1].n for sol in sols for s in sol if s[0] != 0]))

from ore_algebra.analytic.differential_operator import DifferentialOperator
from ore_algebra.analytic.path import Point
def get_order_of_expansion(ode,point):
    """
    Get the maximum order used by ore_algebra to expand local basis series of this ode by default. This is not necessarily the maximum valuation of terms.
    """
    # This initial part of this code is taken from ore_algebra
    # to compute their order in the same way they do
    mypoint = Point(point, ode)
    dop = DifferentialOperator(ode)
    ldop = dop.shift(mypoint)
    ind = ldop.indicial_polynomial(ldop.base_ring().gen())
    order = max(dop.order(), ind.dispersion()) + 3
    return order

def logarithm_of_monodromy(T):
    """
    For T a quasi-unipotent matrix return log(T), u, m where (T^m-1)^u = 0 and u,m are minimal.
    Note: Mondromy of polarized pure Hodge structures are quasi-unipotent (Borel).
    """
    poly=T.minimal_polynomial()
    if not poly.is_cyclotomic_product():
        raise ValueError("T is not quasi-unipotent.")
    mult=lcm([fact[0].is_cyclotomic(certificate=True) for fact in list(poly.factor())])
    varx=poly.parent().gen()
    unipotency=1
    while not poly.divides((varx**mult-1)**unipotency):
        unipotency +=1
    logT=1/mult*sum([(-1)**(k+1)/k*(T**mult-1)**k for k in range(1,unipotency)])
    return T.parent()(logT), unipotency, mult
        
def weight_filtration_of_nilpotent_matrix(N,k):
    """
    For N nilpotent matrix, k any positive integer such that N^(k+1) = 0, compute the weight filtration of Deligne and Schmid. 

    The changing k only shifts the indexing of the filtration.
    """
    if not (k in PositiveIntegers()):
        raise ValueError("k must be a positive integer")
    if not N**(k+1) == 0:
        raise ValueError("N must be nilpotent with N^(k+1)=0")

    V = VectorSpace(QQ,N.ncols()) 
    LN = linear_transformation(N)

    # the order in which we will compute the weight filtration
    # we zigzag diagonally up and down on this table:
    # W{-1} W0      W1      W2      ... W{k-1}
    # W{2k} W{2k-1} W{2k-2} W{2k-3} ... W{k}
    indices = flatten([[n,m] for n,m in zip(range(-1,k-1),range(2*k,k,-1))])
    # initial conditions
    W = {-1:V.zero_subspace(), (2*k):V}
    for a in indices:
        if a < k-1:
            W[2*k-a-2] = (LN**(k-a-1)).inverse_image(W[a])
        elif a > k:
            W[2*k-a] = (LN**(a-k))(W[a])
        
    basis = []
    for a in range(-1,2*k):
        if W[a].dimension() < W[a+1].dimension():
            Q = W[a+1]/W[a]
            basis.append([Q.lift(q) for q in Q.basis()])
            
    return matrix(flatten(basis)), [W[a].dimension() for a in range(-1,2*k+1)]


