from ore_algebra import *
from ore_algebra.analytic.path import Point
from sage.all import *

class LocalMonodromy:
    def __init__(self,ode,point,field):
        """
        Computes the local monodromy matrix on the reduced Frobenius basis 
        of solutions of an ode at a regular point. 
        The matrix acts on the left on column vectors.

        EXAMPLES::

        from ore_algebra import *
        Dops, x, D = DifferentialOperators()

        # try any of these
        # ode = x*D-1/2; point = 0
        # ode=(x*D+1)*D; point = 0
        # ode = (D^3 + ((24*x^2 - 4*x - 12)/(8*x^3 - 8*x))*D^2 + ((32*x^2 + 32*x - 16)/(32*x^4 + 32*x^3 - 32*x^2 - 32*x))*D); point=0
        ode = (x^2*(x^2-34*x+1)*D^3 + 3*x*(2*x^2-51*x+1)*D^2 + (7*x^2-112*x+1)*D + (x-5)); point = 0

        print(ode.local_basis_expansions(point))
        lm = LocalMonodromy(ode,point,ComplexBallField(100))
        pretty_print(lm.monodromy.apply_map(lambda x : x.mid()))
        """

        self.ode = ode
        self.point = point
        self.tpi = 2*field.pi()*field.gens()[0] # 2*Pi*I
        self.exponents = self.exponents_of_leading_monomials()
        self.sols = ode.local_basis_expansions(point)
        self.monodromy = self.monodromy_matrix()

    def exponents_of_leading_monomials(self):
        pt = Point(self.point, self.ode)
        if not (pt.is_ordinary() or pt.is_regular()):
            raise ValueError("Not a regular singularity.")
        struct = pt.local_basis_structure()
        return [(s.valuation,s.log_power) for s in struct]            

    def monomial_contribution(self,exponent1,exponent2):
        # Watch out for the Frobenius basis conventions, no b! but there is a d!
        a,b = exponent1; # monomial1 = x^a*log(x)^b
        c,d = exponent2; # monomial2 = 1/d! * x^c*log(x)^d
        # Compute the contribution of monomial1 under monodromy to monomial2.
        if a == c and d <= b:
            return factorial(d)*exp(self.tpi*a)*binomial(b,d)*self.tpi**(b-d)
        else:
            return 0
    
    def total_monomial_contribution(self,sol,exponent):
        # sol comes from local_basis_expansions
        # we find total contribution to the monomial corresponding to given exponent
        return sum(s[0]*self.monomial_contribution((s[1].n,s[1].k),exponent) for s in sol if s[0] != 0)

    def monodromy_action(self,sol):
        # sol comes from local_basis_expansions
        # we find the new coordinates of sol under monodromy
        return [self.total_monomial_contribution(sol,expo) for expo in self.exponents]

    def monodromy_matrix(self):
        # Compute the local monodromy matrix of going counterclockwise around point.
        # The matrix acts on the left on column vectors.
        return Matrix([self.monodromy_action(sol) for sol in self.sols]).transpose()
