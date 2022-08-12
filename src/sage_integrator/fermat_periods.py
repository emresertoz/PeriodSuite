# Class to compute periods of Fermat type hypersurfaces in projective space
import time

class FermatPeriods:
    def __init__(self,degree,coefficients,digit_precision):
        """degree (d) and coefficients (c0,...,c{n+1}) determine the Fermat polynomial
        c0*x0^d+c1*x1^d+...+cn*x{n+1}^d. This class computes the period matrix for the 
        resulting hypersurface in P^{n+1} with precision 10^-digit_precision. The details are
        worked out in arXiv:1803.08068."""

        self.d = ZZ(degree)
        self.n = ZZ(len(coefficients)-2)
        self.coefs = coefficients
        self.R = PolynomialRing(QQ, self.n+2, 'x', order='degrevlex')
        self.xi = exp(2*pi*I/d) #d-th root of unity
        self.inverse_roots = [c**(-1/d) for c in coefficients[:-1]] + [(-c)**(-1/d) for c in coefficients[-1:]]

        self.bit_precision = ceil(digit_precision*log(10)/log(2))+50
        self.field = ComplexField(self.bit_precision)

        # maybe you don't have to do this now!
        print("Computing the exact period matrix of Fermat.")
        t0=time.time()
        self.compute_cohomology()
        self.compute_homology()
        self.compute_intersection_matrix()
        self.compute_period_matrix()
        print("Done in %f seconds."%(time.time()-t0))

        print("Computing the numerical approximation of the period matrix of Fermat.")
        t0=time.time()
        self.approximate_period_matrix = self.period_matrix.change_ring(self.field)
        print("Done in %f seconds."%(time.time()-t0))
        #DEBUG
        # print(self.cohomology_basis)
        # print(self.homology_basis)
        # print(self.intersection_matrix)
        # print(self.period_matrix)

    def compute_cohomology(self):
        """Compute the grevlex Griffiths basis for the cohomology of this Fermat hypersurface. We store the exponent vectors of the monomials defining the basis.
        """
        n = self.n; d = self.d; R = self.R;
        # x^a*Omega/f^l --- given the pole order l and all but last entry of the exponent a, we can determine the last entry of a. This is because of the rectangular nature of the Jacobian ideal of Fermat.
        is_admissable = lambda exp,l : l*d-(n+2)-sum(exp) in range(d-1); 
        exps = cartesian_product([range(d-1) for _ in range(n+1)])
        exps = [ list(exp)+[l*d-(n+2)-sum(exp)] for exp in exps for l in range(1,n+2) if is_admissable(exp,l)]
        mons = [self.R.monomial(*exp) for exp in exps]
        self.cohomology_basis=[list(m.exponents()[0]) for m in sorted(mons)]

    def compute_homology(self):
        """ 
        Compute the Pham basis for the primitive homology of this Fermat hypersurface. 
        """
        # The Pham basis is obtained by translating one standard cycle under the group action that is coordinatewise scaling by roots of unity. The group elements are viewed as monomials (a power of root of unity acting on a coordinate) and their exponent vectors are stored. See arXiv:1803.08068 for more details.
        ## 
        # we drop one variable now, but account for it later
        N=self.n+1; d=self.d; 
        P=PolynomialRing(QQ,N,'y',order="degrevlex")
        V=list(P.gens())+[prod(P.gens())]
        I=Ideal([sum(v**r for r in range(d)) for v in V])
        mons=sorted(I.normal_basis())
        # we are adding a 0 to the exponent vectors, i.e., the extra variable has exponent 0
        self.homology_basis=[list(m.exponents()[0])+[0] for m in mons]

    # TODO: This function is called a lot, and per call the eps function is called a lot, each creating the little dictionary. Both eps and the permutation g can be created outside of cap_product for efficiency.
    def cap_product(self,hom_cycle_1,hom_cycle_2):
        """ Compute the intersection product of two homology cycles in a Pham basis. """
        # See Theorem 1.3 of Shimada
        # note that beta's always have 0 at the last entry
        # and this entry should be dropped when using these formulas
        d = self.d; N = self.n+2;
        # equivariance under group action
        # we translate one of them back to the standard cycle
        beta = [mod(a-b,d) for a,b in zip(hom_cycle_1,hom_cycle_2)]
        # to use Shimada's results, we need to change back to their conventions
        # this is a simple permutation of coordinates
        g = Permutation([N]+list(range(1,N)))
        beta = permutation_action(g,beta)
        # now get the last coordinate to be 0 again and drop it from notation
        beta = [ b - beta[-1] for b in beta[0:-1] ]
        # auxillary values
        eps = lambda x : {0:ZZ(1),1:ZZ(-1)}.get(x,ZZ(0))
        term1 = prod(eps(b) for b in beta)
        term2 = prod(eps(b+1) for b in beta) # b+1 is computed mod d
        sign = (-1)**(((N-2)*(N-1)/2))
        return sign*(term1-term2)

    def compute_intersection_matrix(self):
        """Compute the intersection product of the Pham basis for the primitive homology of this Fermat hypersurface."""
        B = self.homology_basis
        cap = self.cap_product
        self.intersection_matrix = matrix([[cap(b1,b2) for b1 in B] for b2 in B]).change_ring(ZZ)

    def integration_pairing(self,cohom_class,hom_class):
        """ Compute the exact integration pairing of the given cohomology class over the homology class. """
        d = self.d; xi = self.xi; scalars = self.inverse_roots
        beta=hom_class;
         # increment exponents to simplify formulas (and for compatibility with literature)
        alpha=[ZZ(c+1) for c in cohom_class];
        l=ZZ(sum(alpha)/d) ## pole order of cohom class
        scalar0=prod([(1-alpha[-1]/(j*d)) for j in range(1,l)])
        scalar1=prod([scalars[i]**alpha[i] for i in range(0,len(alpha))])
        scalar2=prod([(1-xi**-a)/d for a in alpha[0:-1]])
        scalar3=xi**sum([alpha[i]*beta[i] for i in range(0,len(alpha))])
        scalar4=prod(gamma(a/d) for a in alpha[0:-1])/gamma(sum([a/d for a in alpha[0:-1]]))
        # print(scalar0,scalar1,scalar2,scalar3,scalar4)
        return scalar0*scalar1*scalar2*scalar3*scalar4

    def compute_period_matrix(self):
        self.period_matrix = Matrix([[self.integration_pairing(a,b) for b in self.homology_basis] for a in self.cohomology_basis])

## test

# fp=FermatPeriods(4,[1,1,1,1],100)
# d=4; n=2; N=4;
# hom_cycle_1,hom_cycle_2=[2, 1, 0, 0], [1, 1, 0, 0]
# fp.cap_product(hom_cycle_1,hom_cycle_2)
# print(fp.bit_precision)

# print(fp.period_matrix)
