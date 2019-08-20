

from sage.geometry.voronoi_diagram import VoronoiDiagram

######
# Takes a list of polynomials, and constructs a rectangular path from 0 to 1 avoiding
# the (non 0 or 1) roots of the polynomials
def voronoi_path(polylist):
    pts=set()
    P.<s>=PolynomialRing(QQ)
    print "Solving for dangerous points..."
    for p in polylist:
        p=P(p).radical()
        # numerical instability messes up the roots at 0 and 1, we will add them manually
        if p(0) == 0:
            p=P(p/s)
        if p(1) == 0:
            p=P(p/(s-1))
        roots=[r[0] for r in p.roots(ring=ComplexField(300))]
        pts = pts.union(set((real(r).simplest_rational(),imag(r).simplest_rational()) for r in roots))
    
    R=QQ
    pts.add((R(0),R(0)))
    pts.add((R(1),R(0)))
    if len(pts) == 2:
        return [0,1]
    # to discourage paths going in off directions
    [[pts.add((R(i),R(j))) for i in {-1,1}] for j in {-1,1}]

    print "Constructing Voronoi diagram..."
    vd = VoronoiDiagram(pts)

    print "Finding shortest path..."
    gr=Graph()
    for pt_, reg in vd.regions().iteritems():
        pt = tuple(pt_)
        pt = pt[0] + I*pt[1]
        ptc = CDF(pt)
        for edge in reg.bounded_edges():
            u = edge[0][0]+I*edge[0][1]
            v = edge[1][0]+I*edge[1][1]
            uc, vc = CDF(u), CDF(v)
            ratio = (uc-vc).abs()/min((ptc-uc).abs(), (ptc-vc).abs(), (ptc-(uc+vc)/2).abs())
            gr.add_edge(u, v, label=float(ratio))

        if pt == 0 or pt == 1:
            for v in reg.vertices():
                u = v[0]+I*v[1]
                gr.add_edge(pt, u, 1.0)
    
    return gr.shortest_path(QQ(0), QQ(1), by_weight=True), pts
