from sage.features.latte import Latte_count, Latte_integrate
from matplotlib.pyplot import figure
from sage.graphs.connectivity import connected_components_subgraphs
from numpy import concatenate as sumarr
from sage.matrix.special import ith_to_zero_rotation_matrix as rotate
import csv

"""
Computes the flow polytope \\mathcal{F}_G(a) for a 
directed graph G and a nonegative integer vector a.
This is an adaptation of Sage's flow_polytope method, 
which has a fixed at the origin. The source code is
from Sage, and writing credit goes to the developers.
"""
def weighted_flow_polytope(digraph, weight, edges=None, ends=None, backend=None):
    from sage.geometry.polyhedron.constructor import Polyhedron
    
    edges = list(digraph.edges())
    m = len(edges)
    ineqs = [[0] * (i + 1) + [1] + [0] * (m - i - 1) for i in range(m)]

    eqs = []
    for u in digraph.vertices():
        ins = set(digraph.incoming_edge_iterator(u))
        outs = set(digraph.outgoing_edge_iterator(u))
        eq = [Integer(j in ins) - Integer(j in outs) for j in edges]

        const = 0
        if not outs:  # sinks (outdegree 0)
            const = -sum(weight)
        else:
            const = weight[u-1]

        eq = [const] + eq
        eqs.append(eq)

    return Polyhedron(ieqs=ineqs, eqns=eqs, backend=backend)

def reduce_dim(pol):
    pass

"""
Performs 2D rotations on a given polytope to reduce 
its ambient space to the dimension of the polytope itself.
"""
def ambient_reduction(pol, nonneg=False, basis=False, symb=False):
    amb = pol.ambient_dim()
    dim = pol.dim()

    # Selecting a cornered simplex
    chain = list(pol.face_lattice().maximal_chains()[0])
    simp_verts = []
    for i in range(1, len(chain)):
        for vert in chain[i].vertices():
            if vert not in chain[i-1].vertices():
                simp_verts.append(vert)
                break
    
    # Translating polytope to have a corner at 0
    t = [-coor for coor in list(simp_verts[0])]
    pol = pol.translation(t)

    # Turn every vertex object to a vector
    simp_verts = [vector(v)+vector(t) for v in simp_verts]
    all_verts = [vector(v) for v in list(pol.vertices())]
    base_verts = False
    if basis:
        base_verts = [vector([0]*(i) + [1] + [0]*(dim-1-i)) for i in range(amb)]

    # Rotations
    for i in range(dim):
        for j in range(amb-1, i, -1):
            vert = simp_verts[i+1]
            rot = rotate(vert, j)
            simp_verts = [rot*v for v in simp_verts]
            all_verts = [rot*v for v in all_verts]
            if basis:
                base_verts = [rot*v for v in base_verts]

    # Reductions
    all_verts = [v[:dim] for v in all_verts]

    # Create new polytope:
    numeric_verts = [[numerical_approx(coor) for coor in v] for v in all_verts]
    new_pol = Polyhedron(vertices=numeric_verts)

    # Flags
    if nonneg:
        t = [0]*dim
        for vert in all_verts:
            for j in range(dim):
                if vert[j] < t[j]:
                    t[j] = vert[j]
        t = [-coor for coor in t]
        new_pol = new_pol.translation(t)

    ret = [new_pol]
    if basis:
        ret = ret + [base_verts]
    if symb:
        ret = ret + [all_verts]
    return ret


"""
Generates blocks of values over which have arcs in the meander.
"""
def blocks(L):
    psums = [sum(L[:i]) for i in range(len(L)+1)]
    blocks = [[i+1 for i in range(psums[k], psums[k+1])] for k in range(len(L))]
    return blocks

"""
Generates the flow meander digraph corresponding to the 
seaweed algebra defined by the compositions L1 and L2 written as lists.
flow=False returns an undirected meander as described in DK2000.
"""
def meander(L1, L2):
    n = sum(L1)

    fmeander=[]
    tb=blocks(L1)
    bb=blocks(L2)

    for t in tb:
        fmeander += list(zip(t, reversed(t), ['T']*len(t)))[:len(t)//2] 
    for t in bb:
        fmeander += list(zip(t, reversed(t), ['B']*len(t)))[:len(t)//2]
    return Graph([list(range(1, n+1)), fmeander], multiedges=True, format="vertices_and_edges")

"""
Takes the meander graph corresponding to the seaweed algebra 
defined by the compositions L1 and L2 written as lists, and
constructs a flow meander by adding a path through the vertices, 
starting at 1 and ending at n.
"""
def flow_meander(L1, L2):
    M = meander(L1, L2)
    edges = list(M.edges())
    edges += [(i+1, i+2, 'S') for i in range(M.order()-1)]
    return DiGraph(edges, multiedges=True)

"""
Takes the meander graph corresponding to the seaweed algebra 
defined by the compositions L1 and L2 written as lists, and
constructs a sink flow meander by connecting every vertex 
to a new sink vertex n+1.
"""
def sink_flow_meander(L1, L2):
    M = meander(L1, L2)
    edges = list(M.edges())
    edges += [(i+1, M.order()+1, 'G') for i in range(M.order())] # G stands for gutter, since S is taken.
    return DiGraph(edges, multiedges=True)

"""
Returns all the directed paths of a flow meander as lists of labelled edges.
Labels only tell us the type of edge (top, bottom, spine), not the component.
draw=True gives a string of Unicode characters representing the path.
"""
def directed_paths(L1, L2, draw=False):
    fmdr = flow_meander(L1, L2)
    paths = fmdr.all_paths(1, fmdr.order(), use_multiedges=True, report_edges=True, labels=True)
    if not draw:
        return paths
    
    spine = chr(8212)
    top = chr(8254)
    bottom = '_'
    drawings = []
    for path in paths:
        drawing = ""
        for edge in path:
            if edge[2] == 'S':
                drawing += spine
            elif edge[2] == 'B':
                drawing += bottom*(edge[1]-edge[0])
            elif edge[2] == 'T':
                drawing += top*(edge[1]-edge[0])
        drawings.append(drawing)
    return drawings
    
    # Generate all of them, count will come and get proven later
    # num=False flag for just returning the number (we can put a formula here later if needed)
    # Use something from the graph library

"""
Generates the polytope corresponding to the seaweed algebra 
defined by the compositions L1 and L2 written as lists.
"""
def flow_poly(L1,L2):
    G=flow_meander(L1,L2)
    return G.flow_polytope()

"""
Calculates the normalized volume of a polytope.
"""
def PVolume(P):
    poly=P.ehrhart_polynomial()
    d=poly.degree()
    v=poly.coefficient(d)

    return v*factorial(P.dimension())
    
def print_meander(M, L1, L2):
    n = M.order()
    plot = M.plot(layout="circular", color_by_label={'T':"blue", 'B':"green", 'S':"red", 'G':"orange"})
    plot.save_image("./" + " ".join(["meander", str(L1),  str(L2)]) + ".png")

"""
Draws a meander graph using matplotlib plots.
(The usual graph plotting functions do not support 
curved edges, so the graphs they produce are not adequate.)
"""
def draw_meander(M, L1, L2):
    n = M.order()
    fm = list(M.edges())
    plots = []
    vsize = 0.2
    curves = {'T':0.5, 'B':-0.5}
    heads = {'T':(0, vsize), 'B':(0, -vsize)}

    for e in fm:
        p = [[(e[0],0), ((e[1]+e[0])/2,(e[1]-e[0])*curves[e[2]]), (e[1]+heads[e[2]][0], 0+heads[e[2]][1])]]
        plots.append(sage.plot.plot.plot(bezier_path(path=p, color=(0,0,0), axes=False)))

    for v in range(n):
        plots.append(sage.plot.plot.plot(circle((v+1, 0), vsize, fill=True, rgbcolor=(0,0,0), axes=False)))
        plots.append(text(str(v+1), (v+1, 0), fontsize="medium", rgbcolor=(1,1,1), zorder=5, axes=False))

    plot = sum(plots)
    plot = plot.matplotlib(axes=False, figsize=(n,n/4))
    plot.savefig(" ".join(["meander", str(L1), str(L2)]) + ".png", dpi=480)
    plot.clear()

"""
Draws a flow meander using matplotlib plots.
(The usual graph plotting functions do not support 
curved edges, so the graphs they produce are not adequate.)
"""
def draw_flow_meander(M, L1, L2):
    n = M.order()
    fm = list(M.edges())
    plots = []
    vsize = 0.2
    curves = {'S':0, 'T':0.5, 'B':-0.5}
    heads = {'S':(-vsize,0), 'T':(0, vsize), 'B':(0, -vsize)}

    for e in fm:
        p = [[(e[0],0), ((e[1]+e[0])/2,(e[1]-e[0])*curves[e[2]]), (e[1]+heads[e[2]][0], 0+heads[e[2]][1])]]
        plots.append(sage.plot.plot.plot(arrow2d(path=p, color=(0,0,0), axes=False)))

    for v in range(n):
        plots.append(sage.plot.plot.plot(circle((v+1, 0), vsize, fill=True, rgbcolor=(0,0,0), axes=False)))
        plots.append(text(str(v+1), (v+1, 0), fontsize="medium", rgbcolor=(1,1,1), zorder=5, axes=False))

    plot = sum(plots)
    plot = plot.matplotlib(axes=False, figsize=(n,n/4))
    plot.savefig(" ".join(["fmeander", str(L1), str(L2)]) + ".png", dpi=480)
    plot.clear()

"""
Draws a sink flow meander using matplotlib plots.
(The usual graph plotting functions do not support 
curved edges, so the graphs they produce are not adequate.)
"""
def draw_sink_flow_meander(M, L1, L2):
    n = M.order()-1
    fm = list(M.edges())
    plots = []
    vsize = 0.2
    curves = {'S':0, 'T':0.5, 'B':-0.5}
    heads = {'S':(-vsize,0), 'T':(0, vsize), 'B':(0, -vsize)}
    sdepth = -0.5*max(L2)

    for e in fm:
        # Draw arrows to sink
        if e[2] == 'G':
            p = [[(e[0],0), (e[0], sdepth)]]
            plots.append(sage.plot.plot.plot(arrow2d(path=p, color=(0,0,0), axes=False)))
        # Draw arrows between non-sink vertices
        else:
            p = [[(e[0],0), ((e[1]+e[0])/2,(e[1]-e[0])*curves[e[2]]), (e[1]+heads[e[2]][0], 0+heads[e[2]][1])]]
            plots.append(sage.plot.plot.plot(arrow2d(path=p, color=(0,0,0), axes=False)))

    # Draw vertices with labels
    for v in range(n):
        plots.append(sage.plot.plot.plot(circle((v+1, 0), vsize, fill=True, rgbcolor=(0,0,0), axes=False)))
        plots.append(text(str(v+1), (v+1, 0), fontsize="medium", rgbcolor=(1,1,1), zorder=5, axes=False))
    # Draw sink with label
    plots.append(sage.plot.plot.plot(polygon2d([(1-vsize, sdepth), (1-vsize, sdepth-2*vsize), (n+vsize, sdepth-2*vsize), (n+vsize, sdepth)], fill=True, rgbcolor=(0,0,0), axes=False)))
    plots.append(text(str(n+1), ((n+1)/2, sdepth-vsize), fontsize="medium", rgbcolor=(1,1,1), zorder=5, axes=False))

    plot = sum(plots)
    plot = plot.matplotlib(axes=False, figsize=(n,n/2))
    plot.savefig(" ".join(["sfmeander", str(L1), str(L2)]) + ".png", dpi=480)
    plot.clear()

"""
Computes the covolume of a flow meander polytope associated with two compositions. 
This method is based on a polytope basis conjecture that is not proven as of 23/10/2026.
"""
def covolume(L1, L2):
    verts = []
    mdr = meander(L1, L2)
    n = mdr.order()

    top_edges = [e for e in mdr.edges() if e[2] == 'T']
    top_edges.sort()
    top_offset = n-1
    top_starts = {top_edges[i][0]:i+top_offset for i in range(len(top_edges))}
    
    bot_edges = [e for e in mdr.edges() if e[2] == 'B']
    bot_edges.sort()
    bot_offset = n-1 + len(top_edges)
    bot_starts = {bot_edges[i][0]:i+bot_offset for i in range(len(bot_edges))}
    
    for edge in mdr.edges():
        v = [0] * (mdr.size()+n-1)
        if edge[2] == 'T':
            v[top_starts[edge[0]]] += 1
        if edge[2] == 'B':
            v[bot_starts[edge[0]]] += 1
        for i in range(edge[0], edge[1]):
            v[i-1] -= 1
        verts.append(vector(v))
    
    V = matrix(verts)
    V_T = column_matrix(verts)
    G = V * V_T
    return(sqrt(det(G)))


def data_dictionary(L1, L2, vec):
    data = dict()
    M = meander(L1, L2)
    F = flow_meander(L1, L2)
    S = sink_flow_meander(L1, L2)
    data["top_comp"] = "|".join([str(c) for c in L1])
    data["bottom_comp"] = "|".join([str(c) for c in L2])
    pieces = connected_components_subgraphs(M)
    data["num_pieces"] = len(pieces)
    data["cycles"] = 0
    data["paths"] = 0
    data["points"] = 0
    for pc in pieces:
        if pc.order() == 1:
            data["points"] += 1
        elif pc.is_cycle():
            data["cycles"] += 1
        else:
            data["paths"] += 1
    data["index"] = 2*data["cycles"] + data["paths"] + data["points"]

    FP = weighted_flow_polytope(F, vec[:-1])
    fpoly = FP.ehrhart_polynomial()
    data["fm_dimension"] = fpoly.degree()
    data["fm_volume"] = fpoly.lc()
    data["fm_nvolume"] = data["fm_volume"]*factorial(data["fm_dimension"])
    data["fm_latpts"] = fpoly(1)
    data["fm_ehrhart"] = str(fpoly)
    data["fm_fvector"] = FP.f_vector()
    
    SFP = weighted_flow_polytope(S, vec)
    sfpoly = SFP.ehrhart_polynomial()
    data["sfm_dimension"] = sfpoly.degree()
    data["sfm_volume"] = sfpoly.lc()
    data["sfm_nvolume"] = data["sfm_volume"]*factorial(data["sfm_dimension"])
    data["sfm_latpts"] = sfpoly(1)
    data["sfm_ehrhart"] = str(sfpoly)
    data["sfm_fvector"] = SFP.f_vector()

    return data

params = ["No.", "top_comp", "bottom_comp", "num_pieces", "cycles", "paths", "points", "index",
        "fm_dimension", "fm_volume", "fm_nvolume", "fm_ehrhart", "fm_fvector", 
        "sfm_dimension", "sfm_volume", "sfm_nvolume", "sfm_ehrhart", "sfm_fvector"]

curr_params = ["sfm_dimension", "sfm_volume", "sfm_nvolume", "sfm_ehrhart", "sfm_fvector"]

bound = 11

for n in range(5, bound):
    comp1 = [n-3,1,1,1]
    comp2 = [n]
    print(n, covolume(comp1, comp2))

# bound = 8
# bar = 1
# for n in range(2, bound+1):
#     comps = list(Compositions(n))
#     vec = tuple([1] + [0]*(n-2))
#     for i in range(len(comps)):
#         for j in range(i):
#             comp1 = list(comps[i])
#             comp2 = list(comps[j])
#             mdr = meander(comp1, comp2)
#             edges = mdr.size()
#             fmdr = flow_meander(comp1, comp2)
#             fpol = weighted_flow_polytope(fmdr, vec)
#             dim = fpol.dim()
#             if edges != dim:
#                 print('|'.join([str[c] for c in comp1]), '|'.join([str[c] for c in comp2]), edges, dim, sep='\t')
#             if bar % 500 == 0:
#                 print(bar)
#             bar += 1