from sage.features.latte import Latte_count, Latte_integrate
from matplotlib.pyplot import figure

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
def meander(L1,L2,flow=True):
    n = sum(L1)
    fmeander=[]
    tb=blocks(L1)
    bb=blocks(L2)

    for t in tb:
        fmeander += list(zip(t, reversed(t), ['T']*len(t)))[:len(t)//2] 
    for t in bb:
        fmeander += list(zip(t, reversed(t), ['B']*len(t)))[:len(t)//2]
    
    if flow:
        fmeander += [(i+1, i+2, 'S') for i in range(n-1)]
        return DiGraph(fmeander, multiedges=True)
    else:
        return Graph(fmeander, multiedges=True)

"""
Takes the meander graph corresponding to the seaweed algebra 
defined by the compositions L1 and L2 written as lists, and
constructs a sink flow meander by directing every edge low-to-high
and connecting every vertex to a new sink vertex n+1.
"""
def sink_flow_meander(M):
    edges = list(M.edges())
    edges += [(i+1, M.order()+1, 'G') for i in range(M.order())] # G stands for gutter, since S is taken.
    return DiGraph(edges, multiedges=True)

"""
Generates the polytope corresponding to the seaweed algebra 
defined by the compositions L1 and L2 written as lists.
"""
def flow_poly(L1,L2):
    G=meander(L1,L2)
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
    plot.savefig(" ".join(["flowm", str(L1), str(L2)]) + ".png", dpi=480)
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
    plot.savefig(" ".join(["sflowm", str(L1), str(L2)]) + ".png", dpi=480)
    plot.clear()

# B1 = [2,4]
# B2 = [1,2,3]
B1 = [6,4,2]
B2 = [2,6,4]
M = meander(B1, B2)
S = sink_flow_meander(M)
draw_flow_meander(M, B1, B2)
draw_sink_flow_meander(S, B1, B2)
# print_meander(G, B1, B2)
# P = flow_poly(B1,B2)
# print(PVolume(P))
# po=P.ehrhart_polynomial()
# print(po)
# print(po.coefficient(4))


