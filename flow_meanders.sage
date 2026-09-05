from sage.features.latte import Latte_count, Latte_integrate

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

# Generates blocks of values over which have arcs in the meander
def blocks(L):
    psums = [sum(L[:i]) for i in range(len(L)+1)]
    blocks = [[i+1 for i in range(psums[k], psums[k+1])] for k in range(len(L))]
    return blocks

# Generates the flow meander digraph corresponding to the 
# seaweed algebra defined by the compositions L1 and L2 written as lists.
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

def gutter_flow_meander(M):
    edges = list(M.edges())
    edges += [(i+1, M.order()+1, 'G') for i in range(M.order())]
    return DiGraph(edges, multiedges=True)

def sink_flow_meander(M):
    edges = list(M.edges())
    edges.append((M.order(), M.order()+1, 'S'))
    return DiGraph(edges, multiedges=True)

# Generates the polytope corresponding to the seaweed algebra 
# defined by the compositions L1 and L2 written as lists
def flow_poly(L1,L2):
    G=meander(L1,L2)
    return G.flow_polytope()

# Calculates volume of polytope
def PVolume(P):
    poly=P.ehrhart_polynomial()
    d=poly.degree()
    v=poly.coefficient(d)

    return v*factorial(P.dimension())
    
def print_meander(M, L1, L2):
    n = M.order()
    plot = M.plot(layout="circular", color_by_label={'T':"blue", 'B':"green", 'S':"red", 'G':"orange"})
    plot.save_image("./" + " ".join(["meander", str(L1),  str(L2)]) + ".png")

def draw_flow_meander(M, L1, L2):
    n = M.order()
    fm = list(M.edges())
    plots = []
    vsize = 0.1
    curves = {'S':0, 'T':0.5, 'B':-0.5}
    heads = {'S':(-vsize,0), 'T':(0, vsize), 'B':(0, -vsize)}

    for e in fm:
        p = [[(e[0],0), ((e[1]+e[0])/2,(e[1]-e[0])*curves[e[2]]), (e[1]+heads[e[2]][0], 0+heads[e[2]][1])]]
        plots.append(sage.plot.plot.plot(arrow2d(path=p, color=(0,0,0), axes=False)))

    for v in range(n):
        plots.append(sage.plot.plot.plot(circle((v+1, 0), vsize, fill=True, rgbcolor=(0,0,0), axes=False)))
        plots.append(text(str(v+1), (v+1, 0), fontsize="large", rgbcolor=(1,1,1), zorder=5, axes=False))

    plot = sum(plots)
    plot.save_image(" ".join(["flowm", str(L1), str(L2)]) + ".png")

# dg = DiGraph([(1,2), (1,3), (2,3)])
# w = (1,1)
# P = weighted_flow_polytope(dg, w)
# print(P.vertices(), P.dimension(), P.volume())

B1=[2,4]
B2=[1,2,3]
M = meander(B1, B2)
S = sink_flow_meander(M)
G = gutter_flow_meander(M)
draw_flow_meander(M, B1, B2)
draw_flow_meander(S, B1, B2)
print_meander(G, B1, B2)
# P = flow_poly(B1,B2)
# print(PVolume(P))
# po=P.ehrhart_polynomial()
# print(po)
# print(po.coefficient(4))


