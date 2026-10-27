### Most of this code file relies on techniques elaborated on in 
### "Meander Graphs and Frobenius Seaweed Lie Algebras II" by
### Coll, Hyatt, Magnant and Wang, referred to as CHMW here.

def flip_down(c, d):
    if c[0] < d[0]:
        return (d,c)
    else:
        print("Error occured in a flip move on compositions", c, "and", d)
        return ("ERROR", "ERROR")

def component_elim_down(c, d):
    if c[0] == d[0]:
        elimcomp = c[0]
        c = c[1:]
        d = d[1:]
        return (c,d,elimcomp)
    else:
        print("Error occured in a component elimination move on compositions", c, "and", d)
        return ("ERROR", "ERROR")

def block_elim_down(c, d):
    if c[0] == d[0]*2:
        c[0] = d[0]
        d = d[1:]
        return (c,d)
    else:
        print("Error occured in a block elimination move on compositions", c, "and", d)
        return ("ERROR", "ERROR")

def rot_contraction_down(c, d):
    if d[0] < c[0] and c[0] < d[0]*2:
        temp = c[0]
        c[0] = d[0]
        d[0] = 2*d[0]-temp
        return (c,d)
    else:
        print("Error occured in a rotational contraction move on compositions", c, "and", d)
        return ("ERROR", "ERROR")

def pure_contraction_down(c, d):
    if c[0] > d[0]*2:
        c = [c[0] - (2*d[0])] + c
        c[1] = d[0]
        d = d[1:]
        return (c,d)
    else:
        print("Error occured in a pure contraction move on compositions", c, "and", d)
        return ("ERROR", "ERROR")

def flip_up(c, d):
    return (d,c)

def component_elim_up(c, d, e):
    return ([e] + c, [e] + d)

def block_elim_up(c, d):
    d = [c[0]] + d
    c[0] *= 2
    return (c,d)

def rot_contraction_up(c, d):
    temp = c[0]
    c[0] = 2*temp-d[0]
    d[0] = temp
    return (c,d)

def pure_contraction_up(c, d):
    temp = c[1]
    c[1] = c[0] + 2*temp
    c = c[1:]
    d = [temp] + d
    return (c,d)

"""
Given a pair of compositions and a set of winding down moves,
computes the resulting (not necessarily empty!) compositions.
We call the sequence of moves the signature here, even though
it might not be the whole signature, or not even a part of it.
"""
def winding_down(signature, c, d, elimcomp=False):
    moves = list(signature.upper())
    elim_comps = []
    for m in moves:
        if c == "ERROR":
            return
        
        # Flip of compositions
        elif m == 'F':
            c,d = flip_down(c,d)

        # Component elimination
        elif m == 'C':
            c,d,e = component_elim_down(c,d)
            elim_comps.append(e)

        # Block elimination
        elif m == 'B':
            c,d = block_elim_down(c,d)

        # Rotation contraction
        elif m == 'R':
            c,d = rot_contraction_down(c,d)

        # Pure contraction
        elif m == 'P':
            c,d = pure_contraction_down(c,d)

        # None of the above - raise error
        else:
            print("Error occured in reading the signature: ")
            print("\tat move", m)
            print("\twith top composition", c)
            print("\tand bottom composition", d)
            return "ERROR"
    
    if elimcomp:
        return (c,d,e)
    return (c,d)

"""
Given a pair of (possibly empty!) compositions and a set of 
winding up moves, computes the resulting compositions.
We call the sequence of moves the signature here, even though
it might not be the whole signature. Note that the given
signature in this method will have numbers for every 
component elimination move, meaning that it is not compatible
with the winding_down method. [CODE SEQUENCE TRANSLATOR METHOD]
"""
def winding_up(signature, c=[], d=[]):
    signature = signature.upper()
    moves = []
    valid_moves = ['F', 'C', 'B', 'R', 'P']
    elim_comps = []
    
    # Preprocessing the signature to get component sizes
    number = False
    digits = []
    for m in signature:
        if m in valid_moves:
            number = False
            moves.append(m)
        elif m.isnumeric():
            number = True
        else:
            return "ERROR"
        
        if number:
            digits.append(m)
        elif not number and len(digits) > 0:
            elim_comps.append(int("".join(digits)))
            digits = []
    if number:
        elim_comps.append(int("".join(digits)))

    # The input is a signature of some meander, so reconstructing 
    # it from empty requires reading it right-to-left.
    moves.reverse()
    elim_comps.reverse()

    for m in moves:
        # Flip of compositions
        if m == 'F':
            c,d = flip_up(c,d)

        # Component elimination
        elif m == 'C':
            c,d = component_elim_up(c,d,elim_comps[0])
            elim_comps = elim_comps[1:]

        # Block elimination
        elif m == 'B':
            c,d = block_elim_up(c,d)

        # Rotation contraction
        elif m == 'R':
            c,d = rot_contraction_up(c,d)

        # Pure contraction
        elif m == 'P':
            c,d = pure_contraction_up(c,d)

    return (c,d)

"""
Given two compositions, returns the signature of the 
associated meander. Note that this method is similar
to the winding_down method - the only difference is
that we do not supply a signature, but return it. 
It might be better to split every move into its own
method, to be reused by both methods. 
"""
def signature(c, d, elimcomp=False):
    moves = []
    elim_comps = []
    while c:
        if c[0] > d[0]:
            c,d = flip_down(c,d)
            moves.append('F')

        elif c[0] == d[0]:
            c,d,e = component_elim_down(c,d)
            elim_comps.append(e)
            moves.append('C')

        elif c[0] == d[0]*2:
            c,d = block_elim_down(c,d)
            moves.append('B')

        elif d[0] < c[0] and c[0] < d[0]*2:
            c,d = rot_contraction_down(c,d)
            moves.append('R')

        elif c[0] > d[0]*2:
            c,d = pure_contraction_down(c,d)
            moves.append('P')

    if elimcomp:
        return "".join(moves), elim_comps
    return "".join(moves)

"""
Detects if a meander is Frobenius by its signature.
This algorithm comes from Theorem 5 in CHMW
"""
def is_frobenius(c,d):
    moves, elim_comps = signature(c,d)
    return (moves.count('C') == 1 and moves[-1] == 'C' and elim_comps[0] == 1)

def index(c,d):
    moves, elim_comps = signature(c,d)
    return sum(elim_comps)

sgn = "PFRBFBFBC1"
print(winding_up(sgn))

