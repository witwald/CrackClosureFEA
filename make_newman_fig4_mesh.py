"""
Centre-cracked panel mesh following Newman (1974), Figure 4, for program FINAL.

Quarter model (symmetry about x = 0 and y = 0):
    half-width b = 9.0 in, height H = 10.0 in, initial crack half-length a = 1.075 in.
Crack-tip block (Figure 4, lower view): x = 0.70 .. 1.50 in, y = 0 .. 0.40 in,
    built on a 32h x 16h grid with h = 0.025 in.  Uniform h x h squares (split by
    "/" diagonals) cover x = 12h..22h, y = 0..3h, with node A (initial crack tip)
    at x = 15h = 1.075 in and nodes B..F at 16h..20h.  Elements coarsen by
    factors of two (h, 2h, 4h, 8h) through transition cells that use a centre
    node, as in the figure.
Outer region (Figure 4, upper view): graded triangles from 0.2 in at the block
    to about 3 in at the top edge, with boundary nodes placed where Newman placed
    them (measured from the figure).

Usage:  python make_newman_fig4_mesh.py [Smax] [ncycles] [out.dat]
Writes the input deck and a PNG picture of the mesh (same name, .png).
"""
import sys, math
import numpy as np
import triangle as tr
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

SMAX  = float(sys.argv[1]) if len(sys.argv) > 1 else 15.0
NCYC  = int(sys.argv[2])   if len(sys.argv) > 2 else 4
OUT   = sys.argv[3]        if len(sys.argv) > 3 else 'ccp_newman_fig4.dat'

B, H, A0 = 9.0, 10.0, 1.075
h  = 0.025
X0 = A0 - 15*h                    # 0.70: left edge of crack-tip block
NX, NY = 32, 16                   # block size in units of h

# ---------------------------------------------------------------- block (quadtree)
def want(xc, yc, s):
    """Desired cell size (units of h) at cell centre."""
    # distance (units of h) from the fine box [12,22] x [0,3]
    dx = max(12 - xc, 0, xc - 22); dy = max(yc - 3, 0)
    d = max(dx, dy)
    if d <= 0:  return 1
    if d <= 2:  return 2
    if d <= 6:  return 4
    return 8

cells = {(i, j, 8) for i in range(0, NX, 8) for j in range(0, NY, 8)}  # (x0,y0,size)
def split(c):
    x, y, s = c; t = s // 2
    return [(x, y, t), (x+t, y, t), (x, y+t, t), (x+t, y+t, t)]
changed = True
while changed:
    changed = False
    for c in list(cells):
        x, y, s = c
        if s > 1 and want(x + s/2, y + s/2, s) < s:
            cells.remove(c); cells.update(split(c)); changed = True
# 2:1 balance
def leaf_at(px, py):
    for c in cells:
        x, y, s = c
        if x <= px < x+s and y <= py < y+s: return c
    return None
changed = True
while changed:
    changed = False
    for c in list(cells):
        if c not in cells: continue
        x, y, s = c
        if s == 1: continue
        # neighbours along each edge, sampled at quarter points
        for px, py in [(x-0.5, y+s*0.25), (x-0.5, y+s*0.75), (x+s+0.5, y+s*0.25), (x+s+0.5, y+s*0.75),
                       (x+s*0.25, y-0.5), (x+s*0.75, y-0.5), (x+s*0.25, y+s+0.5), (x+s*0.75, y+s+0.5)]:
            n = leaf_at(px, py)
            if n and n[2] * 2 < s:
                cells.remove(c); cells.update(split(c)); changed = True
                break

# node set: all cell corners
pts = set()
for x, y, s in cells:
    pts.update([(x, y), (x+s, y), (x, y+s), (x+s, y+s)])
def on_edge_mid(p0, p1):
    """Return list of existing nodes strictly between p0 and p1 (axis-aligned), sorted."""
    (x0, y0), (x1, y1) = p0, p1
    res = []
    if y0 == y1:
        res = sorted([p for p in pts if p[1] == y0 and min(x0, x1) < p[0] < max(x0, x1)])
        if x1 < x0: res.reverse()
    else:
        res = sorted([p for p in pts if p[0] == x0 and min(y0, y1) < p[1] < max(y0, y1)], key=lambda p: p[1])
        if y1 < y0: res.reverse()
    return res

block_tris = []          # in units of h, CCW
for x, y, s in sorted(cells):
    corners = [(x, y), (x+s, y), (x+s, y+s), (x, y+s)]
    poly = []
    for k in range(4):
        p0, p1 = corners[k], corners[(k+1) % 4]
        poly.append(p0); poly += on_edge_mid(p0, p1)
    if len(poly) == 4:                                   # plain cell: "/" diagonal
        block_tris += [((x, y), (x+s, y), (x+s, y+s)), ((x, y), (x+s, y+s), (x, y+s))]
    else:                                                # transition cell: centre fan
        c = (x + s/2, y + s/2)
        for k in range(len(poly)):
            block_tris.append((poly[k], poly[(k+1) % len(poly)], c))

def toxy(p): return (round(X0 + p[0]*h, 6), round(p[1]*h, 6))
block_tris = [tuple(toxy(p) for p in t) for t in block_tris]

# block boundary nodes (left x=0, top y=16, right x=32) in physical coordinates, CCW order
bpts = sorted(pts)
left  = sorted([toxy(p) for p in pts if p[0] == 0],  key=lambda q: q[1])
top   = sorted([toxy(p) for p in pts if p[1] == NY], key=lambda q: q[0])
right = sorted([toxy(p) for p in pts if p[0] == NX], key=lambda q: q[1])

# ---------------------------------------------------------------- outer region
# boundary node positions measured from Newman's Figure 4 (upper view)
yleft  = [0, .15, .30, .55, .80, 1.05, 1.35, 1.65, 2.15, 2.75, 3.35, 4.15, 5.10, 6.20, 7.45, 10.0]
xbot_l = [0, .2333, .4667]                       # crack line, x = 0 .. block
xbot_r = [1.71, 2.01, 2.51, 3.0, 3.6, 4.2, 4.85, 5.76, 7.01, 9.0]
xtop   = [0, 3.0, 6.0, 9.0]
yright = [0, 2.5, 5.0, 7.5, 10.0]

loop = []                                        # CCW around the outer region
loop += [(x, 0.0) for x in xbot_l]               # bottom, left of block
loop += left                                     # block left side (x = 0.70), going up
loop += top                                      # block top (y = 0.40), going right
loop += list(reversed(right))                    # block right side (x = 1.50), going down
loop += [(x, 0.0) for x in xbot_r]               # bottom, right of block
loop += [(B, y) for y in yright[1:]]             # right edge, going up
loop += [(x, H) for x in reversed(xtop[:-1])]    # top edge, going left
loop += [(0.0, y) for y in reversed(yleft[1:-1])]# left edge, going down
# remove consecutive duplicates
L2 = []
for p in loop:
    p = (round(p[0], 6), round(p[1], 6))
    if not L2 or p != L2[-1]: L2.append(p)
if L2[0] == L2[-1]: L2.pop()
V = np.array(L2); nV = len(V)
S = np.array([[i, (i+1) % nV] for i in range(nV)])

def target_size(x, y):
    d = math.hypot(x - 1.1, y)
    return 0.20 + 0.18*d + 0.012*y*y

mesh = tr.triangulate({'vertices': V, 'segments': S}, 'pq28Y')
for it in range(8):
    P, T = mesh['vertices'], mesh['triangles']
    cen = P[T].mean(axis=1)
    amax = np.array([math.sqrt(3)/4*target_size(*c)**2 for c in cen])
    e1 = P[T[:, 1]] - P[T[:, 0]]; e2 = P[T[:, 2]] - P[T[:, 0]]
    area = 0.5*np.abs(e1[:, 0]*e2[:, 1] - e1[:, 1]*e2[:, 0])
    if np.all(area <= amax*1.05): break
    mesh = tr.triangulate({'vertices': P, 'segments': mesh['segments'], 'triangles': T,
                           'triangle_max_area': amax}, 'rpq28Ya')
P, T = mesh['vertices'], mesh['triangles']
outer_tris = [tuple((round(P[i][0], 4), round(P[i][1], 4)) for i in t) for t in T]

# ---------------------------------------------------------------- assemble, orient, order
tris = []
for t in block_tris + outer_tris:
    (x1, y1), (x2, y2), (x3, y3) = t
    ar = 0.5*((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    if ar < 0: t = (t[0], t[2], t[1]); ar = -ar
    assert ar > 1e-6, t
    tris.append(t)
key = lambda p: (round(p[0], 4), round(p[1], 4))
nodes = {}
for t in tris:
    for p in t: nodes.setdefault(key(p), len(nodes))
conn = np.array([[nodes[key(p)] for p in t] for t in tris])
nn = len(nodes)
r, c = [], []
for e in conn:
    for i in e:
        for j in e: r.append(i); c.append(j)
G = coo_matrix((np.ones(len(r)), (r, c)), shape=(nn, nn)).tocsr()
perm = reverse_cuthill_mckee(G, symmetric_mode=True)
rank = np.empty(nn, int); rank[perm] = np.arange(nn)
order = sorted(range(len(tris)), key=lambda e: (rank[conn[e]].min(), rank[conn[e]].max()))
tris = [tris[e] for e in order]

# elements per node (program limit MAXEPN)
cnt = np.bincount(conn.ravel(), minlength=nn)
print(f"{len(tris)} elements, {nn} nodes, max elements at a node = {cnt.max()}", file=sys.stderr)

# ---------------------------------------------------------------- input deck
NE = len(tris)
L = []
L.append("%5d%5d%5d%5d%5d%5d" % (0, NE, 0, 4, 4, 9999))
L.append("%10.4f%10.4f" % (A0, B))
bnd = [((0, 0, 0, 1), (0, 0, H, 0), (1, 0), (0, 0)),        # 1: x = 0   symmetry, u = 0
       ((0, H, 0, 0), (B, 0, H, H), (0, 2), (0, 1)),        # 2: y = H   uniform stress S
       ((0, B, 0, 1), (B, B, H, 0), (0, 0), (0, 0)),        # 3: x = b   free
       ((0, 0, 0, 0), (B, 0, 0, 0), (0, 1), (0, 0))]        # 4: y = 0   crack line
for cf, mm, nb, sb in bnd:
    L += ["".join("%10.4f" % v for v in cf), "".join("%10.4f" % v for v in mm),
          "%5d%5d" % nb, "%10.4f%10.4f" % sb]
L += ["%10.4f%10.4f" % (10000.0, 0.3333), "%5d" % 0]
for k, t in enumerate(tris, 1):
    L.append("%5d" % k + "".join("%10.4f" % p[0] for p in t) + "".join("%10.4f" % p[1] for p in t))
L += ["%10.4f%10.4f" % (50.0, 0.0),                           # SIGYS, AM
      "%10.3E" % 1.5,                                         # G (relaxation)
      "%10.3E%10.3E%5d%5d%5d%5d%5d" % (0.05, 0.001, 500, 1, 1, 1, NE)]
for _ in range(NCYC):
    L.append("%10.3E %-4s" % (SMAX, "GROW"))
    L.append("%10.3E %-4s" % (0.001, ""))
L.append("%10.3E %-4s" % (SMAX, ""))
L.append("%10.3E %-4s" % (0.001, "HALT"))
open(OUT, 'w').write("\n".join(L) + "\n")

# ---------------------------------------------------------------- picture
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 11.5), gridspec_kw={'height_ratios': [10, 5.2]})
for ax, (xl, yl) in ((a1, ((0, B), (0, H))), (a2, ((X0, X0+32*h), (0, 16*h)))):
    ax.add_collection(PolyCollection([list(t) for t in tris], facecolors='none', edgecolors='#0b0b0b', linewidths=0.5 if ax is a1 else 0.8))
    ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect('equal')
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
a1.plot([0, A0], [0, 0], color='#eb6834', lw=3)
a1.add_patch(plt.Rectangle((X0, 0), 32*h, 16*h, fill=False, ec='#2a78d6', lw=1.2))
a1.set_title(f"Centre-crack panel, quarter model: {NE} elements, {nn} nodes", loc='left', fontsize=10)
a1.set_xlabel('x, in.'); a1.set_ylabel('y, in.')
a2.plot([X0, A0], [0, 0], color='#eb6834', lw=4)
for k, ch in enumerate('ABCDEF'):
    xa = A0 + k*h
    a2.plot(xa, 0, 'o', color='#0b0b0b', ms=4); a2.text(xa, -0.022, ch, ha='center', va='top', fontsize=9)
a2.set_title("Crack-tip region (blue box above), h = 0.025 in.", loc='left', fontsize=10)
a2.set_xlabel('x, in.'); a2.set_ylabel('y, in.')
fig.tight_layout(); fig.savefig(OUT.replace('.dat', '_mesh.png'), dpi=150)
