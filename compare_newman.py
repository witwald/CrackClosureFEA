import re, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ------------------------------------------------------------ published values (digitised)
# Figures 7 and 8: closure (c) / opening (o) stresses, ksi
NEW = {15: dict(yield_=5.73, c=[3.97, 4.24, 4.17, 4.19], o=[3.70, 4.24, 4.27]),
       20: dict(yield_=None, c=[7.84, 8.35, 8.33, 8.33], o=[6.10, 7.79, 8.50])}
# Figure 9: v (1e-3 in.) at Smax = 20 ksi, crack tip at A..E;  x in inches
def f9(px): return 1.075 + (px-955)/1250.0
def v9(py): return (926-py)/102.5
F9cols = {690:[604.5,617.5,631,642.5,653], 750:[621.5,636.5,649.5,662,676], 800:[635.5,651.5,667,678.5,695],
          850:[653.5,670.5,687.5,700.5,719], 890:[669.5,686,704.5,719.5,740.5], 922:[684,703.5,722.5,739.5,769],
          938:[734.5,756,778,808,839.5], 955:[790.5,814,835,882.5,926], 971:[808.5,831,862.5,905],
          987:[823,847,891,926], 1003:[836,869,908], 1018:[849,888.5,926], 1034:[869.5,906], 1050:[890,926],
          1065:[906], 1080:[926]}
NEW9 = {t: [] for t in 'EDCBA'}
for px, ys in sorted(F9cols.items()):
    for t, py in zip('EDCBA', ys): NEW9[t].append((f9(px), v9(py)))
# Figure 10: v during unloading, tip at E (x = 1.175)
def f10(px): return 1.075 + (px-938)/1340.0
def v10(py): return (953.5-py)/111.2
F10 = {20.0: [(640,602),(770,636),(850,667),(895,687),(902,688),(938,809),(971,839.5),(1003,865.5),(1034,904.5),(1072,953.5)],
       13.0: [(640,703),(770,722.5),(850,738.5),(895,750),(905,765.5),(938,865.5),(971,886),(1003,902),(1034,930),(1072,953.5)],
        8.5: [(640,762),(770,776.5),(850,786.5),(895,793.5),(905,801.5),(938,902),(971,919),(1003,931),(1034,949),(1072,953.5)],
        5.0: [(640,831),(770,833),(850,837),(895,839),(905,845),(938,937),(971,946.5),(1003,949.5),(1072,953.5)],
        0.0: [(640,887.5),(770,877),(850,868),(895,860.5),(905,866.5),(938,950.5),(971,953.5),(1072,953.5)]}
NEW10 = {s: [(f10(a), v10(b)) for a, b in pts] for s, pts in F10.items()}

# ------------------------------------------------------------ this program's results
def parse(fn):
    txt = open(fn).read()
    ev = re.findall(r'NODE\s+(\d+)\s+(CLOSED|OPENED|BROKE) AT\s+(-?[\d.]+)', txt)
    blocks = []
    for b in txt.split('APPLIED GROSS STRESS =')[1:]:
        head = b.splitlines()[0]
        S = float(head.split()[0]); a = float(re.search(r'CRACK =\s*([\d.]+)', head).group(1))
        cod = []
        for l in b.split('DISPLACEMENTS')[0].splitlines()[1:]:
            m = re.match(r'\s+(\d+)\s+([\d.]+)\s+(-?[\d.]+E[-+]\d+)\s*$', l)
            if m: cod.append((float(m.group(2)), float(m.group(3))*1e3))
        blocks.append(dict(S=S, a=a, cod=sorted(cod)))
    return ev, blocks
def events(ev):
    """closure/opening stress of the node that was the crack tip, per cycle"""
    order = [n for n, k, s in ev if k == 'BROKE']
    c, o = [], []
    for n in order:
        cl = [float(s) for m, k, s in ev if m == n and k == 'CLOSED']
        op = [float(s) for m, k, s in ev if m == n and k == 'OPENED']
        if cl: c.append(cl[0])
        if op: o.append(op[0])
    return c, o
OURS = {}
for S in (15, 20):
    ev, bl = parse(f'trace_S{S}.out')
    c, o = events(ev)
    y0 = float(re.search(r'INITIAL YIELD IS\s+([\d.]+)', open(f'trace_S{S}.out').read()).group(1))
    OURS[S] = dict(c=c[:4], o=o[:3], yield_=y0, blocks=bl)

# ------------------------------------------------------------ styling
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
PUB, NOW = '#0b0b0b', '#eb6834'
RAMP = ['#86b6ef', '#5598e7', '#2a78d6', '#1c5cab', '#104281']
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK2, 'xtick.color': INK2,
                     'ytick.color': INK2, 'figure.facecolor': SURF, 'axes.facecolor': SURF})
def clean(ax):
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    ax.grid(axis='y', color=GRID, lw=0.6); ax.set_axisbelow(True)

fig = plt.figure(figsize=(12, 10.2))
gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.22)

# ---- (a),(b) closure / opening on the load history, like Figures 7 and 8
for k, S in enumerate((15, 20)):
    ax = fig.add_subplot(gs[0, k])
    t = [0]; s = [0]
    for i in range(4): t += [i + 0.5, i + 1]; s += [S, 0]
    ax.plot(t, s, color=INK2, lw=1)
    lab = 'ABCD'
    def put(vals, phase, pub):
        for i, v in enumerate(vals):
            if v is None: continue
            # unloading leg of cycle i+1 (closure) or loading leg of cycle i+2 (opening)
            if phase == 'c': x = i + 0.5 + 0.5*(1 - v/S)
            else:            x = i + 1 + 0.5*(v/S)
            dx = -0.035 if pub else 0.035
            ax.plot(x + dx, v, 'o' if pub else 'D', ms=7 if pub else 6.5,
                    mfc=(PUB if pub else NOW) if phase == 'c' else SURF,
                    mec=PUB if pub else NOW, mew=1.5, zorder=5)
            if pub: ax.text(x - 0.11 if phase == 'c' else x + 0.1, v + 0.9, lab[i], ha='center', fontsize=8.5, color=INK)
    put(NEW[S]['c'], 'c', True); put(NEW[S]['o'], 'o', True)
    put(OURS[S]['c'], 'c', False); put(OURS[S]['o'], 'o', False)
    if NEW[S]['yield_']:
        ax.plot(0.5*NEW[S]['yield_']/S - 0.02, NEW[S]['yield_'], 'x', color=PUB, ms=7, mew=1.5)
    ax.plot(0.5*OURS[S]['yield_']/S + 0.02, OURS[S]['yield_'], 'x', color=NOW, ms=7, mew=1.5)
    so_pub = NEW[S]['c'][-1]; so_now = OURS[S]['o'][-1]
    ax.set_xlim(-0.05, 4.6); ax.set_ylim(0, S*1.25); ax.set_xticks([])
    ax.set_xlabel('Time'); ax.set_ylabel('S, ksi')
    ax.set_title(f"({'ab'[k]}) Smax = {S} ksi, R = 0  (Newman Fig. {7 if S == 15 else 8})", loc='left', color=INK, fontsize=10)
    clean(ax)
    txt = (f"Stabilised opening/closure:\n  Newman  ≈ {so_pub:.1f} ksi ({100*so_pub/S:.0f}% Smax)\n"
           f"  This run = {OURS[S]['o'][-1]:.1f} / {OURS[S]['c'][-1]:.1f} ksi")
    ax.text(4.55, S*1.22, txt, ha='right', va='top', fontsize=8, color=INK2, family='monospace')

# ---- (c) crack-surface displacement at Smax = 20 ksi, tip at A..E (Figure 9)
ax = fig.add_subplot(gs[1, 0])
bl = OURS[20]['blocks']
smax = max(b['S'] for b in bl)
for i, t in enumerate('ABCDE'):
    a = 1.075 + 0.025*i
    cand = [b for b in bl if abs(b['S'] - smax) < 1e-6 and abs(b['a'] - a) < 1e-4]
    b = cand[-1]
    xs = [x for x, v in b['cod'] if x >= 0.84]; vs = [v for x, v in b['cod'] if x >= 0.84]
    ax.plot(xs, vs, '-', color=RAMP[i], lw=2, marker='o', ms=3.5, zorder=3)
    nx, nv = zip(*NEW9[t])
    ax.plot(nx, nv, '--', color=RAMP[i], lw=1.4, zorder=2)
    ax.text(a, -0.13, t, ha='center', va='top', fontsize=8.5, color=INK)
ax.set_xlim(0.84, 1.19); ax.set_ylim(-0.25, 4); ax.set_xlabel('x, in.'); ax.set_ylabel('v × 10³, in.')
ax.set_title('(c) Crack-surface displacement at Smax = 20 ksi, tip at A…E  (Fig. 9)', loc='left', color=INK, fontsize=10)
clean(ax)
ax.text(0.85, 3.75, 'curves: crack tip at A (lightest) … E (darkest)', fontsize=8, color=INK2)

# ---- (d) unloading with tip at E (Figure 10)
ax = fig.add_subplot(gs[1, 1])
unl = [b for b in bl if abs(b['a'] - 1.175) < 1e-4]
k0 = next(k for k, b in enumerate(unl) if abs(b['S'] - smax) < 1e-6)
k1 = next(k for k in range(k0, len(unl)) if unl[k]['S'] < smax - 1e-6)
leg = unl[k1 - 1:]                       # Smax state with tip at E, then unloading
kend = next((k for k in range(1, len(leg)) if leg[k]['S'] > leg[k-1]['S']), len(leg))
leg = leg[:kend]
for i, St in enumerate((20.0, 13.0, 8.5, 5.0, 0.0)):
    b = leg[0] if St == 20.0 else min(leg, key=lambda q: abs(q['S'] - St))
    xs = [x for x, v in b['cod'] if x >= 0.84]; vs = [v for x, v in b['cod'] if x >= 0.84]
    col = RAMP[4 - i]
    ax.plot(xs, vs, '-', color=col, lw=2, marker='o', ms=3.5, zorder=3)
    nx, nv = zip(*NEW10[St])
    ax.plot(nx, nv, '--', color=col, lw=1.4, zorder=2)
    ax.text(0.845, vs[0] + 0.06, f"S = {St:g} ksi" + (f"  (run: {b['S']:.2f})" if abs(b['S']-St) > 0.05 else ''), fontsize=8, color=INK, va='bottom')
for i, t in enumerate('ABCDE'):
    ax.text(1.075 + 0.025*i, -0.13, t, ha='center', va='top', fontsize=8.5, color=INK)
ax.set_xlim(0.84, 1.19); ax.set_ylim(-0.25, 4); ax.set_xlabel('x, in.'); ax.set_ylabel('v × 10³, in.')
ax.set_title('(d) Unloading from Smax = 20 ksi with tip at E  (Fig. 10)', loc='left', color=INK, fontsize=10)
clean(ax)

h = [Line2D([], [], ls='', marker='o', mfc=PUB, mec=PUB, ms=7, label='Newman 1974: closure'),
     Line2D([], [], ls='', marker='o', mfc=SURF, mec=PUB, mew=1.5, ms=7, label='Newman 1974: opening'),
     Line2D([], [], ls='', marker='D', mfc=NOW, mec=NOW, ms=6.5, label='This F77 run: closure'),
     Line2D([], [], ls='', marker='D', mfc=SURF, mec=NOW, mew=1.5, ms=6.5, label='This F77 run: opening'),
     Line2D([], [], ls='', marker='x', color=INK2, mew=1.5, label='initial yield'),
     Line2D([], [], color=INK2, lw=1.4, ls='--', label='Newman 1974 (digitised)'),
     Line2D([], [], color=INK2, lw=2, marker='o', ms=3.5, label='This F77 run')]
fig.legend(handles=h, loc='upper center', ncol=7, frameon=False, fontsize=8, bbox_to_anchor=(0.5, 0.995))
fig.savefig('newman_comparison.png', dpi=160, bbox_inches='tight')

# table for the record
for S in (15, 20):
    print(S, 'closure ours', OURS[S]['c'], 'pub', NEW[S]['c'])
    print(S, 'opening ours', OURS[S]['o'], 'pub', NEW[S]['o'])
