# Graded centre-cracked panel mesh for Newman's FINAL program.
# Quarter model: x = 0..1 (half-width b = 1), y = 0..1, crack a0 = 0.30
xs=[0,.10,.17,.22,.25,.27]+[round(.28+.01*i,4) for i in range(11)]+[.41,.46,.54,.66,.82,1.0]
ys=[0,.01,.02,.035,.06,.10,.17,.30,.55,1.0]
a0=0.30; E=10000.0; nu=0.3333; sys_=50.0
Smax=27.5; ncyc=4
els=[]
for i in range(len(xs)-1):            # column by column -> small bandwidth
    for j in range(len(ys)-1):
        x0,x1=xs[i],xs[i+1]; y0,y1=ys[j],ys[j+1]
        els.append(((x0,x1,x1),(y0,y0,y1)))
        els.append(((x0,x1,x0),(y0,y1,y1)))
NE=len(els)
L=["%5d%5d%5d%5d%5d%5d"%(0,NE,0,4,4,999),"%10.4f%10.4f"%(a0,1.0)]
W=1.0;H=1.0
bnd=[((0,0,0,1),(0,0,H,0),(1,0),(0,0)),     # 1: x=0, symmetry (u=0)
     ((0,H,0,0),(W,0,H,H),(0,2),(0,1)),     # 2: y=H, uniform stress
     ((0,W,0,1),(W,W,H,0),(0,0),(0,0)),     # 3: x=W, free
     ((0,0,0,0),(W,0,0,0),(0,1),(0,0))]     # 4: y=0, crack plane
for c,mm,nb,sb in bnd:
    L+=["".join("%10.4f"%v for v in c),"".join("%10.4f"%v for v in mm),"%5d%5d"%nb,"%10.4f%10.4f"%sb]
L+=["%10.4f%10.4f"%(E,nu),"%5d"%0]
for k,(x,y) in enumerate(els,1):
    L.append("%5d"%k+"".join("%10.4f"%v for v in x+y))
L+=["%10.4f%10.4f"%(sys_,0.0),"%10.3E"%1.5,
    "%10.3E%10.3E%5d%5d%5d%5d%5d"%(0.05,0.001,300,1,10,1,NE)]
for c in range(ncyc):
    L.append("%10.3E %-4s"%(Smax,"GROW"))
    L.append("%10.3E %-4s"%(0.001,""))
L.append("%10.3E %-4s"%(Smax,""))
L.append("%10.3E %-4s"%(0.001,"HALT"))
print("\n".join(L))
import sys; print(NE,"elements,",len(xs)*len(ys),"grid nodes",file=sys.stderr)
