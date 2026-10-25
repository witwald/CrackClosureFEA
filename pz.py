import re, json, sys
def load(S):
    lines=open(f'ccp_fig4_S{S}.dat').read().split('\n')
    NE=int(lines[0][5:10]); els=[]
    start=2+4*4+2
    for l in lines[start:start+NE]:
        v=[float(l[5+10*k:15+10*k]) for k in range(6)]
        els.append((v[0:3],v[3:6]))
    txt=open(f'ccp_fig4_S{S}.out').read()
    st=[]
    for b in txt.split('APPLIED GROSS STRESS =')[1:]:
        head=b.splitlines()[0]; s=float(head.split()[0]); a=float(re.search(r'CRACK =\s*([\d.]+)',head).group(1))
        d={}
        for l in b.splitlines():
            m=re.match(r'\s+(\d+)(\*?)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+\S+\s+\S+\s+\S+$',l)
            if m:
                e=int(m.group(1)); sx,sy,t=map(float,m.group(3,4,5))
                d[e]=dict(sx=sx,sy=sy,t=t,se=(sx*sx+sy*sy-sx*sy+3*t*t)**.5,flag=bool(m.group(2)))
        cod=[(float(x),float(v)) for x,v in re.findall(r'^\s+\d+\s+([\d.]+)\s+(-?[\d.]+E[-+]\d+)\s*$',b.split('DISPLACEMENTS')[0],re.M)]
        st.append(dict(S=s,a=a,el=d,cod=cod))
    return els,st
for S in (15,20):
    els,st=load(S)
    for k,s in enumerate(st):
        y=[e for e,v in s['el'].items() if v['se']>=49.5]
        c=[e for e in y if s['el'][e]['sy']<0]
        xm=max([max(els[e-1][0]) for e in y if min(els[e-1][1])<1e-6 and s['el'][e]['sy']>0]+[s['a']])
        print(S,k,f"S={s['S']:6.2f} a={s['a']:.3f} yield={len(y)} comp={len(c)} flag={sum(v['flag'] for v in s['el'].values())} fwd_ext={xm-s['a']:.3f}")
for S in (15,20):
    els,st=load(S)
    for k in (0,10,11):
        s=st[k]; y=[e for e,v in s['el'].items() if v['se']>=49.5]
        xs=[x for e in y for x in els[e-1][0]]; ys=[v for e in y for v in els[e-1][1]]
        print(S,k,min(xs),max(xs),max(ys))
