import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch
from matplotlib.lines import Line2D
exec(open('pz.py').read().split('for S in (15,20):')[0])
INK,INK2,SURF='#0b0b0b','#52514e','#fcfcfb'; FWD,REV,ELA='#eb6834','#2a78d6','#efeeea'
plt.rcParams.update({'font.size':9,'axes.edgecolor':INK2,'axes.labelcolor':INK2,'xtick.color':INK2,'ytick.color':INK2,'figure.facecolor':SURF,'axes.facecolor':SURF})
XL,YL=(1.0,1.45),(0,0.20)
fig,axs=plt.subplots(2,3,figsize=(14,5.6),sharex=True,sharey=True)
titles=['First peak, tip at A','4th peak, after release of D (tip at E)','Unloaded to 0, tip at E']
for r,S in enumerate((15,20)):
    els,st=load(S)
    for c,k in enumerate((0,10,11)):
        ax=axs[r,c]; s=st[k]; ny=0; nr=0
        for e,(x,y) in enumerate(els,1):
            if min(x)>XL[1]+0.1 or min(y)>YL[1]+0.1 or max(x)<XL[0]-0.1: continue
            v=s['el'][e]
            if v['se']>=49.5 and v['sy']>=0: fc,h=FWD,None; ny+=1
            elif v['se']>=49.5: fc,h=REV,'////'; nr+=1
            else: fc,h=ELA,None
            ax.add_patch(Polygon(list(zip(x,y)),closed=True,fc=fc,ec='#ffffff' if fc!=ELA else '#c9c8c3',lw=0.6,hatch=h))
        a=s['a']
        ax.plot([XL[0],a],[0,0],color=INK,lw=3.5,solid_capstyle='butt',zorder=5)
        if c==2:
            cl=[x for x,v in s['cod'] if abs(v)<1e-8 and x<a-1e-4 and x>=XL[0]]
            ax.plot(cl,[0]*len(cl),'o',ms=6,mfc=SURF,mec=INK,mew=1.5,zorder=7)
        for i,ch in enumerate('ABCDEF'):
            ax.text(1.075+0.025*i,-0.005,ch,ha='center',va='top',fontsize=7.5,color=INK)
        ax.tick_params(axis='x',pad=13); ax.set_xticks([1.0,1.1,1.2,1.3,1.4]); ax.set_xlim(*XL); ax.set_ylim(*YL); ax.set_aspect('equal')
        for sp in ('top','right'): ax.spines[sp].set_visible(False)
        n=(f"{ny} element{'s'*(ny!=1)} at yield" if c<2 else f"{nr} element{'s'*(nr!=1)} in compressive yield")
        ax.set_title(f"Smax = {S} ksi — {titles[c]}\nS = {abs(s['S']) if abs(s['S'])>0.05 else 0:.1f} ksi, a = {a:.3f} in.;  {n}",loc='left',fontsize=9,color=INK)
        if r==1: ax.set_xlabel('x, in.')
        if c==0: ax.set_ylabel('y, in.')
h=[Patch(fc=FWD,ec='w',label='at yield, tensile σy'),Patch(fc=REV,ec='w',hatch='////',label='at yield, compressive σy (reverse yielding)'),
   Patch(fc=ELA,ec='#c9c8c3',label='elastic'),Line2D([],[],color=INK,lw=3.5,label='crack'),
   Line2D([],[],marker='o',ls='',mfc=SURF,mec=INK,mew=1.5,label='closed crack-face node')]
fig.legend(handles=h,loc='upper center',ncol=5,frameon=False,fontsize=8.5,bbox_to_anchor=(0.5,1.0))
fig.tight_layout(rect=(0,0,1,0.95))
fig.savefig('fig4_mesh_plastic_zones.png',dpi=160)
