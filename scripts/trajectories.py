import math, random, json
from collections import deque

PI=math.pi; TAU=2*PI; P=7
APO=math.acosh(math.cos(PI/3)/math.sin(PI/7))
HE=math.acosh(math.cos(PI/7)/math.sin(PI/3))
RC=math.acosh(math.cosh(APO)*math.cosh(HE))
RHO_CC=math.tanh(APO)

def pm(zr,zi,th,rho):
    pr,pi_=rho*math.cos(th),rho*math.sin(th); nr,ni=pr+zr,pi_+zi
    dr=1+zr*pr+zi*pi_; di=zr*pi_-zi*pr; d2=dr*dr+di*di
    return (nr*dr+ni*di)/d2,(ni*dr-nr*di)/d2
def hd(x1,y1,x2,y2):
    nr,ni=x2-x1,y2-y1; dr=1-x1*x2-y1*y2; di=-(x1*y2-y1*x2); d2=dr*dr+di*di
    return math.atan2((ni*dr-nr*di)/d2,(nr*dr+ni*di)/d2)

cells=[]; nbrs=[]
def add_cell(cr,ci,orient,depth,pm_):
    cid=len(cells)
    cells.append({'id':cid,'cr':cr,'ci':ci,'orient':orient,'depth':depth,'state':0,'age':0})
    nbrs.append([-1]*P); pm_[(round(cr*7000),round(ci*7000))]=cid; return cid

pm_={}; add_cell(0,0,0,0,pm_); q=deque([0])
while q:
    cid=q.popleft(); c=cells[cid]
    if c['depth']>=4: continue
    for j in range(P):
        if nbrs[cid][j]!=-1: continue
        th=c['orient']+TAU*j/P; cr,ci=pm(c['cr'],c['ci'],th,RHO_CC)
        if cr*cr+ci*ci>0.971**2: continue
        key=(round(cr*7000),round(ci*7000))
        if key in pm_:
            eid=pm_[key]; nbrs[cid][j]=eid; ec=cells[eid]
            bd=hd(ec['cr'],ec['ci'],c['cr'],c['ci']); bj,bd2=0,9.0
            for jj in range(P):
                d=bd-(ec['orient']+TAU*jj/P); d-=TAU*round(d/TAU)
                if abs(d)<bd2: bd2,bj=abs(d),jj
            if nbrs[eid][bj]==-1: nbrs[eid][bj]=cid
        else:
            bd=hd(cr,ci,c['cr'],c['ci']); nid=add_cell(cr,ci,bd,c['depth']+1,pm_)
            nbrs[nid][0]=cid; nbrs[cid][j]=nid; q.append(nid)
N=len(cells)

SECTOR=TAU/P; BASELINE=(SECTOR/4)*180/PI
GEO_CUTOFF=BASELINE*0.55

def parse_rule(s):
    b,sv=s.split('_'); return set(map(int,b)),set(map(int,sv))
def step_ca(birth,surv):
    nx=[]
    for cid in range(N):
        a=sum(1 for j in range(P) if nbrs[cid][j]>=0 and cells[nbrs[cid][j]]['state'])
        c=cells[cid]; nx.append((1 if surv.__contains__(a) else 0) if c['state'] else (1 if birth.__contains__(a) else 0))
    for cid in range(N):
        c=cells[cid]; c['age']=(c['age']+1 if c['state'] else 0) if nx[cid] else 0; c['state']=nx[cid]
def compute_align():
    far=[c for c in cells if c['state'] and (c['cr']**2+c['ci']**2)>0.002]
    if not far: return None
    total=0
    for c in far:
        theta=math.atan2(c['ci'],c['cr']); t=((theta%SECTOR)+SECTOR)%SECTOR
        total+=min(t,SECTOR-t)
    return (total/len(far))*180/PI
def seed_random(d,s):
    rng=random.Random(s)
    for c in cells: c['state']=1 if rng.random()<d else 0; c['age']=0

# Generate full 100-gen trajectories for 5 representative seeds
birth,surv=parse_rule('2_23')
seeds_to_plot=[0,5,7,11,17]  # diverse: includes the two extreme minima
trajectories={}
for seed in seeds_to_plot:
    seed_random(0.28,seed)
    traj=[]
    for g in range(100):
        step_ca(birth,surv); traj.append(compute_align())
    trajectories[seed]=traj

# Also generate B1/S12 trajectory (null result) for seed 0
birth2,surv2=parse_rule('1_12')
seed_random(0.28,0)
null_traj=[]
for g in range(100):
    step_ca(birth2,surv2); null_traj.append(compute_align())

# Export cell positions for figure 1
cell_data=[]
for c in cells:
    if c['cr']**2+c['ci']**2<0.002: continue
    theta=math.atan2(c['ci'],c['cr'])
    t=((theta%SECTOR)+SECTOR)%SECTOR
    ang_dist=min(t,SECTOR-t)*180/PI
    is_geo=ang_dist<GEO_CUTOFF*0.5  # tighter for figure clarity
    cell_data.append({'cr':c['cr'],'ci':c['ci'],'depth':c['depth'],'geo':is_geo,'ang':ang_dist})

data={
    'baseline':BASELINE,
    'threshold':GEO_CUTOFF,
    'trajectories':trajectories,
    'null_traj':null_traj,
    'cells':cell_data,
    'N':N
}
with open('/home/claude/fig_data.json','w') as f:
    json.dump(data,f)
print(f"Generated trajectory data for {len(seeds_to_plot)} seeds + null")
print(f"Cell data: {len(cell_data)} cells (non-center)")
