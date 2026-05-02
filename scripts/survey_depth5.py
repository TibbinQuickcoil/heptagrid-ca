import math, random
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
    nbrs.append([-1]*P)
    # higher precision key for depth 5
    pm_[(round(cr*15000),round(ci*15000))]=cid; return cid

# Build to depth 5 with tighter cutoff
DEPTH_LIMIT = 5
RADIUS_CUTOFF = 0.989  # higher cutoff for depth 5

pm_={}; add_cell(0,0,0,0,pm_); q=deque([0])
while q:
    cid=q.popleft(); c=cells[cid]
    if c['depth']>=DEPTH_LIMIT: continue
    for j in range(P):
        if nbrs[cid][j]!=-1: continue
        th=c['orient']+TAU*j/P; cr,ci=pm(c['cr'],c['ci'],th,RHO_CC)
        if cr*cr+ci*ci>RADIUS_CUTOFF**2: continue
        key=(round(cr*15000),round(ci*15000))
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

# Verify graph integrity
depth_counts = {}
for c in cells: depth_counts[c['depth']] = depth_counts.get(c['depth'],0)+1
total_edges = sum(1 for cid in range(N) for j in range(P) if nbrs[cid][j]>=0)
print(f"Depth-{DEPTH_LIMIT} grid built: {N} cells")
print(f"Cells per depth: {dict(sorted(depth_counts.items()))}")
print(f"Total directed edges: {total_edges}  (avg {total_edges/N:.2f} per cell)")
print()

def parse_rule(s):
    b,sv=s.split('_'); return set(map(int,b)),set(map(int,sv))
def step_ca(birth,surv):
    nx=[]
    for cid in range(N):
        a=sum(1 for j in range(P) if nbrs[cid][j]>=0 and cells[nbrs[cid][j]]['state'])
        c=cells[cid]; nx.append((1 if surv.__contains__(a) else 0) if c['state'] else (1 if birth.__contains__(a) else 0))
    for cid in range(N):
        c=cells[cid]; c['age']=(c['age']+1 if c['state'] else 0) if nx[cid] else 0; c['state']=nx[cid]

SECTOR=TAU/P; BASELINE=(SECTOR/4)*180/PI
def compute_align():
    far=[c for c in cells if c['state'] and (c['cr']**2+c['ci']**2)>0.002]
    if not far: return None
    total=0
    for c in far:
        theta=math.atan2(c['ci'],c['cr']); t=((theta%SECTOR)+SECTOR)%SECTOR
        total+=min(t,SECTOR-t)
    return (total/len(far))*180/PI

def seed_random(density,seed_val):
    rng=random.Random(seed_val)
    for c in cells: c['state']=1 if rng.random()<density else 0; c['age']=0

# === PRIMARY TEST: B2/S23 across 20 seeds at depth 5 ===
print("="*72)
print(f"DEPTH-{DEPTH_LIMIT} VALIDATION: B2/S23 transient geodesic alignment")
print(f"Grid: {N} cells  |  Baseline: {BASELINE:.3f}°  |  Threshold: 7.071°")
print("="*72)

birth,surv=parse_rule('2_23')
mins=[]; mins_at=[]
for seed in range(20):
    seed_random(0.28, seed)
    aligns=[]
    for g in range(150):  # longer run for larger grid
        step_ca(birth,surv); aligns.append(compute_align())
    valid=[(g,a) for g,a in enumerate(aligns) if a is not None]
    if valid:
        mn = min(valid, key=lambda x: x[1])
        mins.append(mn[1]); mins_at.append(mn[0]+1)
        print(f"  Seed {seed:2d}: min={mn[1]:.2f}° at gen {mn[0]+1:>3}")

print()
mean_min = sum(mins)/len(mins)
below_base = sum(1 for m in mins if m < BASELINE)
below_thresh = sum(1 for m in mins if m < 7.071)
print(f"  Mean minimum: {mean_min:.3f}°")
print(f"  Reduction from baseline: {((BASELINE-mean_min)/BASELINE*100):.1f}%")
print(f"  Below baseline:  {below_base}/20")
print(f"  Below threshold: {below_thresh}/20")
print()
print(f"  DEPTH-4 COMPARISON:")
print(f"    depth-4 mean min: 7.998°  →  depth-5 mean min: {mean_min:.3f}°")
print(f"    depth-4 reduction: 37.8%  →  depth-5 reduction: {((BASELINE-mean_min)/BASELINE*100):.1f}%")

# === STRUCTURAL VALIDATION: does the 7.5x/10x ratio persist? ===
print()
print("="*72)
print("STRUCTURAL VALIDATION: geodesic neighborhood concentration at depth 5")
print("="*72)

GEO_CUTOFF = (SECTOR/4)*180/PI * 0.5

def angular_dist(cr,ci):
    if cr**2+ci**2<0.002: return None
    theta=math.atan2(ci,cr); t=((theta%SECTOR)+SECTOR)%SECTOR
    return min(t,SECTOR-t)*180/PI

for c in cells:
    d = angular_dist(c['cr'],c['ci'])
    c['geo'] = (d is not None and d < GEO_CUTOFF)
    c['ang_dist'] = d

print(f"\n{'Depth':>5} {'Class':>13} {'n':>4} {'GeoNbrs':>9} {'Ratio vs off':>14}")
print("─"*55)
for depth in [1,2,3,4,5]:
    on  = [c for c in cells if c['depth']==depth and c['geo']]
    off = [c for c in cells if c['depth']==depth and not c['geo'] and c['ang_dist'] is not None]
    if not on or not off: continue
    def gnc(c): return sum(1 for j in range(P) if nbrs[c['id']][j]>=0 and cells[nbrs[c['id']][j]]['geo'])
    on_m  = sum(gnc(c) for c in on)/len(on)
    off_m = sum(gnc(c) for c in off)/len(off)
    ratio = on_m/off_m if off_m > 0 else float('inf')
    print(f"{depth:>5} {'on-geodesic':>13} {len(on):>4} {on_m:>8.3f}  {ratio:>10.2f}x")
    print(f"{'':>5} {'off-geodesic':>13} {len(off):>4} {off_m:>8.3f}")

# === B34/S234 oscillator at depth 5 ===
print()
print("="*72)
print("B34/S234 PERIOD CHECK at depth 5")
print("="*72)
birth,surv=parse_rule('34_234')
seed_random(0.28, 42)
counts=[]
for g in range(120):
    step_ca(birth,surv); counts.append(sum(c['state'] for c in cells))
print(f"Counts gens 80-110: {counts[79:110]}")
late = counts[80:]
for period in [2,3,4,5,6,7,8,12,16]:
    if all(late[i]==late[i+period] for i in range(len(late)-period-1)):
        print(f"Period {period} confirmed. Cycle: {late[:period]}")
        break
else:
    print("No simple periodicity detected within tested range — may not converge to oscillator at this scale")
