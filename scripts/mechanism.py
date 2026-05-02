import math
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
    cells.append({'id':cid,'cr':cr,'ci':ci,'orient':orient,'depth':depth})
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

SECTOR=TAU/P

def angular_dist(cr,ci):
    if cr**2+ci**2<0.002: return None
    theta=math.atan2(ci,cr)
    t=((theta%SECTOR)+SECTOR)%SECTOR
    return min(t,SECTOR-t)*180/PI

# ── Label every cell as geodesic or off-geodesic ─────────────────
# Geodesic: angular distance < half a sector / 2 = π/14 * 0.5
# i.e. within the inner quarter of its sector
GEO_CUTOFF = (SECTOR/4)*180/PI * 0.5  # 6.43° — inner half of inner half

for c in cells:
    d = angular_dist(c['cr'],c['ci'])
    c['geo'] = (d is not None and d < GEO_CUTOFF)
    c['ang_dist'] = d

# ── KEY COMPUTATION ───────────────────────────────────────────────
# For each non-center cell, count how many of its 7 neighbors
# are themselves on-geodesic.
# Compare: on-geodesic cells vs off-geodesic cells.

print("="*70)
print("MECHANISM TEST: Neighbor Geodesic Density")
print("="*70)
print()
print("Question: Do on-geodesic cells have more on-geodesic neighbors")
print("than off-geodesic cells of the same depth?")
print("If yes, correlated birth under sparse random IC is more likely")
print("along geodesics — this is the structural mechanism.")
print()

for depth in [1,2,3,4]:
    on_geo  = [c for c in cells if c['depth']==depth and c['geo']]
    off_geo = [c for c in cells if c['depth']==depth and not c['geo']
               and c['ang_dist'] is not None]

    if not on_geo or not off_geo:
        continue

    # For each cell, count neighbors that are on-geodesic
    def geo_neighbor_count(c):
        return sum(1 for j in range(P)
                   if nbrs[c['id']][j]>=0
                   and cells[nbrs[c['id']][j]]['geo'])

    on_counts  = [geo_neighbor_count(c) for c in on_geo]
    off_counts = [geo_neighbor_count(c) for c in off_geo]

    on_mean  = sum(on_counts)/len(on_counts)
    off_mean = sum(off_counts)/len(off_counts)

    print(f"Depth {depth}:")
    print(f"  On-geodesic  cells: n={len(on_geo):3d}  "
          f"mean geo-neighbors = {on_mean:.3f}  counts={sorted(set(on_counts))}")
    print(f"  Off-geodesic cells: n={len(off_geo):3d}  "
          f"mean geo-neighbors = {off_mean:.3f}  counts={sorted(set(off_counts))}")
    print(f"  Ratio on/off: {on_mean/off_mean:.3f}x  "
          f"Δ = {on_mean-off_mean:+.3f}")
    print()

# ── BIRTH PROBABILITY UNDER SPARSE RANDOM IC ─────────────────────
print("="*70)
print("BIRTH PROBABILITY ANALYSIS: B2/S23 under density p")
print("="*70)
print()
print("For B2/S23, a dead cell is born if exactly 2 of 7 neighbors are alive.")
print("P(birth | density p) = C(7,2) * p^2 * (1-p)^5")
print()
print("BUT: if neighbors are positively correlated (co-located on geodesic),")
print("P(exactly 2 alive) is HIGHER than the independent-cell prediction.")
print("This is the structural birth advantage.")
print()

import math as m

# For a cell with k geo-neighbors, under density p,
# expected P(exactly 2 alive) if geo-neighbors are correlated
# Model: geo-neighbors alive together with prob q_geo,
#        non-geo neighbors alive independently with prob p

p = 0.28  # initial density

def p_birth_independent(p, total=7, required=2):
    return m.comb(total, required) * p**required * (1-p)**(total-required)

base = p_birth_independent(p)
print(f"P(birth) independent model, p={p}: {base:.6f}")
print()

# Now compute for on-geodesic cells vs off-geodesic cells
# using actual neighbor geo-counts from our graph
print("Depth-2 cells — actual neighbor structure:")
depth2 = [c for c in cells if c['depth']==2 and c['ang_dist'] is not None]
for c in depth2[:6]:  # show a few examples
    nb_ids = [nbrs[c['id']][j] for j in range(P) if nbrs[c['id']][j]>=0]
    nb_geo = [cells[n]['geo'] for n in nb_ids]
    nb_depths = [cells[n]['depth'] for n in nb_ids]
    print(f"  Cell {c['id']:3d} ang={c['ang_dist']:5.2f}° geo={c['geo']}  "
          f"neighbors: {sum(nb_geo)} geo-adjacent  depths={sorted(nb_depths)}")

print()

# ── SHARED NEIGHBOR ANALYSIS ──────────────────────────────────────
print("="*70)
print("SHARED NEIGHBOR ANALYSIS")
print("="*70)
print()
print("Key structural question: for two cells A,B on the same geodesic")
print("at adjacent depths, how many neighbors do they share?")
print("Shared neighbors = correlated exposure = coordinated birth/survival.")
print()

# For each pair of cells (c1 at depth d, c2 at depth d+1)
# where c2 is a neighbor of c1 and both are on-geodesic:
# count shared neighbors

def shared_neighbors(id1, id2):
    n1 = set(nbrs[id1][j] for j in range(P) if nbrs[id1][j]>=0)
    n2 = set(nbrs[id2][j] for j in range(P) if nbrs[id2][j]>=0)
    return len(n1 & n2)

# On-geodesic pairs (adjacent depths, connected)
geo_pairs_shared = []
offgeo_pairs_shared = []

for c1 in cells:
    if c1['depth'] >= 4: continue
    for j in range(P):
        nid = nbrs[c1['id']][j]
        if nid < 0: continue
        c2 = cells[nid]
        if c2['depth'] != c1['depth']+1: continue
        sh = shared_neighbors(c1['id'], c2['id'])
        if c1['geo'] and c2['geo']:
            geo_pairs_shared.append(sh)
        elif not c1['geo'] and not c2['geo']:
            offgeo_pairs_shared.append(sh)

if geo_pairs_shared:
    gm = sum(geo_pairs_shared)/len(geo_pairs_shared)
    print(f"On-geodesic adjacent pairs:  n={len(geo_pairs_shared)}  "
          f"mean shared neighbors = {gm:.4f}")
    print(f"  Distribution: {sorted(set(geo_pairs_shared))} "
          f"counts={[geo_pairs_shared.count(v) for v in sorted(set(geo_pairs_shared))]}")

if offgeo_pairs_shared:
    om = sum(offgeo_pairs_shared)/len(offgeo_pairs_shared)
    print(f"Off-geodesic adjacent pairs: n={len(offgeo_pairs_shared)}  "
          f"mean shared neighbors = {om:.4f}")
    print(f"  Distribution: {sorted(set(offgeo_pairs_shared))} "
          f"counts={[offgeo_pairs_shared.count(v) for v in sorted(set(offgeo_pairs_shared))]}")

if geo_pairs_shared and offgeo_pairs_shared:
    print(f"\n  Ratio: {gm/om:.4f}x more shared neighbors on-geodesic vs off")
    print(f"  This is a STRUCTURAL property of the {{7,3}} tiling.")
    print(f"  It does not depend on any rule or initial condition.")

print()
print("="*70)
print("CONCLUSION")
print("="*70)
