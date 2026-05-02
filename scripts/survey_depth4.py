import math
import json
from collections import deque

PI = math.pi
TAU = 2 * PI
P = 7

# {7,3} hyperbolic constants
APO = math.acosh(math.cos(PI/3) / math.sin(PI/7))
HE  = math.acosh(math.cos(PI/7) / math.sin(PI/3))
RC  = math.acosh(math.cosh(APO) * math.cosh(HE))
RHO_CC = math.tanh(APO)
RHO_V  = math.tanh(RC / 2)

def pm(zr, zi, th, rho):
    pr, pi_ = rho*math.cos(th), rho*math.sin(th)
    nr, ni = pr+zr, pi_+zi
    dr = 1 + zr*pr + zi*pi_
    di = zr*pi_ - zi*pr
    d2 = dr*dr + di*di
    return (nr*dr + ni*di)/d2, (ni*dr - nr*di)/d2

def hd(x1, y1, x2, y2):
    nr, ni = x2-x1, y2-y1
    dr = 1 - x1*x2 - y1*y2
    di = -(x1*y2 - y1*x2)
    d2 = dr*dr + di*di
    return math.atan2((ni*dr - nr*di)/d2, (nr*dr + ni*di)/d2)

# BFS cell generation
cells = []  # {id, cr, ci, orient, depth, state, age}
nbrs  = []

def add_cell(cr, ci, orient, depth, pos_map):
    cid = len(cells)
    cells.append({'id':cid,'cr':cr,'ci':ci,'orient':orient,
                  'depth':depth,'state':0,'age':0})
    nbrs.append([-1]*P)
    key = (round(cr*7000), round(ci*7000))
    pos_map[key] = cid
    return cid

def build_tiling():
    pos_map = {}
    add_cell(0, 0, 0, 0, pos_map)
    q = deque([0])
    while q:
        cid = q.popleft()
        c = cells[cid]
        if c['depth'] >= 4:
            continue
        for j in range(P):
            if nbrs[cid][j] != -1:
                continue
            th = c['orient'] + TAU*j/P
            cr, ci = pm(c['cr'], c['ci'], th, RHO_CC)
            if cr*cr + ci*ci > 0.971**2:
                continue
            key = (round(cr*7000), round(ci*7000))
            if key in pos_map:
                eid = pos_map[key]
                nbrs[cid][j] = eid
                ec = cells[eid]
                bd = hd(ec['cr'], ec['ci'], c['cr'], c['ci'])
                bj, bd2 = 0, 9.0
                for jj in range(P):
                    d = bd - (ec['orient'] + TAU*jj/P)
                    d -= TAU * round(d/TAU)
                    if abs(d) < bd2:
                        bd2, bj = abs(d), jj
                if nbrs[eid][bj] == -1:
                    nbrs[eid][bj] = cid
            else:
                bd = hd(cr, ci, c['cr'], c['ci'])
                nid = add_cell(cr, ci, bd, c['depth']+1, pos_map)
                nbrs[nid][0] = cid
                nbrs[cid][j] = nid
                q.append(nid)

build_tiling()
N = len(cells)
print(f"Tiling built: {N} cells\n")

# CA logic
def parse_rule(s):
    b, sv = s.split('_')
    return set(map(int, b)), set(map(int, sv))

def step_ca(birth, surv):
    nexts = []
    for cid in range(N):
        a = sum(1 for j in range(P) if nbrs[cid][j]>=0 and cells[nbrs[cid][j]]['state'])
        c = cells[cid]
        n = (1 if surv.__contains__(a) else 0) if c['state'] else (1 if birth.__contains__(a) else 0)
        nexts.append(n)
    for cid in range(N):
        c = cells[cid]
        c['age'] = (c['age']+1 if c['state'] else 0) if nexts[cid] else 0
        c['state'] = nexts[cid]

# Geodesic alignment
SECTOR = TAU / P
BASELINE_DEG = (SECTOR/4) * 180/PI   # 12.857...
GEO_THRESHOLD = BASELINE_DEG * 0.55  # 7.071...

def compute_align():
    far = [c for c in cells if c['state'] and (c['cr']**2+c['ci']**2) > 0.002]
    if not far:
        return None
    total = 0
    for c in far:
        theta = math.atan2(c['ci'], c['cr'])
        t = ((theta % SECTOR) + SECTOR) % SECTOR
        total += min(t, SECTOR - t)
    return (total / len(far)) * 180/PI

def count_alive():
    return sum(c['state'] for c in cells)

def seed_default():
    for c in cells:
        c['state'] = 0; c['age'] = 0
    for i in [0,1,2,8]:
        if i < N:
            cells[i]['state'] = 1

def classify(final_alive, alive_series, avg_align):
    if final_alive < 2:
        return "extinction", False
    if avg_align is not None and avg_align < GEO_THRESHOLD and final_alive >= 3:
        return "geodesic chains", True
    last10 = alive_series[50:]
    variance = sum((v-last10[0])**2 for v in last10)/len(last10) if last10 else 999
    if variance < 4 and final_alive > 2:
        return "stable structures", False
    if alive_series[-1] > alive_series[30]*2.5 and alive_series[-1] > 40:
        return "rapid expansion", False
    return "diffuse propagation", False

# --- SURVEY ---
rule_list = [
    ('3_34',   'B3/S34'),
    ('2_34',   'B2/S34'),
    ('23_34',  'B23/S34'),
    ('3_23',   'B3/S23'),
    ('2_23',   'B2/S23'),
    ('34_234', 'B34/S234'),
    ('2_345',  'B2/S345'),
    ('1_12',   'B1/S12'),
]

print(f"{'Rule':<12} {'Alive@60':>8} {'AvgAlign':>9} {'MinAlign':>9} {'vs Baseline':>12}  Behavior")
print("─"*75)

results = []
for key, label in rule_list:
    birth, surv = parse_rule(key)
    seed_default()

    align_series = []
    alive_series = []
    for g in range(60):
        step_ca(birth, surv)
        align_series.append(compute_align())
        alive_series.append(count_alive())

    final_alive = alive_series[-1]
    valid_a = [a for a in align_series if a is not None]
    late_a  = valid_a[max(0,len(valid_a)-35):]
    avg_align = sum(late_a)/len(late_a) if late_a else None
    min_align = min(valid_a) if valid_a else None

    behavior, is_geo = classify(final_alive, alive_series, avg_align)

    if avg_align is not None:
        bias_pct = round((1 - avg_align/BASELINE_DEG)*100)
        bias_str = f"{bias_pct:+d}%"
    else:
        bias_str = "—"

    avg_str = f"{avg_align:.2f}°" if avg_align is not None else "—"
    min_str = f"{min_align:.2f}°" if min_align is not None else "—"
    geo_mark = " ◉" if is_geo else ""

    print(f"{label:<12} {final_alive:>8} {avg_str:>9} {min_str:>9} {bias_str:>12}  {behavior}{geo_mark}")

    results.append({
        'rule': label,
        'alive_at_60': final_alive,
        'avg_align_deg': round(avg_align, 3) if avg_align else None,
        'min_align_deg': round(min_align, 3) if min_align else None,
        'geodesic_bias': is_geo,
        'behavior': behavior,
    })

print()
print(f"Random baseline: {BASELINE_DEG:.3f}°")
print(f"Geodesic threshold (<{GEO_THRESHOLD:.3f}°): {sum(r['geodesic_bias'] for r in results)}/{len(results)} rules")
print()

# Save JSON
out = {
    'experiment': 'Geodesic Alignment Survey — {7,3} Heptagrid',
    'tiling': '{7,3}',
    'neighbors': 7,
    'seed': 'cells [0,1,2,8]',
    'generations': 60,
    'total_cells': N,
    'baseline_deg': round(BASELINE_DEG, 4),
    'geodesic_threshold_deg': round(GEO_THRESHOLD, 4),
    'results': results
}
with open('/home/claude/survey_results.json','w') as f:
    json.dump(out, f, indent=2)
print("Results saved to survey_results.json")

print("\n" + "="*75)
print("DEEP ANALYSIS — per-generation breakdown for flagged rules")
print("="*75)

# Re-run with per-generation detail for B3/S23 and B34/S234
for key, label in [('3_23','B3/S23'),('34_234','B34/S234'),('23_34','B23/S34'),('3_34','B3/S34')]:
    birth, surv = parse_rule(key)
    seed_default()
    print(f"\n{label}:")
    print(f"  {'Gen':>4}  {'Alive':>6}  {'Align':>8}  {'Depths of alive cells'}")
    for g in range(60):
        step_ca(birth, surv)
        alive = count_alive()
        align = compute_align()
        alive_cells = [c for c in cells if c['state']]
        depths = sorted(set(c['depth'] for c in alive_cells))
        if g < 15 or g % 10 == 9 or g == 59:
            align_str = f"{align:.2f}°" if align is not None else "  none"
            print(f"  {g+1:>4}  {alive:>6}  {align_str:>8}  depths={depths}")

print()
print("NOTE: Cells at depth=1 sit exactly on geodesic directions.")
print("      0.00° alignment = all alive cells are depth-1 neighbors of center.")
print("      Metric is most meaningful when population spans depth 2+")
