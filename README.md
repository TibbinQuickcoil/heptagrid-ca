# heptagrid-ca

**Transient Geodesic Alignment in Life-Like Cellular Automata on the {7,3} Heptagrid**

Clinton Foster — Independent Researcher

> *Preprint in preparation for arXiv:nlin.CG*

---

## Overview

This repository contains all code, data, and figures accompanying the paper *Transient Geodesic Alignment in Life-Like Cellular Automata on the {7,3} Heptagrid: Empirical Findings and a Structural Mechanism.*

We investigate whether the curvature of the hyperbolic plane leaves measurable fingerprints on the dynamics of Life-like cellular automaton rules — even when those rules are themselves geometrically blind. Working on the {7,3} heptagrid (the regular tiling of the hyperbolic plane by heptagons meeting three to a vertex) in the Poincaré disk model, we introduce a geodesic alignment metric and survey eight Life-like birth-survival rule families across multiple random initial conditions at two grid scales (depth 4: 190 cells; depth 5: 519 cells).

**Primary finding:** Under rule B2/S23, all twenty tested random initial conditions at both grid sizes produce minimum alignment scores below the 12.857° random baseline. The effect is transient — it appears during early growth and relaxes toward baseline as grid occupancy increases — and is explained by a two-component structural mechanism derived directly from the {7,3} neighbor graph.

**Retraction:** An initially observed period-4 oscillator under B34/S234 does not persist at depth 5 and is a finite-grid artifact.

---

## Repository Structure

```
heptagrid-ca/
│
├── README.md                      — this file
│
├── explorer/
│   └── poincare_ca_analysis.html  — interactive browser-based CA explorer
│                                    with live geodesic alignment metric,
│                                    history chart, and rule survey
│
├── scripts/
│   ├── survey_depth4.py           — depth-4 rule survey (190 cells, 20 seeds)
│   ├── survey_depth5.py           — depth-5 validation (519 cells, 20 seeds)
│   ├── mechanism.py               — structural mechanism computation
│   │                                (geodesic neighborhood concentration,
│   │                                shared neighbor analysis)
│   └── trajectories.py            — alignment trajectory data for figures
│
├── data/
│   ├── survey_depth4.json         — depth-4 survey results (all rules)
│   ├── survey_depth5_b2s23.json   — depth-5 B2/S23 results (20 seeds)
│   └── mechanism_tables.json      — Tables 2 and 3 from the paper
│
├── figures/
│   ├── fig1_geodesic.svg          — geodesic classification of the depth-4 grid
│   └── fig2_trajectory.svg        — alignment trajectories under B2/S23
│
└── paper/
    └── heptagrid_paper_v3.md      — current draft manuscript
```

---

## Quick Start

### Interactive Explorer

Open `explorer/poincare_ca_analysis.html` in any modern browser — no server required, no dependencies. The file is self-contained.

Features:
- Real-time geodesic alignment score with visual bar and verdict
- Alignment history chart with baseline and threshold lines
- Full rule survey across all 8 rule families (generates the Table 1 data)
- CSV and JSON export of survey results
- Click any cell to toggle its state

### Python Scripts

Requirements: Python 3.8+, standard library only (no numpy, no scipy).

```bash
# Run the depth-4 survey
python scripts/survey_depth4.py

# Run the depth-5 validation
python scripts/survey_depth5.py

# Compute the structural mechanism (Tables 2 and 3)
python scripts/mechanism.py

# Generate figure data
python scripts/trajectories.py
```

All scripts are self-contained and will print results to stdout.

---

## Key Results

### Geodesic Alignment Metric

For each alive cell at Poincaré disk position (c_r, c_i), angular distance to the nearest of 7 primary geodesic directions:

```
δ(θ) = min(t, 2π/7 - t)   where t = θ mod (2π/7)
```

Alignment score A = mean of δ across all non-center alive cells, in degrees.
Random baseline: E[δ] = π/14 ≈ **12.857°**

### B2/S23 Survey Results (depth 4, 20 seeds, 100 generations)

| Metric | Value |
|---|---|
| Mean minimum alignment | 7.998° |
| Reduction from baseline | 37.8% |
| Seeds below baseline | 20/20 |
| Seeds below 7.071° threshold | 2/20 |

### Depth-5 Validation (519 cells, 20 seeds, 150 generations)

| Metric | Value |
|---|---|
| Mean minimum alignment | 9.516° |
| Reduction from baseline | 26.0% |
| Seeds below baseline | 20/20 |

### Structural Mechanism (depth-4 neighbor graph)

| Depth | On-geodesic mean geo-neighbors | Off-geodesic mean | Ratio |
|---|---|---|---|
| 3 | 3.000 | 0.400 | 7.5× |
| 4 | 1.667 | 0.167 | 10.0× |

Off-geodesic adjacent pairs share 30% more neighbors than on-geodesic pairs (1.043 vs 0.800). Both properties are rule-independent and intensify at depth 5.

---

## Tiling Construction

The {7,3} tiling is constructed by BFS from the origin in the Poincaré disk using two geometric constants derived from the Schläfli symbol:

```
R_c = arccosh(cos(π/3) / sin(π/7))  ≈ 0.574   # circumradius (center to vertex)
r   = arccosh(cos(π/7) / sin(π/3))  ≈ 0.284   # inradius (center to edge midpoint)

ρ_cc = tanh(R_c)                               # cell-to-cell Euclidean displacement
ρ_v  = tanh(D/2)  where D = arccosh(cosh(R_c)·cosh(r))  # vertex placement
```

Cell positions are verified to resolve exactly 7 neighbors per interior cell.

---

## Reproducing the Paper Results

To reproduce Table 1 (rule survey):
```bash
python scripts/survey_depth4.py
```

To reproduce Table 4 (depth-5 validation):
```bash
python scripts/survey_depth5.py
```

To reproduce Tables 2 and 3 (structural mechanism):
```bash
python scripts/mechanism.py
```

Raw output matches the paper values exactly. All random seeds are fixed and specified in each script.

---

## Notes on the Retracted Finding

An earlier version of this work reported a period-4 oscillator under B34/S234, with alive-cell count cycling (39, 37, 40, 36) within the 190-cell depth-4 grid. This finding does not persist at depth 5. The oscillator is real within the 190-cell grid; it is not a property of the rule on the infinite {7,3} tiling.

The depth-4 oscillator data is preserved in `data/survey_depth4.json` for completeness.

---

## Citation

If you use this code or data, please cite the accompanying preprint (arXiv link forthcoming):

```
Foster, C. (2025). Transient geodesic alignment in Life-like cellular automata
on the {7,3} heptagrid: Empirical findings and a structural mechanism.
arXiv preprint. https://github.com/TibbinQuickcoil/heptagrid-ca
```

---

## Acknowledgments

The interactive explorer and Python scripts were developed with AI-assisted coding tools (Anthropic Claude). The analysis, experimental design, multi-scale validation, and decision to retract the oscillator finding were developed through human-AI collaboration and verified computationally. Peer review by Gemini 2.5 and Kimi contributed to the revision.

---

## License

MIT License. See LICENSE file.

The paper draft (`paper/`) is licensed under CC BY 4.0.
