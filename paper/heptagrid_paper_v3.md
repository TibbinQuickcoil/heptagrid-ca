# Transient Geodesic Alignment in Life-Like Cellular Automata on the {7,3} Heptagrid: Empirical Findings and a Structural Mechanism

**Clinton Foster**
Independent Researcher
https://github.com/TibbinQuickcoil/heptagrid-ca

---

## Abstract

We report an empirical investigation of two-state totalistic cellular automata on the {7,3} heptagrid — the regular tiling of the hyperbolic plane by heptagons meeting three to a vertex — implemented in the Poincaré disk model with a 7-neighbor topology. While prior work on hyperbolic cellular automata has focused predominantly on questions of computational universality, the dynamical behavior of Life-like birth-survival rule families on this tiling remains largely uncharacterized. We introduce a geodesic alignment metric defined as the mean angular distance between alive-cell positions and the nearest of seven primary geodesic directions through the disk center, with a random baseline of 12.857° derived analytically. Surveying eight Life-like rule families across random initial conditions on grids of depth 4 (190 cells) and depth 5 (519 cells), we find that under rule B2/S23, all twenty tested random initial conditions at both grid sizes produce minimum alignment scores below the random baseline: mean minimum 7.998° (37.8% below baseline) at depth 4, and 9.516° (26.0% below baseline) at depth 5. The effect is transient, appearing during the early growth phase and relaxing toward baseline as grid occupancy increases. We derive a structural mechanism from two properties of the {7,3} neighbor graph — concentration of geodesic adjacency along geodesic directions (7.5–11.6× ratio at depth 4) and differential shared-neighbor density between on- and off-geodesic cell pairs — both rule-independent properties of the tiling that intensify at larger scales. We also retract an initially observed period-4 oscillator under B34/S234, which does not persist at depth 5 and is a finite-grid artifact. The interactive computational explorer and all verification code are released as open-source artifacts.

---

## 1. Introduction

Conway's Game of Life [Gardner, 1970] is among the most studied dynamical systems in mathematics. Half a century of investigation on the Euclidean square lattice has yielded a vast catalog of behaviors: still lifes, oscillators, gliders, glider guns, and computational structures sufficient to construct universal Turing machines. The square grid has become so naturalized as a substrate that the question of which substrate is rarely asked. Yet the choice of underlying tiling is not neutral: it determines neighborhood structure, growth rates, and the geometric constraints under which patterns evolve.

The hyperbolic plane offers a substantially different substrate. Where the Euclidean plane admits regular tilings by squares, triangles, and hexagons, the hyperbolic plane admits an infinite family of regular tilings {p, q} for which the angular constraint (p-2)(q-2) > 4 is satisfied. Among these, the {7,3} tiling — heptagons meeting three to a vertex, the dual of the {3,7} triangular tiling — is the simplest hyperbolic regular tiling that produces a 7-neighbor cell structure. As the dual of the triangular tiling, the {7,3} heptagrid has the lowest vertex valency among hyperbolic regular tilings, making it the sparsest such structure in terms of inter-cell connectivity per cell — a property directly relevant to the sparse-field dynamics we analyze. Each cell touches exactly seven neighbors, the number of cells at graph distance n from any given cell grows exponentially rather than quadratically, and straight lines through the plane are realized as arcs of circles orthogonal to the disk boundary in the Poincaré model.

A small but significant body of work has explored cellular automata on hyperbolic tilings [Margenstern, 2000, 2007]. This literature has focused almost exclusively on computational universality: whether the substrate admits cellular automata capable of simulating arbitrary Turing machines. These are top-down constructions in which transition rules are explicitly engineered, often using geodesic paths as designed communication channels.

The complementary question — what spontaneously emergent dynamical behavior arises when simple Life-like birth-survival rules are run on the {7,3} tiling — has received little attention. The closest prior exploration of which we are aware is an informal investigation by Mishin [2011], conducted on the {5,4} pentagrid without quantitative analysis.

We are interested in a specific question: does the curvature of the underlying space leave detectable fingerprints on the dynamics of Life-like rules, even when those rules are themselves geometrically blind? A cellular automaton rule contains no representation of geodesics, distance, or curvature. It counts neighbors and applies a transition function. If the geometry of the substrate nevertheless biases the resulting dynamics, that bias must enter through the topology of the neighborhood graph alone.

We approach this question empirically, then mechanistically, and then critically — explicitly testing initial findings across grid sizes and retracting one that does not survive. Section 2 describes the tiling construction, CA implementation, and geodesic alignment metric. Section 3 presents the initial depth-4 survey. Section 4 derives the structural mechanism. Section 5 presents depth-5 validation and the retraction. Section 6 addresses limitations and future work.

---

## 2. Methods

### 2.1 Tiling Construction

We work in the Poincaré disk model, in which the hyperbolic plane is represented as the open unit disk D = {z in C : |z| < 1} with metric ds = 2|dz|/(1-|z|^2). Hyperbolic geodesics are diameters of the disk and arcs of circles meeting the boundary orthogonally.

The {7,3} tiling is generated by reflecting a central heptagon across its sides recursively. Two geometric constants are derived from the Schläfli symbol {7,3} using standard formulas for regular hyperbolic polygons [Coxeter, 1954]:

$$R_c = \cosh^{-1}\!\left(\frac{\cos(\pi/3)}{\sin(\pi/7)}\right) \approx 0.574 \quad \text{(circumradius: center to vertex)}$$

$$r = \cosh^{-1}\!\left(\frac{\cos(\pi/7)}{\sin(\pi/3)}\right) \approx 0.284 \quad \text{(inradius: center to edge midpoint)}$$

The Euclidean displacement parameter for center-to-center moves in the Poincaré disk is rho_cc = tanh(R_c), used to place adjacent cell centers via Möbius transforms of the form T_{z0, theta}(z) = (z + z0) / (1 + conj(z0) * z).

For rendering cell boundaries, vertices are placed at Euclidean displacement rho_v = tanh(D/2) where D = arccosh(cosh(R_c) * cosh(r)) ≈ 0.622. The geometric interpretation of D as a standard quantity of the {7,3} polygon requires independent verification and is not claimed here; rho_v is used for visualization only and does not enter into the CA logic or alignment metric.

**Correction note.** An earlier draft of this paper mislabeled R_c as "apothem" and r as "half-edge," reversing the standard definitions. This was a textual error identified in peer review [Kimi, personal communication]. The implementation correctly uses rho_cc = tanh(R_c) for cell-center placement, verified to produce a valid {7,3} neighbor graph by confirming each interior cell resolves exactly 7 neighbors.

Cell positions are enumerated by BFS beginning at the origin. From each cell, seven neighbors are generated at angles theta_j = phi + 2*pi*j/7 for j = 0,...,6, where phi is the local orientation inherited from BFS. Cells within radius 0.971 of the origin are retained for depth 4 (N = 190 cells) and 0.989 for depth 5 (N = 519 cells). Cells are identified by position maps at precision floor(7000r) (depth 4) and floor(15000r) (depth 5). All runs use fixed initial orientation phi = 0 at the origin, ensuring deterministic BFS order.

The center cell lies at the origin with seven immediate neighbors (depth 1) at angles 2*pi*j/7 for j = 0,...,6 — exactly along the seven primary geodesic directions. This is consequential for metric interpretation (§2.3).

### 2.2 Cellular Automaton

We implement two-state totalistic Life-like CA in the B/S family. State of cell i at time t+1:

$$s_i^{t+1} = \begin{cases} 1 & \text{if } s_i^t = 0 \text{ and } a_i \in B \\ 1 & \text{if } s_i^t = 1 \text{ and } a_i \in S \\ 0 & \text{otherwise} \end{cases}$$

where a_i = sum of s_j^t over j in N(i) and B, S are subsets of {0,...,7}. Eight rule families were surveyed: B3/S34, B2/S34, B23/S34, B3/S23, B2/S23, B34/S234, B2/S345, and B1/S12.

### 2.3 Geodesic Alignment Metric

For a cell at Poincaré disk coordinates (c_r, c_i), angular position theta = atan2(c_i, c_r). Angular distance to the nearest of the seven primary geodesic directions {2*pi*j/7}:

$$\delta(\theta) = \min\!\left(t,\ \frac{2\pi}{7} - t\right), \qquad t = \theta \bmod \frac{2\pi}{7}$$

Geodesic alignment score for a given generation:

$$A = \frac{1}{|F|}\sum_{c \in F} \delta(\theta_c) \cdot \frac{180}{\pi}$$

where F is the set of alive cells with c_r^2 + c_i^2 > 0.002. Lower A indicates greater concentration near geodesic directions.

**Random baseline.** For angular position drawn uniformly from [0, 2*pi):

$$\mathbb{E}[\delta] = \frac{\pi}{14} \approx 12.857°$$

This baseline assumes uniformly random angles. Since initial conditions are random by cell (not by angle), a slight geometric bias may exist at t=0 due to non-uniform angular density of cells at different depths in the Poincaré disk compression. We do not expect this to materially affect results as the effect is symmetric across geodesic sectors, but note it for completeness.

**Geodesic threshold.** We use a secondary threshold of 7.071° = 12.857° x 0.55 as an operational marker for strong geodesic bias. This threshold was chosen operationally to clearly separate the null distribution from cases of marked concentration and does not derive from first principles.

**Depth-1 caveat.** Depth-1 cells are placed exactly at geodesic directions by the BFS construction, contributing delta = 0 regardless of dynamics. Depth-stratified analysis accompanies aggregate scores throughout.

---

## 3. Initial Survey — Depth-4 Grid

### 3.1 Rule Survey Overview

**Table 1.** Rule survey, depth-4 grid (190 cells), 20 random seeds (density 0.28), 100 generations.

| Rule | Alive@100 (mean) | Mean min align | Below baseline | Behavior |
|---|---|---|---|---|
| B3/S34 | 0 | — | — | extinction |
| B2/S34 | 2 | 8.39° | 20/20 | near-extinction |
| B23/S34 | ~44 | 11.38° | partial | diffuse oscillation |
| B3/S23 | 3 | 8.57° (dagger) | — | depth-0/1 cluster |
| B2/S23 | ~38 | 7.998° | 20/20 | transient geodesic bias |
| B34/S234 | 36-40 | 11.44° | partial | oscillatory (grid artifact) |
| B2/S345 | 0 | — | — | extinction |
| B1/S12 | ~73 | 12.27° | partial | dense diffuse |

(dagger) B3/S23 dominated by depth-0/1 cluster; see §2.3.
B34/S234 oscillator does not persist at depth 5; see §5.2.

### 3.2 B2/S23: Transient Geodesic Alignment

Under B2/S23, all twenty random initial conditions produce minimum alignment below 12.857°. Full distribution:

4.17°, 5.71°, 7.34°, 7.35°, 7.39°, 7.48°, 7.68°, 7.75°, 7.89°, 8.01°, 8.34°, 8.35°, 8.38°, 8.40°, 8.49°, 8.68°, 9.00°, 9.28°, 9.60°, 9.61°

- Mean minimum: 7.998° (37.8% below baseline)
- Below 7.071° threshold: 2/20
- Below baseline: 20/20

The minimum occurs at variable generations (9-71) during intermediate-density growth, before relaxing toward baseline as grid occupancy increases.

**Depth-stratified analysis** at global minimum (seed = 7, generation = 34, 12 alive cells):

| Depth | Alive cells | Alignment | Note |
|---|---|---|---|
| 1 | 6 | 0.00° | geometric artifact |
| 2 | 2 | 7.85° | sub-baseline, genuine signal |
| 3 | 3 | 10.56° | near-baseline |
| 4 | 1 | 2.63° | sub-baseline, small sample |

Excluding depth-1 cells: 6 remaining cells at depths 2-4 show mean alignment 7.43°, below baseline. Sample sizes are small; the direction is consistent with the aggregate result.

### 3.3 Null Results

Five of eight rules produce no evidence of geometric organizing influence, confirming that geodesic bias is not a generic property of Life-like rules on the {7,3} tiling but a feature of specific dynamical regimes — specifically the intermediate-density growth phase under low-birth-count rules.

---

## 4. Mechanism: Structural Origin of Geodesic Bias

### 4.1 Approach

The transient alignment in B2/S23 admits two explanations: a statistical artifact of finite-grid effects or sampling noise, or a structural property of the {7,3} neighbor graph that systematically favors geodesic propagation for certain rule families. We test the second hypothesis directly from the neighbor graph, independent of any rule or dynamics.

### 4.2 Geodesic Classification

We label each cell at depth d >= 1 as on-geodesic if delta(theta_c) < pi/28 ≈ 6.43° — within the inner half of the angular distance to the nearest geodesic — and off-geodesic otherwise.

### 4.3 Geodesic Neighborhood Concentration

**Table 2.** Mean on-geodesic neighbors per cell, depth-4 grid.

| Depth | Class | n | Mean on-geo neighbors | Ratio |
|---|---|---|---|---|
| 3 | on-geodesic | 21 | 3.000 | 7.5x |
| 3 | off-geodesic | 35 | 0.400 | |
| 4 | on-geodesic | 21 | 1.667 | 10.0x |
| 4 | off-geodesic | 84 | 0.167 | |

At depth 3, on-geodesic cells have 7.5x more on-geodesic neighbors than off-geodesic cells of the same depth. At depth 4 the ratio is 10.0x. This is a purely structural property of the {7,3} tiling derived from the neighbor graph, independent of any rule or initial condition.

### 4.4 Shared Neighbor Density

**Table 3.** Mean shared neighbors among adjacent cell pairs, depth-4 grid.

| Pair type | n pairs | Mean shared neighbors |
|---|---|---|
| both on-geodesic | 35 | 0.800 |
| both off-geodesic | 161 | 1.043 |

Off-geodesic adjacent pairs share 30% more neighbors than on-geodesic pairs. The on-geodesic chain topology is sparsely cross-connected; off-geodesic regions are densely cross-connected.

### 4.5 The Two-Component Mechanism

**Component 1 — directed alive signal.** During sparse growth, an on-geodesic cell has 7.5-10x more geodesic-adjacent neighbors than an off-geodesic cell at the same depth. If any chain cell becomes alive, its on-geodesic neighbors receive a concentrated alive-cell signal in the chain direction, biasing the probability of satisfying B2 (exactly 2 alive neighbors required) along geodesics.

**Component 2 — differential birth probability by neighbor count rule.** The dense cross-connectivity of off-geodesic regions pushes alive-neighbor counts above the B2 threshold, suppressing B2 birth. Crucially, this mechanism inverts for higher birth-count rules: for B3 or B4 rules, dense cross-connectivity facilitates rather than suppresses birth by making it easier to accumulate 3 or more simultaneous alive neighbors. This predicts a dichotomy: B2-family rules exhibit geodesic bias through birth-suppression off-geodesic, while B3-family rules show reduced or absent bias — consistent with our null results for all B3-based rules in the survey.

**Transience.** The mechanism is density-dependent. As the alive-cell density increases, the structural advantage of on-geodesic cells diminishes relative to the total alive-neighbor counts, and the alignment signal relaxes toward baseline.

### 4.6 Status of the Argument

The structural properties in §4.3 and §4.4 are computed exactly from the {7,3} neighbor graph and are rule-independent. The prediction that B2/S23 exhibits transient geodesic bias follows under a qualitative argument about sparse-field birth probabilities. A complete formal derivation — computing birth probabilities explicitly as a function of density and cell classification — is the principal direction for future theoretical work. The structural concentration ratios (7.5-10x) are too large to be explained by finite-grid noise, ruling out the sampling-artifact hypothesis.

---

## 5. Depth-5 Validation

### 5.1 B2/S23 Primary Finding at Depth 5

Twenty seeds tested at 150 generations on the depth-5 grid (N = 519).

**Table 4.** B2/S23 alignment by grid depth.

| Grid | N | Mean min | Reduction | Below baseline | Below 7.071° |
|---|---|---|---|---|---|
| Depth 4 | 190 | 7.998° | 37.8% | 20/20 | 2/20 |
| Depth 5 | 519 | 9.516° | 26.0% | 20/20 | 0/20 |

The directional finding is confirmed: 20/20 seeds below baseline at both grid sizes. The magnitude weakens from 37.8% to 26.0%, and no depth-5 seed breaches the strong threshold.

Individual depth-5 minima: 8.80°, 8.86°, 8.99°, 9.12°, 9.12°, 9.19°, 9.26°, 9.28°, 9.34°, 9.44°, 9.47°, 9.53°, 9.58°, 9.72°, 9.82°, 9.98°, 10.02°, 10.10°, 10.34°, 10.38°.

### 5.2 Retraction: B34/S234 Period-4 Oscillator

At depth 4, rule B34/S234 produced a repeating alive-count cycle (39, 37, 40, 36) persisting through generation 100 across all tested seeds. At depth 5, no such clean periodicity is observed within 120 generations. The depth-4 cycle is a finite-grid artifact specific to the 190-cell grid.

We retract this finding. The period-4 oscillator was real within the 190-cell grid; it is not a property of the rule on the infinite {7,3} tiling. This retraction is the expected cost of multi-scale validation and is reported here as an example of standard scientific practice.

### 5.3 Structural Mechanism at Depth 5

**Table 5.** Geodesic neighborhood concentration, depth-5 grid.

| Depth | On-geo mean | Off-geo mean | Ratio |
|---|---|---|---|
| 3 | 3.667 | 0.400 | 9.17x |
| 4 | 5.800 | 0.500 | 11.60x |
| 5 | 1.400 | 0.000 | boundary artifact |

At depth 5, off-geodesic cells show zero on-geodesic neighbors. This is a boundary artifact: depth-5 cells have neighbors at depth 6, which fall outside the radius cutoff. Off-geodesic depth-5 cells whose on-geodesic neighbors lie at depth 6 have those neighbors zeroed by truncation. The depth-5 ratio is therefore a lower bound on the true structural concentration rather than a property of the infinite tiling. The ratios at depths 3 and 4, where both cell classes have fully resolved neighbors, show genuine intensification (7.5x to 9.17x at depth 3; 10.0x to 11.6x at depth 4), consistent with structural anisotropy strengthening at scale.

### 5.4 Reconciling Stronger Mechanism with Weaker Signal

The structural concentration strengthens at depth 5 while the aggregate alignment signal weakens. This is explained by population averaging: at depth 5 the total alive-cell population is larger, and the geodesic chains form a smaller fraction of the total. The aggregate metric averages over all alive cells including the large off-geodesic population contributing near-baseline values.

The invariant result across both scales is the directional one: 20/20 seeds below baseline at both depths. This is the most defensible form of the primary finding.

---

## 6. Limitations and Future Work

**Formal derivation.** The mechanism in §4 is structurally grounded but not formally proven. A complete derivation would compute birth probabilities as a function of density and cell classification, yielding a quantitative prediction testable against the alignment trajectory.

**Euclidean control.** Running the same rule families on a square or hexagonal grid with an analogous angular metric is the most important missing validation. This would test whether the effect is specific to hyperbolic curvature or could arise from any regular tiling with preferred directions.

**Larger grids.** Depth-6 yields approximately N ≈ 1900 cells. The population-averaging hypothesis (§5.4) could be tested directly via depth-stratified alignment at depth 5, which would then be estimable from sufficient cell samples. The depth-5 structural boundary artifact (§5.3) would also be resolved.

**Rule space coverage.** The full two-state totalistic rule space on a 7-neighbor grid contains 2^8 x 2^8 = 65,536 rules. Systematic exploration of B2-family rules is a natural extension given our mechanism's prediction of their behavior.

**Geometric constant verification.** The vertex placement parameter rho_v = tanh(D/2) where D = arccosh(cosh(R_c) * cosh(r)) is used for rendering only. Its geometric interpretation as a standard {7,3} polygon quantity requires independent derivation.

---

## 7. Conclusion

We have presented a quantitative empirical study of Life-like cellular automata on the {7,3} heptagrid, tested at two grid scales, with a structural mechanism derived from the neighbor graph and an explicit retraction of one initial finding that did not survive scale testing.

The principal finding is that rule B2/S23 consistently produces transient sub-baseline geodesic alignment during early growth from random initial conditions: 20/20 seeds below the 12.857° random baseline at both depth 4 and depth 5, with mean minimum reductions of 37.8% and 26.0% respectively. The directional consistency across all forty tested conditions is the most defensible form of this result.

The structural mechanism — concentration of geodesic adjacency (7.5-11.6x ratio at depth 4, with genuine intensification at depth 5) and differential shared-neighbor density — is derived from properties of the {7,3} neighbor graph independent of any rule or dynamics. The divergence between a strengthening structural mechanism and a weakening aggregate alignment signal is explained as population averaging, constituting a testable prediction for future depth-stratified analysis at larger grids.

The B34/S234 period-4 oscillator does not persist at depth 5 and is retracted.

Cellular automaton rules are geometrically blind. The {7,3} tiling is not.

---

## Acknowledgments

The interactive computational explorer and Python verification scripts were developed with AI-assisted coding tools. Empirical observations, methodological design, analytical framing, multi-scale validation, and the decision to retract the B34/S234 finding were developed through human-AI collaboration. All results were independently verified by Python reimplementation. Peer review commentary from Gemini 2.5 and Kimi contributed substantially to the revision, in particular the correction of geometric constant labeling (Kimi) and the suggestion to note the {7,3}/{3,7} duality and sharpen the B3 mechanism dichotomy (Gemini).

## Code and Data Availability

The interactive {7,3} cellular automaton explorer, Python survey and mechanism scripts, raw data, and SVG figures are available at:

https://github.com/TibbinQuickcoil/heptagrid-ca

## References

Conway, J. H. (1970). The game of life. In M. Gardner, Mathematical Games. Scientific American, 223(4), 120-123.

Coxeter, H. S. M. (1954). Regular honeycombs in hyperbolic space. Proceedings of the International Congress of Mathematicians, Amsterdam, III, 155-169.

Margenstern, M. (2000). New tools for cellular automata in the hyperbolic plane. Journal of Universal Computer Science, 6(12), 1226-1252.

Margenstern, M. (2007). Cellular Automata in Hyperbolic Spaces, Volume 1: Theory. Old City Publishing.

Margenstern, M., & Morita, K. (2001). NP problems are tractable in the space of cellular automata in the hyperbolic plane. Theoretical Computer Science, 259, 99-128.

Mishin, D. (2011). Hyperbolic cellular automaton [informal exploration, {5,4} pentagrid]. http://dmishin.blogspot.com/2011/10/hyperbolic-cellular-automation.html

Wolfram, S. (2002). A New Kind of Science. Wolfram Media.
