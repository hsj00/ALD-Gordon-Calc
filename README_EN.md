# ALD Gordon Model Calculator v5 — User Manual

> Gordon et al. *Chem. Vap. Deposition* 9, 73 (2003) · Cremers et al. *Appl. Phys. Rev.* 6, 021302 (2019)

---

## Table of Contents

1. [Introduction](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#1-introduction)
2. [Installation &amp; Launch](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#2-installation--launch)
3. [Core Physical Model](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#3-core-physical-model)
4. [Input Parameter Selection Guide](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#4-input-parameter-selection-guide)
5. [Practical Examples](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#5-practical-examples)
6. [Interpreting Results](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#6-interpreting-results)
7. [Using the Graph Tabs](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#7-using-the-graph-tabs)
8. [Fill Tank Model Details](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#8-fill-tank-model-details)
9. [Running Unit Tests](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#9-running-unit-tests)
10. [Limitations &amp; Caveats](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#10-limitations--caveats)
11. [References](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#11-references)

---

## 1. Introduction

This calculator predicts the  **minimum precursor exposure (in Langmuirs) and dose time required to achieve conformal ALD coating in high-aspect-ratio (HAR) structures** .

 **Problem it solves** : "For an AR=50:1 cylindrical hole coated with TMA/H₂O Al₂O₃, how many Langmuirs of exposure are needed to reach 100% step coverage at the bottom?"

 **Key capabilities** :

* Required exposure and pulse time via the Gordon model
* Actual exposure from Fill Tank ODE integration
* Time-resolved penetration depth tracking
* Coverage profile θ(z) visualization
* Structure geometry comparison (Hole, Trench, Pillar)

 **Target users** : ALD process engineers, semiconductor device integration engineers, thin film researchers

---

## 2. Installation & Launch

### Prerequisites

```bash
pip install streamlit numpy pandas plotly scipy
```

### Running the App

```bash
streamlit run ald_gordon_calculator_v5.py
```

Open your browser at `http://localhost:8501`.

### Unit Tests

```bash
python run_tests.py
```

Or with pytest:

```bash
pytest test_gordon_calculator.py -v
```

---

## 3. Core Physical Model

### 3.1 Gordon Required Exposure (Eq. 14)

The minimum precursor exposure for 100% conformal coating to the bottom of a structure:

```
E_required = K_max × √(2π·m·k_B·T) × (1/19 + 4a/3 + 2a²)
```

where `a` is the Equivalent Aspect Ratio (EAR).

### 3.2 Equivalent Aspect Ratio (EAR)

Structures with the same geometric AR can have very different coating difficulty depending on shape:

| Structure           | EAR Formula      | Example (L=5 μm, w=100 nm) |
| ------------------- | ---------------- | --------------------------- |
| Cylindrical Hole    | a = L/w          | 50                          |
| Square Hole         | a = L/w          | 50                          |
| Infinite Trench     | a = L/(2w)       | 25                          |
| Elongated Hole      | a = L(w+z)/(2wz) | Geometry dependent          |
| Square Pillar Array | a = L/(2√2·w)  | 17.7                        |

 **Key insight** : A trench is open on both sides, so EAR = AR/2. The same AR is twice as easy to coat in a trench as in a hole.

### 3.3 Penetration Depth (Eq. 24)

The depth to which the precursor has penetrated at a given time:

```
l(t) = (4w/3) × [√(1 + (3/8) × E(t)/scale) − 1]
```

* Early stage (small E): l ∝ E (linear growth)
* Late stage (large E): l ∝ √E (diffusion-limited)
* Doubling penetration depth requires 4× the exposure

### 3.4 Knudsen Number

```
Kn = λ/w,   λ = k_B·T / (√2·π·d²·P)
```

| Kn Range      | Flow Regime    | Physical Picture                                | Gordon Model   |
| ------------- | -------------- | ----------------------------------------------- | -------------- |
| Kn > 10       | Molecular Flow | Molecules bounce off walls like ping-pong balls | Accurate       |
| 0.1 < Kn < 10 | Transition     | Both wall and intermolecular collisions         | Approximate    |
| Kn < 0.1      | Viscous Flow   | Fluid-like behavior                             | Not applicable |

---

## 4. Input Parameter Selection Guide

### 4.1 Structure Parameters: Depth (L) and Width (w)

 **Where to get these values** :

* Design rules for target L and w
* Actual etch results from TEM or SEM cross-section analysis
* CD-SEM for top/bottom CD measurements (use bottom CD for w — most conservative)

 **Representative values** :

| Structure                  | Typical Depth (L) | Typical Width (w) | AR         |
| -------------------------- | ----------------- | ----------------- | ---------- |
| DRAM Capacitor (DDR4)      | 3–5 μm          | 60–100 nm        | 40–80:1   |
| DRAM Capacitor (DDR5+)     | 5–8 μm          | 40–70 nm         | 80–150:1  |
| 3D NAND 96L Channel Hole   | 6–7 μm          | 100–140 nm       | 50–70:1   |
| 3D NAND 128L Channel Hole  | 8–9 μm          | 100–120 nm       | 70–90:1   |
| 3D NAND 200L+ Channel Hole | 12–15 μm        | 90–110 nm        | 120–160:1 |
| 3D NAND Slit (Trench)      | 6–10 μm         | 100–200 nm       | 30–100:1  |
| FinFET Gate Trench         | 40–60 nm         | 7–20 nm          | 2–8:1     |
| GAA Nanosheet Gap          | 8–12 nm          | 8–12 nm          | ~1:1       |

 **Practical tip** : Use the narrowest point (bottleneck) in the structure for w. For tapered etch profiles, the minimum CD gives the most conservative (safest) result.

### 4.2 Sticking Coefficient (s₀)

 **Physical meaning** : The probability that a precursor molecule reacts upon hitting the surface. Ranges from 0 to 1.

 **How to measure** :

1. **Lateral HARS test structures** : Deposit ALD on structures with known AR, measure step coverage profile by TEM, then back-calculate s₀ using the Gordon model
2. **QCM (Quartz Crystal Microbalance)** : Directly measure surface reaction rate in-situ
3. **Monte Carlo simulation fitting** : Fit simulation to experimental coverage profiles

**Literature values** (from Gordon 2003, Cremers 2019):

| Precursor            | Film          | s₀        | Source                  |
| -------------------- | ------------- | ---------- | ----------------------- |
| TMA [Al(CH₃)₃]     | Al₂O₃       | ~0.01      | Cremers (2019)          |
| TiCl₄               | TiO₂/TiN     | ~0.006     | Cremers (2019)          |
| TEMAHf [Hf(NEtMe)₄] | HfO₂         | ~0.1       | Gordon (2003)           |
| DEZ [Zn(C₂H₅)₂]   | ZnO           | ~0.007     | Gordon (2003)           |
| BDEAS                | SiO₂ (PEALD) | ~3×10⁻⁵ | Literature estimate     |
| MoCl₅               | Mo            | ~0.05      | Limited data — measure |
| MoO₂Cl₂            | Mo            | ~0.04      | Limited data — measure |

 **Caution** : s₀ is highly dependent on process temperature, surface termination, and co-reactant type. For accurate process optimization, experimental measurement is essential. Literature values should only be used for initial design estimates.

### 4.3 GPC (Growth Per Cycle)

 **Where to get it** :

* Measure on a **flat substrate** under saturation conditions (within the ALD window)
* Methods: Ellipsometry (most common), XRR (X-ray Reflectometry), TEM cross-section
* Always verify by measuring the GPC saturation curve (GPC vs. dose time) and use the saturated value

 **Representative values** :

| Film                              | Typical GPC (nm/cycle) |
| --------------------------------- | ---------------------- |
| Al₂O₃ (TMA+H₂O)                | 0.08–0.12             |
| TiO₂ (TiCl₄+H₂O)               | 0.04–0.06             |
| HfO₂ (TEMAHf+H₂O)               | 0.08–0.12             |
| ZnO (DEZ+H₂O)                    | 0.15–0.22             |
| SiO₂ (BDEAS+O₂ plasma)          | 0.10–0.15             |
| Mo (MoCl₅+H₂, or MoO₂Cl₂+H₂) | 0.03–0.06             |

### 4.4 K_max (Maximum Surface Adsorption Density)

**Automatic estimation** (Film DB mode):

```
K_max = GPC [m] × ρ_film [kg/m³] × N_A / MW_film [g/mol]
```

 **Uncertainty** : ±20–50%. The actual K_max depends on surface –OH density, steric hindrance of ligands, and deposition temperature.

 **More accurate method** : In-situ QCM measurement of mass gain per cycle (Δm), then compute K_max = Δm × N_A / MW_film directly.

### 4.5 Process Temperature (T) and Reference Pressure (P)

 **Temperature** :

* Enter the deposition temperature within the ALD window
* In the model, temperature affects only molecular thermal velocity (∝ √T); reaction activation energy is not considered
* Typical range: 150–400°C

 **Reference pressure** :

* Used in "Constant pressure P×t" mode as the precursor partial pressure during dosing
* Read from the chamber manometer during dosing or from the process recipe
* Units: mTorr (1 mTorr = 0.133 Pa)
* Typical range: 10–1000 mTorr

### 4.6 Fill Tank Parameters

Parameters needed for Fill Tank mode and how to determine each:

#### Fill Tank Volume (V_fill)

* **Where to find** : Equipment manual — look for "canister volume", "precursor cylinder headspace", or "manifold volume"
* **Important** : Include dead volume in the plumbing between the valve and chamber
* **Typical range** : 10–200 cc
* **Measurement method** : Fill the tank with known-pressure gas, open the valve to expand into a known volume, and back-calculate using PV = nRT

#### Pre/Post-fill Pressures (P_before, P_after)

* **P_before** : Pressure in the fill tank after filling with precursor vapor. Determined by the precursor vapor pressure at canister temperature
* **P_after** : Residual pressure in the fill tank after releasing precursor to the chamber
* **Measurement** : Direct reading from the gauge (e.g., Baratron) on the fill tank
* **Example** : TMA at 25°C has a vapor pressure of ~11 Torr. P_before = 10 Torr, P_after = 5 Torr gives ΔP = 5 Torr

#### Pump Speed (S)

* **Caution** : Pump spec sheet values are  **no-load pumping speeds** . Effective S at the chamber is reduced to **50–80%** due to plumbing conductance
* **Measurement** : Stabilize chamber at a steady P → shut off gas supply → record pressure decay curve → S_eff = V_c × (dP/dt) / P
* **Typical range** : 50–500 L/s (turbo pump), 10–50 L/s (dry pump)

#### Valve Conductance (C)

* **When you can skip it (Fast-fill condition)** : ALD valve is fully open and C >> S. The fill tank equilibrates with the chamber nearly instantly
* **When you need it** : Valve is partially open, flow restrictor is present, or the fill tank pressure is observed to change slowly during dosing
* **Decision criterion** : If fill tank pressure drops to P_after within tens of ms, it is fast-fill — C input is unnecessary
* **Typical values** : Swagelok ALD valve ~20–50 L/s (N₂ basis)

---

## 5. Practical Examples

### Example 1: DRAM Capacitor HfO₂ Coating

 **Scenario** : ALD HfO₂ dielectric coating on a DRAM capacitor. Cylindrical hole, depth 4 μm, width 80 nm.

 **Step 1 — Use a preset** :

* Select "DRAM Capacitor (typical)" from the sidebar
* Values auto-fill: L=4 μm, w=80 nm, TEMAHf, HfO₂, T=250°C, P=100 mTorr

 **Step 2 — Adjust to actual structure** :

* TEM shows actual bottom CD = 75 nm → change w to 75
* Actual etch depth = 4.2 μm → change L to 4.2

 **Step 3 — Read results** :

* AR = 56:1, EAR = 56:1 (Cylindrical Hole, so EAR = AR)
* Required exposure ≈ 350 L
* Required pulse time at 100 mTorr ≈ 0.46 s
* Kn = 8.2 → Molecular Flow → Gordon model valid

 **Step 4 — Verify with Fill Tank** :

* Select Fill Tank mode
* V_fill=50 cc, P_before=10 Torr, P_after=5 Torr, V_chamber=10 L
* t_dose=0.5 s, S_pump=100 L/s
* Saturation = 85% → Insufficient! → Increase t_dose to 0.7 s → Saturation = 98%

### Example 2: 3D NAND Al₂O₃ Blocking Oxide

 **Scenario** : Al₂O₃ blocking oxide in a 128-layer 3D NAND channel hole. Cylindrical hole, depth 8 μm, width 110 nm.

 **Step 1** : Select "3D NAND Channel Hole" preset and adjust values

 **Step 2** : Set:

* L=8 μm, w=110 nm, TMA, Al₂O₃
* s₀=0.01 (Cremers 2019), GPC=0.10 nm/cycle
* T=300°C, P=200 mTorr

 **Step 3 — Read results** :

* EAR = 72.7:1
* Required exposure ≈ 1,200 L → 🔴 Extremely challenging
* Kn = 3.8 → Transition flow → Approximate, correction may be needed

 **Step 4 — Optimize Fill Tank** :

* Initial conditions give saturation = 55% → Far insufficient
* Strategy: Increase ΔP by raising fill tank temperature (P_before ↑) or increase V_fill
* Set P_before=20 Torr, V_fill=100 cc → Saturation = 110% → Sufficient

### Example 3: Flash Memory Mo Word Line

 **Scenario** : Mo metal ALD for 3D NAND Flash word line using MoO₂Cl₂ + H₂.

 **Step 1** : Select "Flash WL Mo Fill" preset

 **Step 2** : Settings:

* L=6 μm, w=50 nm → AR=120:1
* MoO₂Cl₂ (MW=198.85), s₀=0.04, GPC=0.04
* T=500°C, P=300 mTorr

 **Step 3** :

* EAR = 120:1 → 🔴 Extreme structure
* Required exposure ≈ 15,000 L
* Multi-pulse dosing strategy is essential at this level

---

## 6. Interpreting Results

### 6.1 Reading the Result Cards

Five metric cards appear at the top after calculation:

| Card              | Meaning                                 | Favorable Direction          |
| ----------------- | --------------------------------------- | ---------------------------- |
| AR                | Geometric aspect ratio                  | Lower = easier to coat       |
| EAR               | Equivalent aspect ratio (model input)   | Lower = easier to coat       |
| Required Exposure | Minimum exposure for 100% step coverage | Lower = easier process       |
| Required Pulse    | Minimum dose time at reference pressure | Shorter = higher throughput  |
| Knudsen Kn        | Flow regime indicator                   | > 10 means model is reliable |

### 6.2 Traffic Light Indicators

Each result shows a 🟢🟡🔴 indicator with a one-line interpretation:

 **EAR interpretation** :

* 🟢 EAR < 10: Easy to coat with most ALD processes
* 🟡 EAR 10–50: Moderate difficulty, ensure sufficient exposure (DRAM class)
* 🔴 EAR 50–100: Very challenging (latest 3D NAND class)
* 🔴 EAR > 100: Extreme, consider precursor redesign or multi-pulse

 **Required exposure interpretation** :

* 🟢 < 100 L: Standard for typical ALD processes
* 🟡 100–1,000 L: Challenging; need extended dose time or higher pressure
* 🔴 > 1,000 L: Extremely difficult; multi-pulse or process strategy change needed

 **Knudsen number interpretation** :

* 🟢 Kn > 10: Molecule-wall collisions dominate. Molecules bounce off walls like ping-pong balls. Gordon model is accurate
* 🟡 Kn 0.1–10: Transition regime. Approximate application
* 🔴 Kn < 0.1: Viscous flow. Fluid-like behavior. Navier-Stokes-based models needed

### 6.3 Saturation Ratio Interpretation

The ratio of actual exposure to required exposure indicates coating completeness:

| Saturation | Meaning                    | Coverage Profile                 | Action                                                |
| ---------- | -------------------------- | -------------------------------- | ----------------------------------------------------- |
| ≥ 100%    | Sufficient                 | Uniform GPC to bottom            | Maintain conditions (excess wastes time)              |
| 90–99%    | Nearly achieved            | Slight thinning at bottom 5–10% | Slightly increase dose time                           |
| 70–89%    | Insufficient               | Significant bare area at bottom  | Increase dose time or pressure by 1.2–1.5×          |
| 50–69%    | Significantly insufficient | > half of structure uncoated     | Fundamental condition change needed                   |
| < 50%      | Critical                   | No meaningful coating            | Redesign fill tank capacity, pump, precursor strategy |

### 6.4 Fill Tank Model A vs. B Difference

| Model              | Definition                                        | Meaning                          |
| ------------------ | ------------------------------------------------- | -------------------------------- |
| A: P_eq × t_dose  | Exposure upper bound (assumes constant P_eq)      | Optimistic estimate              |
| B: ODE integration | Actual exposure accounting for pressure transient | Realistic estimate (recommended) |

 **When the two differ significantly (> 20%)** :

* If t_dose < τ: Dose time is short; fill tank precursor is not fully utilized → **Extend dose time to 2–3× τ** for effective exposure increase
* If t_dose > τ: Pumping removes precursor quickly → **Increase ΔP or V_fill** to raise E_max itself

 **When the two are similar (< 20%)** :

* Fast-fill approximation is valid. Current valve/plumbing design is adequate

### 6.5 Time-Resolved Penetration Depth

 **l(t) graph elements** :

* Blue solid line: Penetration depth vs. time
* Red dashed line: Total structure depth L
* Red dot: Penetration depth at current t_dose
* Green dashed line: Time to achieve 100% coating (t_full)

 **Key metrics** :

* **Penetration depth @ t_dose** : How deep the precursor has reached at the current dose time. Compare with L for completeness assessment
* **t_full** : Minimum time to coat to the bottom. If t_dose > t_full, coating is complete
* **l_max (t→∞)** : Maximum achievable penetration depth from a single fill tank release. If l_max < L, 100% coating is impossible with current fill tank conditions regardless of dose time

---

## 7. Using the Graph Tabs

### Tab 1: Exposure vs. EAR

 **Question** : "If the structure shrinks and EAR increases, how much more exposure is needed?"

* The red dashed line (current EAR) shows the required exposure
* At high EAR, exposure scales as EAR² (quadratic) — confirmed by the 2a² asymptote
* Use this to forecast exposure requirements for next-generation structures

### Tab 2: Pulse Time vs. Pressure

 **Question** : "At what pressure should I dose to achieve the target dose time?"

* Visualizes the P × t = const relationship
* Orange region (Kn < 1) is where Gordon model breaks down — avoid operating there
* Combine with equipment pumping limits to find the optimal operating point

### Tab 3: Pulse Time vs. Temperature

 **Question** : "How does changing the process temperature affect dose time?"

* T ↑ → thermal velocity ↑ → slightly longer pulse needed for same exposure (weak dependence)
* Lower graph: λ vs. T → check if Kn improves
* Combine with ALD window to determine optimal temperature

### Tab 4: EAR Structure Comparison

 **Question** : "For the same AR, is a trench or hole easier to coat?"

* Hole: EAR = AR (most difficult)
* Trench: EAR = AR/2 (open on both sides)
* Pillar: EAR = AR/(2√2) (easiest)
* Quantitatively compare the effect of changing structure geometry

---

## 8. Fill Tank Model Details

### 8.1 Governing Equations

In fill-tank-based ALD, precursor expands from a high-pressure fill tank into a low-pressure reaction chamber:

```
dP_f/dt = −C(P_f − P_c) / V_f        (fill tank pressure decrease)
dP_c/dt =  C(P_f − P_c) / V_c − S·P_c / V_c   (chamber pressure change)
```

### 8.2 Key Parameters

* **P_eq = ΔP × V_fill / V_chamber** : Equilibrium chamber pressure after instantaneous mixing
* **τ = V_chamber / S_pump** : Time constant for precursor pump-out from the chamber
* **E_max = P_eq × τ** : Maximum achievable exposure from a single fill tank release

### 8.3 Practical Decision Criteria

| Condition          | Meaning                            | Action                                   |
| ------------------ | ---------------------------------- | ---------------------------------------- |
| E_max > E_required | 100% coating is possible           | Set t_dose ≥ t_full                     |
| E_max < E_required | 100% coating impossible            | Increase ΔP, V_fill, or decrease S_pump |
| t_dose/τ < 0.5    | Dose time too short                | Extending t_dose is effective            |
| t_dose/τ > 3      | Diminishing returns from more time | Increase ΔP or V_fill instead           |

---

## 9. Running Unit Tests

```bash
# Standalone runner (no streamlit/plotly needed)
python run_tests.py

# Or with pytest
pytest test_gordon_calculator.py -v
```

Test coverage (44 tests total):

* Unit conversion accuracy
* Gordon EAR calculation verification
* K_max estimation formula validation
* Mean free path calculation (including boundary conditions)
* Fill Tank ODE vs. analytical solution agreement (5 dose times, relative error < 2%)
* Model A ≥ Model B inequality
* Cumulative exposure monotonicity and E_max convergence
* Penetration depth √E asymptotic behavior
* find_t_full forward/reverse verification
* Input validation (physically impossible input defense)

---

## 10. Limitations & Caveats

### Model Limitations

1. **Assumes s₀ = 1** : The Gordon model provides results for the fully diffusion-limited case. For actual s₀ < 1, the coverage profile is smoother, so the model provides a **conservative upper bound** for the required exposure
2. **First-order Langmuir adsorption** : Only considers s(θ) = s₀(1−θ). Multi-layer adsorption and complex reaction mechanisms are not included
3. **Isothermal assumption** : Temperature gradients within the structure are not considered
4. **Ideal gas assumption** : Valid at low pressures (< a few Torr). Real gas corrections needed at higher pressures

### Practical Caveats

1. **Kn < 1 regime** : Gordon model is inaccurate in viscous flow. Use CFD (computational fluid dynamics) simulations for these conditions
2. **Fill tank pressure measurement** : P_before and P_after are gauge readings at the fill tank; pressure drops in the plumbing must be considered separately
3. **Precursor decomposition** : Some precursors (especially metalorganics) can undergo thermal decomposition at high temperatures, resulting in a CVD component. The Gordon model considers pure ALD only
4. **Surface chemistry variations** : The surface chemistry inside the structure may vary with depth (e.g., SiO₂ at the top, Si₃N₄ at the bottom). In such cases s₀ varies with position, reducing model accuracy

---

## 11. References

1. R. G. Gordon, D. Hausmann, E. Kim, J. Shepard, "A Kinetic Model for Step Coverage by Atomic Layer Deposition in Narrow Holes or Trenches," *Chem. Vap. Deposition*  **9** , 73–78 (2003)
2. V. Cremers, R. L. Puurunen, J. Dendooven, "Conformality in Atomic Layer Deposition: Current Status Overview of Analysis and Modelling," *Appl. Phys. Rev.*  **6** , 021302 (2019)
3. S. M. Sze, M. K. Lee,  *Semiconductor Devices: Physics and Technology* , 3rd Ed., Wiley (2012)
4. B. A. Neaman,  *Semiconductor Physics and Devices: Basic Principles* , 4th Ed., McGraw-Hill (2012)

---

*ALD Gordon Model Calculator v5 — 2025*
