
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
9. [Unit System &amp; Verification](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#9-unit-system--verification)
10. [Unit Tests](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#10-unit-tests)
11. [Limitations &amp; Caveats](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#11-limitations--caveats)
12. [References](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#12-references)

---

## 1. Introduction

This calculator predicts the  **minimum precursor exposure and dose time required for conformal ALD coating in high-aspect-ratio (HAR) structures** .

 **Key capabilities** :

* Required exposure and pulse time via the Gordon model
* Actual exposure from Fill Tank ODE integration
* Time-resolved penetration depth tracking
* Dose Multiple (HAR vs. flat substrate ratio)
* Coverage profile θ(z) visualization
* 5 presets (DRAM, 3D NAND, FinFET, Flash WL)
* Pressure unit selection: Torr or mTorr

---

## 2. Installation & Launch

```bash
pip install streamlit numpy pandas plotly scipy
streamlit run ald_gordon_calculator_v5.py
```

---

## 3. Core Physical Model

### 3.1 Gordon Required Exposure (Eq. 14)

```
E_required = K_max × √(2π·m·k_B·T) × (1/19 + 4a/3 + 2a²)
             ───────────────────────   ──────────────────────
                    scale [Pa·s]           f(a) [dimensionless]
```

Role of each term:

* **K_max** [molecules/m²]: molecules consumed per m² of wall — "how hungry is the surface"
* **√(2πmk_BT)** [kg·m/s]: momentum scale linking flux to exposure (P×t)
* **f(a)** : geometric diffusion resistance — "how deep must we push"

### 3.2 Equivalent Aspect Ratio (EAR)

| Structure           | EAR Formula   | Example (L=5 μm, w=100 nm) |
| ------------------- | ------------- | -------------------------- |
| Cylindrical Hole    | a = L/w       | 50                         |
| Infinite Trench     | a = L/(2w)    | 25                         |
| Square Pillar Array | a = L/(2√2·w) | 17.7                       |

A trench is open on both sides, so EAR = AR/2 — easier to coat than a hole at the same AR.

### 3.3 K_max — Film Property vs. Precursor Property

K_max is the  **maximum number of precursor molecules adsorbed per unit area per ALD half-cycle** .

```
K_max = GPC [m] × ρ_film [kg/m³] × N_A [/mol] / MW_film [g/mol]
```

Logic: "How many atoms fit in a GPC-thick film layer?" = "How many precursor molecules were consumed."

**Which property belongs to which?**

| Parameter              | Property of             | Why                          |
| ---------------------- | ----------------------- | ---------------------------- |
| MW (molecular weight)  | **Precursor**           | Gas-phase diffusion speed    |
| d (molecular diameter) | **Precursor**           | Knudsen number, MFP          |
| s₀ (sticking coeff.)   | **Precursor + surface** | Surface reaction probability |
| GPC                    | **Process measurement** | Growth per cycle             |
| ρ_film, MW_film        | **Film**                | Used for K_max estimation    |

Example: MoO₂Cl₂ depositing Mo metal

* Precursor MW = 198.85 g/mol (MoO₂Cl₂) → diffusion calculation
* Film ρ = 10.2 g/cm³, MW = 95.95 g/mol (Mo) → K_max calculation

### 3.4 Penetration Depth (Eq. 24)

```
l(t) = (4w/3) × [√(1 + (3/8) × E(t)/scale) − 1]
```

* Early: l ∝ E (linear)
* Late: l ∝ √E (diffusion-limited) — 4× exposure = 2× depth

### 3.5 Knudsen Number

| Kn Range | Flow Regime                     | Gordon Model   |
| -------- | ------------------------------- | -------------- |
| Kn > 10  | Molecular — ping-pong off walls | Accurate       |
| 0.1–10   | Transition — mixed collisions   | Approximate    |
| Kn < 0.1 | Viscous — fluid behavior        | Not applicable |

---

## 4. Input Parameter Selection Guide

### 4.1 Depth (L) and Width (w)

 **Measurement** : SEM/TEM cross-section, CD-SEM (use bottom CD)

| Structure                  | Depth (L)        | Width (w)              | AR        |
| -------------------------- | ---------------- | ---------------------- | --------- |
| DRAM Capacitor (DDR4)      | 3–5 μm           | 60–100 nm              | 40–80:1   |
| DRAM Capacitor (DDR5+)     | 5–8 μm           | 40–70 nm               | 80–150:1  |
| 3D NAND 128L Channel Hole  | 8–9 μm           | 100–120 nm             | 70–90:1   |
| 3D NAND 200L+ Channel Hole | 12–15 μm         | 90–110 nm              | 120–160:1 |
| 3D NAND WL Lateral Cavity  | 2–4 μm (lateral) | 15–20 nm (gate height) | 100–200:1 |
| FinFET Gate Trench         | 40–60 nm         | 7–20 nm                | 2–8:1     |

 **3D NAND Word Line note** : WL metal fill occurs through vertical slits into  **horizontal lateral recess cavities** . This is best modeled as an  **Infinite Trench** , not a cylindrical hole.

* L = lateral distance from slit to channel hole (~2–4 μm)
* w = WL gate height (~15–20 nm, roughly half of z-pitch ~40 nm)
* Current z-pitch: ~40 nm; next-gen: 25–30 nm (gate length 10–15 nm)

### 4.2 Sticking Coefficient (s₀)

| Precursor       | Film         | s₀         | Source                 |
| --------------- | ------------ | ---------- | ---------------------- |
| TMA             | Al₂O₃        | ~0.01      | Cremers (2019)         |
| TiCl₄           | TiO₂/TiN     | ~0.006     | Cremers (2019)         |
| TEMAHf          | HfO₂         | ~0.1       | Gordon (2003)          |
| DEZ             | ZnO          | ~0.007     | Gordon (2003)          |
| BDEAS           | SiO₂ (PEALD) | ~3×10⁻⁵    | Literature estimate    |
| MoCl₅ / MoO₂Cl₂ | Mo           | ~0.04–0.05 | Limited data — measure |

### 4.3 Pressure (P)

The calculator supports **Torr or mTorr** (selectable via sidebar radio button, default: Torr). Choose the unit that matches your equipment gauge.

All internal calculations convert to Pa via `Torr_to_Pa()`. For mTorr input: mTorr → ÷1000 → Torr → `Torr_to_Pa()` → Pa.

### 4.4 Fill Tank Parameters

| Parameter       | How to Determine                                    | Typical Range      |
| --------------- | --------------------------------------------------- | ------------------ |
| V_fill (cc)     | Equipment manual + plumbing dead volume             | 10–200 cc          |
| P_before (Torr) | Fill tank gauge after charging                      | 1–50 Torr          |
| P_after (Torr)  | Fill tank gauge after release                       | 50–90% of P_before |
| V_chamber (L)   | Equipment manual or N₂ purge measurement            | 1–50 L             |
| S_pump (L/s)    | 50–80% of spec (effective S). Measure via dP/dt     | 50–500 L/s         |
| C_valve (L/s)   | Skip if fast-fill. Enter if flow restrictor present | 10–100 L/s         |

---

## 5. Practical Examples

### Example 1: DRAM Capacitor HfO₂

1. Select "DRAM Capacitor" preset → adjust to L=4.2 μm, w=75 nm
2. Result: EAR=56, Required ~350 L, Pulse ~0.46 s @ 0.1 Torr
3. Fill Tank: V_fill=50 cc, ΔP=5 Torr → saturation sufficient

### Example 2: 3D NAND Channel Hole Al₂O₃

1. Select "3D NAND Channel Hole" → L=8 μm, w=110 nm, TMA, 0.2 Torr
2. EAR=72.7, Required ~1,200 L → challenging
3. Optimize: P_before=20 Torr, V_fill=100 cc → saturation 110%

### Example 3: Flash WL Mo Fill

1. Select "Flash WL Mo Fill" → **Infinite Trench** (lateral recess)
2. L=3 μm, w=20 nm, EAR=75:1
3. Dose Multiple ~hundreds → multi-pulse strategy needed

---

## 6. Interpreting Results

### 6.1 Result Cards (6 metrics)

| Card              | Meaning                                        |
| ----------------- | ---------------------------------------------- |
| AR                | Geometric aspect ratio                         |
| EAR               | Equivalent aspect ratio (model input)          |
| Required Exposure | Minimum exposure for 100% step coverage [L]    |
| Required Pulse    | Minimum dose time at reference pressure [s]    |
| Knudsen Kn        | Flow regime (> 10 = model reliable)            |
| **Dose Multiple** | HAR required / flat substrate saturation ratio |

### 6.2 Traffic Light Indicators

 **EAR** : 🟢 < 10 / 🟡 10–50 / 🔴 > 50

 **Required Exposure** : 🟢 < 100 L / 🟡 100–1,000 L / 🔴 > 1,000 L

 **Knudsen** : 🟢 > 10 / 🟡 0.1–10 / 🔴 < 0.1

### 6.3 Dose Multiple — The "100–1000× vs. Flat Substrate" Explained

```
Dose Multiple = E_required(HAR) / E_flat(flat substrate)
```

 **Literature basis** : Cremers (2019) — HfO₂ needed 3–43 L on flat substrate vs. 9,000 L in AR≈43 holes (200–3,000×).

 **Critical distinction** :

| Ratio                | Definition            | What it means                                           |
| -------------------- | --------------------- | ------------------------------------------------------- |
| **Dose Multiple**    | E_req(HAR) / E_flat   | "N× harder than flat" →**already in Gordon Required**   |
| **Saturation Ratio** | E_actual / E_req(HAR) | "Is supply sufficient?" →**≥100% = conformal achieved** |

A Dose Multiple of 500× means "500× harder than flat" — do **not** multiply Required Exposure by 500 again.

### 6.4 Saturation Ratio

| Saturation | State                     | Action                          |
| ---------- | ------------------------- | ------------------------------- |
| ≥ 100%     | Conformal                 | Maintain conditions             |
| 90–99%     | Slight thinning at bottom | Slightly increase dose time     |
| 70–89%     | Significant bare area     | Increase time/pressure 1.2–1.5× |
| < 50%      | No meaningful coating     | Fundamental redesign needed     |

### 6.5 Fill Tank Model A vs. B

| Model              | Definition               | Use                         |
| ------------------ | ------------------------ | --------------------------- |
| A: P_eq × t_dose   | Upper bound (constant P) | Optimistic                  |
| B: ODE integration | Real pressure transient  | **Realistic (recommended)** |

---

## 7. Using the Graph Tabs

| Tab                      | Question Answered                                           |
| ------------------------ | ----------------------------------------------------------- |
| Exposure vs EAR          | How much more exposure if structure shrinks? (EAR² scaling) |
| Pulse Time vs P          | What pressure for target dose time?                         |
| Pulse Time vs T          | How does temperature affect dose time?                      |
| EAR Structure Comparison | Hole vs. Trench coating difficulty at same AR?              |

---

## 8. Fill Tank Model Details

### Governing Equations

```
dP_f/dt = −C(P_f − P_c) / V_f
dP_c/dt =  C(P_f − P_c) / V_c − S·P_c / V_c
```

### Key Parameters

* **P_eq = ΔP × V_fill / V_chamber** — equilibrium chamber pressure
* **τ = V_chamber / S_pump** — time constant
* **E_max = P_eq × τ** — maximum achievable exposure

### Decision Criteria

| Condition          | Action                                     |
| ------------------ | ------------------------------------------ |
| E_max > E_required | Set t_dose ≥ t_full                        |
| E_max < E_required | Increase ΔP, V_fill, or decrease S_pump    |
| t_dose/τ < 0.5     | Extending dose time is effective           |
| t_dose/τ > 3       | Diminishing returns; increase ΔP or V_fill |

---

## 9. Unit System & Verification

### 9.1 Pressure Unit Handling

The calculator supports **Torr or mTorr** (sidebar radio button, default: Torr).

Internal conversion path:

```
Torr mode:  user input [Torr]  → P_Torr          → Torr_to_Pa() → P_Pa [Pa]
mTorr mode: user input [mTorr] → ÷1000 → P_Torr  → Torr_to_Pa() → P_Pa [Pa]
```

**All internal calculations use SI units (Pa).** No hardcoded conversion constants (`0.133322`, `/133.322`) exist in the codebase. All conversions go through `Torr_to_Pa()`, `Pa_to_Torr()`, and `Pa_s_to_L()` functions only.

### 9.2 Full Unit Flow

```
Input Unit    SI Conversion          Calculation           Output Unit
──────────   ──────────────────   ────────────────────   ─────────
μm           → m  (×10⁻⁶)        gordon_a → a           [dimensionless]
nm           → m  (×10⁻⁹)
°C           → K  (+273.15)
Torr/mTorr   → Pa (Torr_to_Pa)   E = scale × f(a)       Pa·s → L
g/mol        → kg (×10⁻³/N_A)    t = E / P              s
Å            → m  (×10⁻¹⁰)       λ = kBT/(..d²P)        m → nm
nm, g/cm³    → K_max [1/m²]      Kn = λ/w               [dimensionless]
cc           → m³ (×10⁻⁶)        Fill Tank ODE
L            → m³ (×10⁻³)          P_c(t)                Pa → Torr
L/s          → m³/s (×10⁻³)        E_cum(t)              Pa·s → L
```

### 9.3 Verification

```bash
python run_tests.py     # 44 automated verification tests
```

Verified items:

* Physical constants (k_B, N_A) match CODATA 2018 exact values
* Unit conversion round-trip (Torr↔Pa, Pa·s↔L)
* K_max dimensional analysis: [m]×[g/m³]×[/mol]/[g/mol] = [molecules/m²]
* Required Exposure: two independent calculation paths → exact match
* Torr vs. mTorr input → identical results confirmed
* Fill Tank ODE vs. analytical solution: relative error < 2%

---

## 10. Unit Tests

```bash
python run_tests.py          # Standalone (no streamlit needed)
pytest test_gordon_calculator.py -v   # With pytest
```

44 tests covering:

* Unit conversion accuracy (3)
* Gordon EAR calculation (4)
* K_max estimation (1)
* Mean free path (2)
* Fill Tank ODE vs. analytical (22)
* Penetration depth (3)
* find_t_full (4)
* Input validation (3)
* Gordon exposure formula (2)

---

## 11. Limitations & Caveats

### Model Limitations

1. **s₀ = 1 assumption** : Gordon model gives diffusion-limited result. For s₀ < 1, the actual profile is smoother — model provides a **conservative upper bound**
2. **First-order Langmuir** : Only s(θ) = s₀(1−θ). Complex reaction mechanisms not included
3. **Isothermal** : No temperature gradients inside the structure
4. **Ideal gas** : Valid at low pressures (< few Torr)

### Practical Caveats

1. **Kn < 1** : Gordon model inaccurate. Use CFD simulations
2. **K_max uncertainty** : GPC-based estimate has ±20–50% error. QCM measurement recommended for precision
3. **Precursor decomposition** : CVD component possible at high T. Model assumes pure ALD
4. **Surface chemistry variation** : If s₀ varies with depth (e.g., SiO₂ top vs. Si₃N₄ bottom), accuracy degrades

---

## 12. References

1. R. G. Gordon et al., "A Kinetic Model for Step Coverage by ALD in Narrow Holes or Trenches," *Chem. Vap. Deposition*  **9** , 73–78 (2003)
2. V. Cremers et al., "Conformality in ALD: Current Status Overview of Analysis and Modelling," *Appl. Phys. Rev.*  **6** , 021302 (2019)
3. H. C. M. Knoops et al., "Conformality of Plasma-Assisted ALD: Physical Processes and Modeling," *J. Electrochem. Soc.*  **157** , G241–G249 (2010)
4. S. M. Sze, M. K. Lee,  *Semiconductor Devices: Physics and Technology* , 3rd Ed., Wiley (2012)

---

*ALD Gordon Model Calculator v5 — 2025*
