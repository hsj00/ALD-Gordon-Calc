# Gordon Model ALD Conformality Calculator

A Streamlit app that estimates the **exposure (dose) and pulse time** needed to conformally coat
high-aspect-ratio (HAR) features by ALD, and conversely the **maximum coatable EAR / penetration
depth** achievable with a given exposure budget.

> **References**
> - Review: V. Cremers, R. L. Puurunen, J. Dendooven, *"Conformality in atomic layer deposition,"* **Appl. Phys. Rev. 6, 021302 (2019).**
> - Original model: R. G. Gordon, D. Hausmann, E. Kim, J. Shepard, *Chem. Vap. Depos.* **9(2), 73–78 (2003).**

> ⚠️ **Scope and interpretation**
> This tool is quantitatively meaningful only where **thermal ALD · molecular flow · diffusion-limited
> (s=1) · irreversible reaction** hold. It does **not** apply to PE/ozone (recombination), viscous flow,
> or strongly reaction-limited chemistry. **Treat absolute values as order-of-magnitude estimates.**
> The largest uncertainty comes from the input `K_max`.

---

## Contents

1. [Quick Start](#1-quick-start)
2. [Gordon Model Equations — Meaning of Each Term](#2-gordon-model-equations--meaning-of-each-term)
3. [Inputs — What to Enter, From Where, and On What Basis](#3-inputs--what-to-enter-from-where-and-on-what-basis)
4. [Feature Overview](#4-feature-overview)
5. [Mo ALD Notes](#5-mo-ald-notes-moo₂cl₂--mocl₅)
6. [Validation & Limitations](#6-validation--limitations)
7. [File Layout · Deployment](#7-file-layout--deployment)
8. [Supplementary — Physical Origin of the Equations](#supplementary--physical-origin-of-the-equations-advanced)

---

## 1. Quick Start

```bash
pip install -r requirements.txt
streamlit run app.py
```

If it's your first time, just run the **defaults (HfO₂ hole, AR≈43)**. Forward mode returns a required
exposure of ≈ 9,600 L, reproducing the paper's HfO₂-hole case (~9,000 L).

---

## 2. Gordon Model Equations — Meaning of Each Term

You only need to understand four equations. Focus on **what each term physically represents**, not on the arithmetic.

> Equation numbers **Eq.(14)·Eq.(24) follow Cremers et al. (2019)** (the original model is Gordon 2003).

### 2.1 Equivalent Aspect Ratio (EAR)

```
a = L·p / (4A)        # general form (p: cross-section perimeter, A: cross-section area)
```

| Geometry | Reduction | Intuition |
|---|---|---|
| circular/square hole | `a = L/w` | baseline |
| trench | `a = L/(2w)` | **half** of a hole at the same L/w — easier to coat |
| elongated hole | `a = L(w+z)/(2wz)` | between hole and trench |
| square pillar | `a = L/(2√2·w)` | pillar array, easiest (~1/2.8). ⚠️ **MC result** |

`a` (= EAR) is the "equivalent aspect ratio normalized to a circular hole." **Deeper and narrower → larger a → harder to coat.**

> ⚠️ The **square-pillar formula is a Monte Carlo result, not analytical**, valid only for `w/w_pillar = 3` and `L/w = 5–50` (otherwise extrapolation). See Supplementary S4.

### 2.2 Flat-Surface Saturation Exposure — the reference for everything

```
(Pt)_flat = K_max · √(2π · m · k_B · T)
```

| term | meaning | when larger |
|---|---|---|
| `(Pt)_flat` | minimum exposure to saturate one layer on a flat surface (AR=0) [Pa·s] | — |
| `K_max` | saturated reactant density per unit area [/m²] | more molecules to deposit → exposure ↑ |
| `√(2π·m·k_B·T)` | denominator of the Hertz–Knudsen flux `Φ = P/√(2πmk_BT)` | heavier molecule (m↑) / higher T → slower arrival at fixed pressure → exposure ↑ |

> Exposure `Pt = partial pressure P × pulse time t`, so `(Pt)_flat` is the "pressure×time to coat a flat surface."
> Higher T **increases** the required exposure because, at fixed partial pressure, higher temperature lowers
> the gas density and thus the surface flux (`Φ ∝ P/√(mT)`).
>
> ⚠️ The model treats `K_max` as temperature-independent, but the real saturation density usually **decreases at
> higher T**, partially offsetting the `(Pt)_flat ∝ √T` effect above (another reason to read results as order-of-magnitude).

### 2.3 Required Exposure for Conformal Coating — Eq.(14)

```
Pt = (Pt)_flat · ( 1 + (19/4)·a + (3/2)·a² )
```

| term | value | physical meaning |
|---|---|---|
| `1` | 1 | flat-surface (entrance) contribution |
| `(19/4)·a` | 4.75·a | first-order (linear) entrance/wall correction |
| `(3/2)·a²` | 1.5·a² | **cost of transporting molecules deep into the feature** — dominates at high AR |

> **Key intuition:** at high AR the `(3/2)a²` term dominates → **required exposure ∝ EAR²**.
> Doubling EAR requires roughly **4×** the exposure.

### 2.4 Penetration Depth for Sub-Saturating Exposure — Eq.(24)

```
l = (4w/3) · ( √(1 + (3/8)·E*) − 1 ),     E* ≡ Pt / (Pt)_flat
```

| term | meaning |
|---|---|
| `l` | depth **fully covered** by a given exposure `Pt` [m] |
| `E*` | dimensionless exposure ratio = how many times the flat-saturation dose was applied |
| `w` | feature width — since `l ∝ w`, narrower features penetrate less at the same exposure |

> **Key intuition:** at high exposure `l ∝ √E*` → **doubling penetration depth needs 4× the exposure** (same √ relation as 2.3).
>
> ⚠️ **The `4w/3` coefficient in Eq.(24) is derived for circular holes.** The app applies the same formula (using width `w`)
> to trenches and square pillars, so **penetration depth is an approximation for those geometries**. (The required exposure
> from Eq.(14), being EAR-based, is geometry-general and unaffected.)

### 2.5 Molecular-Flow Check — Knudsen Number

```
λ = k_B·T / (√2 · π · d² · P)        # mean free path (distance traveled between collisions)
Kn = λ / w
```

| Kn | meaning |
|---|---|
| `Kn ≥ 10` | ✅ molecular flow — molecules collide only with walls → **Gordon model applies** |
| `1 ≤ Kn < 10` | ⚠️ transition regime — reduced quantitative reliability |
| `Kn < 1` | 🔴 viscous flow — gas–gas collisions → **Gordon (molecular flow) does not apply** |

> The app computes this automatically and shows a colored badge.

### 2.6 ⚠️ The Two Equations Are Not Exact Inverses

"Required exposure" comes from **Eq.(14)** while "penetration depth" comes from **Eq.(24)**. Their correction
terms differ slightly, giving a discrepancy of about `1 + 0.75a` at **low-to-moderate AR** (negligible at high AR
where the a² term dominates). The app always labels which equation was used. (Math background in the Supplementary.)

---

## 3. Inputs — What to Enter, From Where, and On What Basis

The sidebar is organized into **three groups by data source**. The same notes appear as tooltips (`?`) on each input.

### ① Structure · from cross-section TEM / drawings

| input | meaning | how to decide |
|---|---|---|
| geometry | sets the EAR reduction | choose from device structure / **cross-section TEM** shape |
| depth L | feature depth/height | prefer **cross-section TEM·SEM measurement**, else drawing depth |
| width/gap w | feature width/diameter (gap between pillars for pillar) | prefer **TEM·SEM measurement**, else drawing CD |
| length z | long-axis length of an elongated hole | plan-view TEM / drawing |

> Criterion: **prefer measured (TEM/SEM) over design values.** Post-etch CD/depth govern the actual EAR.

### ② Process Conditions · from recipe / tool

| input | meaning | how to decide |
|---|---|---|
| temperature T | the `√(2πm·k_BT)` flux term | **chuck (substrate) temperature per recipe step** (MoN seed / Mo bulk) |
| partial pressure P | used for pulse time and the Knudsen check | **chamber pressure × precursor fraction** (MFC flow or ampoule vapor pressure). Uncertain for solid precursors |
| exposure budget | (reverse mode) total exposure available | recipe pulse×pressure, or a target value |

### ③ Precursor · Film — from your own measurements

| input | meaning | how to decide |
|---|---|---|
| precursor preset | auto-fills molar mass M | choose your precursor (reference sticking s is also shown) |
| molar mass M | the `√m` flux term | from chemical formula, auto-filled by preset |
| **K_max** | saturated areal density (sets `(Pt)_flat`) | see box below — **the accuracy bottleneck** |
| molecular diameter d | Knudsen check only | ~6–7 Å approximation, **usually leave as is** |
| sticking s | (optional) reaction-limited bracket | entering s<1 adds an upper bound `×(1/s)`. Default 1 |

#### 🔴 K_max — the input that governs accuracy

Enter `K_max` in one of two ways:

1. **Direct input** `[/nm²]` — simplest if known.
2. **Estimate from measured GPC** (recommended) — `K_max ≈ GPC · ρ · N_A / M_film`

| helper input | how to decide |
|---|---|
| measured GPC [nm/cycle] | **XRR thickness ÷ number of cycles** (your measurement) |
| film density ρ [g/cm³] | **XRR measurement recommended**. Preset Mo-series defaults are lattice-based estimates |
| film molar mass M_film [g/mol] | depends on film phase → **confirm phase by XRD** before choosing |

> ⚠️ **GPC-based K_max is a simplified approximation** because "number of adsorption sites ≠ number of
> deposited atoms" (e.g., for TMA/H₂O, reactive OH sites ~7–9/nm² vs. deposited Al ~4.5/nm² per cycle).
> Interpret results at the **order-of-magnitude** level, and prefer your own measured GPC/density/phase.

---

## 4. Feature Overview

- **Forward / reverse modes** — required exposure/pulse time ↔ coatable EAR/penetration depth from a budget. Partial pressure P is a shared input.
- **Four plots** (each with PNG + data-CSV download)
  - exposure vs EAR (Eq.14, log–log)
  - penetration depth vs exposure (Eq.24)
  - penetration depth vs pulse time (feeding time) — at fixed pressure, `l ∝ √t`
  - geometry comparison (reproduces Fig.17: hole > trench > pillar at the same L/w)
- **Multi-cycle EAR evolution** — wall growth narrows the feature each cycle, raising EAR (reproduces Gordon 36→43); penetration depth falls each cycle.
- **Unit toggles** — temperature (°C/K), pressure (Pa/Torr/mTorr/mbar), length (nm/µm/mm), exposure (L/Pa·s/Torr·s).
- **Reliability badges** — Knudsen molecular-flow check, EAR<30 reaction-limited advisory.
- **Reaction-limited bracket** — entering `s<1` in advanced settings shows a `Gordon (lower) ~ ×(1/s) (upper)` range (see below).
- **Export** — results-summary CSV (Excel-friendly BOM) + plot PNG/CSV.

### Meaning of the reaction-limited bracket (s<1)

Because the Gordon equation assumes `s=1` (diffusion-limited), it can underestimate the required exposure
for low-sticking processes. To make this visible, the app shows the rigorous bound as a range:

```
Pt(s=1)  ≤  actual Pt(s)  ≤  Pt(s=1) / s
   (lower, Gordon)            (upper, reaction-limited)
```

- The upper bound is exact at a=0 (flat), while at **high AR the value converges to the lower bound** (the upper
  bound becomes over-conservative) — hence the note "high AR → lower bound, low AR → upper bound."
- For the exact s-dependence within the range, a Monte Carlo model is required (not yet implemented).

---

## 5. Mo ALD Notes (MoO₂Cl₂ / MoCl₅)

- Example process: **MoN seed (MoO₂Cl₂ + NH₃ + H₂, 3-step) → Mo bulk (MoO₂Cl₂ + H₂, 2-step).**
- With a **thermal H₂/NH₃** co-reactant the process is thermal/diffusion-limited, so **Gordon applies.**
- However, interpret results as **upper-bound / order-of-magnitude** due to (also shown in the app's Mo popover):
  1. **HCl by-product** occupying sites → reduced GPC deeper in the feature.
  2. **CVD component** at high temperature (near decomposition).
  3. If an **H₂/NH₃ step is rate-limiting** in the multi-step cycle, actual demand may differ.
  4. **Vapor-pressure instability** of the solid precursor makes `P·t` ill-defined.
- The calculation models the **MoO₂Cl₂ (metal reactant A)** exposure. Use your own **XRR/XRD** values for K_max
  (preset Mo-series densities are lattice-based estimates).

---

## 6. Validation & Limitations

Run `python validation.py` for self-checks against the paper's reference cases (all PASS).

| Check | Computed | Paper/expected |
|---|---|---|
| flat (Pt)_flat (HfO₂, 200°C) | 3.23 L | 3–43 L |
| hole EAR 36 / 43 required exposure | 6,840 / 9,620 L | ~9,000 L |
| hole/trench exposure ratio (L/w=50) | 3.77× | →4 (high AR) |
| forward↔inverse self-consistency | a=20 round-trip | — |
| multi-cycle EAR (self-test) | 36 → 43 | Gordon 36→43 |

**Limitations (stated honestly)**

- **K_max uncertainty dominates accuracy** — the GPC-based estimate is a simplification.
- **Risk of out-of-scope misuse** — PE/ozone (recombination), viscous flow, strong reaction limitation.
- **s-dependence** — currently only a bracket (range). Exact regime behavior needs Monte Carlo (not implemented).
- **Non-inverse equations** (§2.6) — the app always labels which equation was used.

---

## 7. File Layout · Deployment

| file | role |
|---|---|
| `app.py` | Streamlit UI |
| `physics.py` | Gordon-model pure functions (SI) |
| `units.py` | unit conversions |
| `plots.py` | curve data (`curve_*`) + matplotlib figures (`fig_*`) |
| `export.py` | CSV/PNG export |
| `presets.py` | precursor/film presets |
| `multicycle.py` | multi-cycle EAR evolution |
| `validation.py` | reference-case validation |

**Deploy (Streamlit Community Cloud):** push to GitHub → at `share.streamlit.io` select the repo and `app.py` → Deploy.
Subsequent pushes auto-redeploy.

**Korean chart labels (optional):** chart axis labels are English for portability. For Korean, install a Nanum
font and add to the top of `plots.py`:
```python
import matplotlib
matplotlib.rcParams["font.family"] = "NanumGothic"   # or Malgun Gothic / AppleGothic
matplotlib.rcParams["axes.unicode_minus"] = False
```

---

## Supplementary — Physical Origin of the Equations (advanced)

> Not needed to use §2. Read only if you want to know "why 19/4·a and 3/2·a²."
> For the exact derivation, see Gordon (2003) / Cremers (2019) §IV.

### S1. Where sticking hides in the flat-surface term

The Hertz–Knudsen flux is `Φ = P/√(2πmk_BT)`. The time to fill the sites is `t = K_max/(Φ·s)`, so **making the
reaction probability s explicit**:

```
(Pt)_flat = K_max · √(2πmk_BT) / s
```

The Gordon equation assumes `s=1`, dropping this `1/s`. The app's "reaction-limited bracket" restores this `1/s`
as an upper bound. **But dividing all of Eq.(14) by `1/s` is wrong** — the `(3/2)a²` term that dominates at high AR
is *transport (diffusion) limited* and is nearly s-independent (paper Fig.23). The s-dependence is regime-specific
(*low AR = reaction-limited / high AR = s-independent*), so it cannot be inserted as a single factor.

### S2. Origin of the three terms in Eq.(14)

The Gordon model treats the precursor as moving through the feature by **Knudsen diffusion**, reacting
**irreversibly (s=1)** at the walls, and solves the depth-wise transport (conductance) together with wall
consumption to obtain the required exposure:

```
Pt / (Pt)_flat = 1 + (19/4)·a + (3/2)·a²
```

- `1` : one flat (entrance) layer.
- `(19/4)·a` : first-order (linear) correction near the entrance/walls.
- `(3/2)·a²` : cumulative transport cost with depth (second order); dominates at high AR → exposure ∝ a².

The `a²` dominance reflects the Knudsen-transport property that the cost of delivering molecules to the deepest
point grows with the square of the depth.

### S3. Why Eq.(14) and Eq.(24) are not exact inverses

Inverting Eq.(24) with `l = L` for a circular hole (`a = L/w`) gives:

```
E* = 4·a + (3/2)·a²
```

whereas Eq.(14) gives:

```
E* = 1 + (19/4)·a + (3/2)·a²  =  1 + 4.75·a + 1.5·a²
```

→ the difference is `1 + 0.75·a` (an entrance/end-effect–level difference). It is negligible at high AR (a² dominates)
but not at low-to-moderate AR. Hence "required exposure" uses Eq.(14) and "penetration depth" uses Eq.(24),
with the chosen equation always labeled.

### S4. Origin of the geometry EAR reductions

The equivalent aspect ratio `a = L·p/(4A)` is defined by perimeter `p` and cross-section area `A`.
Circular/square holes give `a = L/w`; a trench, blocked only on two sides, gives `a = L/(2w)`.
The **square-pillar `a = L/(2√2·w)` is a Monte Carlo result**, not analytical, valid for `w/wpillar = 3` and `L/w = 5–50`.

### S5. Multi-cycle EAR evolution

Each cycle grows the walls by GPC, narrowing the width to `w_n = w₀ − 2·GPC·n`, so EAR `= L/w_n` increases.
At a fixed per-cycle exposure the penetration depth `l ∝ w_n` decreases each cycle (consistent with Perez/Gordon).
This is simply the Gordon model applied iteratively to a narrowing feature.

---

*Gordon ALD Conformality Calculator · Model: Gordon et al. (2003) · Review: Cremers et al. (2019).
Results depend on the model assumptions (s=1 · molecular flow · diffusion-limited); interpret absolute values at the order-of-magnitude level.*
