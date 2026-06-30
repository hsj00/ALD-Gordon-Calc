"""plots.py — 그래프 + 곡선 데이터 (matplotlib Figure / numpy 배열).

curve_* : 플롯/CSV 공용 데이터(단일 진실 공급원).
fig_*   : curve_* 를 사용해 Figure 생성.
물리 식은 physics.py 와 동일. 축 라벨은 폰트 호환을 위해 영문(설명은 앱에서 한글).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import physics as phys
import units as u

_HOLE = "#2c3e50"; _TRENCH = "#2980b9"; _PILLAR = "#c0392b"
_CURVE = "#1f3a93"; _PEN = "#16a085"; _TIME = "#8e44ad"; _MARK = "#e74c3c"


def _req_L(a, flat):
    """Eq.(14) 벡터화 → Langmuir."""
    Pt = flat * (1.0 + (19.0 / 4.0) * a + (3.0 / 2.0) * a ** 2)
    return u.from_pa_s(Pt, "L")


# ----------------------------- 곡선 데이터 -----------------------------
def curve_exposure_vs_ear(K_max, m, T, ear_min=1.0, ear_max=100.0, n=400):
    flat = phys.flat_saturation_exposure(K_max, m, T)
    a = np.linspace(ear_min, ear_max, n)
    return a, _req_L(a, flat)


def curve_penetration_vs_dose(K_max, m, T, w, dose_max_L, n=400):
    flat = phys.flat_saturation_exposure(K_max, m, T)
    dose_L = np.linspace(1.0, dose_max_L, n)
    E = u.to_pa_s(dose_L, "L") / flat
    l_um = (4.0 * w / 3.0) * (np.sqrt(1.0 + (3.0 / 8.0) * E) - 1.0) * 1e6
    return dose_L, l_um


def curve_penetration_vs_time(K_max, m, T, w, P, t_max, n=400):
    flat = phys.flat_saturation_exposure(K_max, m, T)
    t = np.linspace(t_max / n, t_max, n)
    E = (P * t) / flat
    l_um = (4.0 * w / 3.0) * (np.sqrt(1.0 + (3.0 / 8.0) * E) - 1.0) * 1e6
    ear_hole = l_um / (w * 1e6)
    return t, l_um, ear_hole


def curve_geometry_comparison(K_max, m, T, lw_min=5.0, lw_max=25.0, n=300):
    flat = phys.flat_saturation_exposure(K_max, m, T)
    lw = np.linspace(lw_min, lw_max, n)
    return (lw, _req_L(lw, flat), _req_L(lw / 2.0, flat),
            _req_L(lw / (2.0 * np.sqrt(2.0)), flat))


# ------------------------------- Figure -------------------------------
def fig_exposure_vs_ear(K_max, m, T, current_a=None, ear_min=1.0, ear_max=100.0,
                        budget_dose_L=None):
    a, Pt_L = curve_exposure_vs_ear(K_max, m, T, ear_min, ear_max)
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.loglog(a, Pt_L, color=_CURVE, lw=2, label="Eq.(14)")
    if budget_dose_L is not None and budget_dose_L > 0:
        ax.axhline(budget_dose_L, color="#7f8c8d", ls="--", lw=1.2,
                   label=f"budget ({budget_dose_L:,.0f} L)")
        max_a = phys.max_aspect_ratio(u.to_pa_s(budget_dose_L, "L"), K_max, m, T)
        if max_a > 0:
            ax.axvline(max_a, color="#7f8c8d", ls=":", lw=1)
            ax.annotate(f"max EAR≈{max_a:.0f}", xy=(max_a, budget_dose_L), fontsize=8, color="#555")
    if current_a is not None and current_a > 0:
        flat = phys.flat_saturation_exposure(K_max, m, T)
        ax.scatter([current_a], [_req_L(current_a, flat)], color=_MARK, zorder=6, s=45,
                   label=f"current (EAR={current_a:.0f})")
    ax.set_xlabel("EAR"); ax.set_ylabel("Required exposure (L)")
    ax.set_title("Required exposure vs EAR  [Eq.14]")
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); return fig


def fig_penetration_vs_dose(K_max, m, T, w, current_dose_L=None, dose_max_L=None,
                            target_depth_m=None):
    if dose_max_L is None:
        dose_max_L = (current_dose_L or 1000.0) * 1.5
    dose_L, l_um = curve_penetration_vs_dose(K_max, m, T, w, dose_max_L)
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.plot(dose_L, l_um, color=_PEN, lw=2, label="Eq.(24)")
    if target_depth_m is not None:
        ax.axhline(target_depth_m * 1e6, color="#7f8c8d", ls="--", lw=1.2, label="target depth L")
    if current_dose_L is not None:
        flat = phys.flat_saturation_exposure(K_max, m, T)
        Ec = u.to_pa_s(current_dose_L, "L") / flat
        lc = (4.0 * w / 3.0) * (np.sqrt(1.0 + (3.0 / 8.0) * Ec) - 1.0) * 1e6
        ax.scatter([current_dose_L], [lc], color=_MARK, zorder=6, s=45, label=f"{current_dose_L:,.0f} L")
    ax.set_xlabel("Exposure (L)"); ax.set_ylabel(r"Penetration depth ($\mu$m)")
    ax.set_title("Penetration depth vs exposure  [Eq.24]")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); return fig


def fig_penetration_vs_time(K_max, m, T, w, P, current_t=None, t_max=None,
                            target_depth_m=None, show_ear=True):
    if t_max is None:
        if target_depth_m is not None:
            t_ref = phys.feeding_time_for_depth(target_depth_m, P, w, K_max, m, T)
        else:
            t_ref = current_t or 10.0
        t_max = max(t_ref * 3.0, (current_t or 0.0) * 1.3, 1.0)
    t, l_um, _ = curve_penetration_vs_time(K_max, m, T, w, P, t_max)
    fig, ax = plt.subplots(figsize=(5.4, 3.7))
    ax.plot(t, l_um, color=_TIME, lw=2, label="Eq.(24),  Pt=P·t")
    if target_depth_m is not None:
        ax.axhline(target_depth_m * 1e6, color="#7f8c8d", ls="--", lw=1.2, label="target depth L")
    if current_t is not None and current_t > 0:
        lc = phys.penetration_depth_vs_time(current_t, P, w, K_max, m, T) * 1e6
        ax.scatter([current_t], [lc], color=_MARK, zorder=6, s=45, label=f"t={current_t:.1f} s")
    ax.set_xlabel("Feeding (pulse) time t (s)"); ax.set_ylabel(r"Penetration depth ($\mu$m)")
    ax.set_title(f"Penetration depth vs feeding time  (P={P:g} Pa)")
    if show_ear:
        w_um = w * 1e6
        secax = ax.secondary_yaxis("right", functions=(lambda d: d / w_um, lambda e: e * w_um))
        secax.set_ylabel("coated EAR (hole: l/w)")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); return fig


def fig_geometry_comparison(K_max, m, T, lw_min=5.0, lw_max=25.0):
    lw, hole, trench, pillar = curve_geometry_comparison(K_max, m, T, lw_min, lw_max)
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.semilogy(lw, hole, color=_HOLE, lw=2, marker="*", markevery=40, markersize=9, label="hole")
    ax.semilogy(lw, trench, color=_TRENCH, lw=2, marker="s", markevery=40, label="trench")
    ax.semilogy(lw, pillar, color=_PILLAR, lw=2, marker="o", markevery=40, label="pillar (w/wpillar=3)")
    ax.set_xlabel("L/w"); ax.set_ylabel("Required exposure (L)")
    ax.set_title("Geometry comparison  (cf. Fig.17)")
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); return fig


# --------------------------- 다중 사이클 (확장) ---------------------------
def fig_multicycle_ear(n, ear, w_nm):
    """EAR(좌축) & 피처 폭(우축) vs ALD 사이클."""
    fig, ax = plt.subplots(figsize=(5.4, 3.7))
    ax.plot(n, ear, color=_PILLAR, lw=2)
    ax.set_xlabel("ALD cycle"); ax.set_ylabel("EAR", color=_PILLAR)
    ax.tick_params(axis="y", labelcolor=_PILLAR)
    ax2 = ax.twinx()
    ax2.plot(n, w_nm, color=_TRENCH, lw=2, ls="--")
    ax2.set_ylabel("Pore width (nm)", color=_TRENCH)
    ax2.tick_params(axis="y", labelcolor=_TRENCH)
    ax.set_title("EAR & pore width vs cycle")
    ax.grid(True, alpha=0.3)
    fig.tight_layout(); return fig


def fig_multicycle_penetration(n, l_um):
    """고정 per-cycle 노출에서 사이클별 침투 깊이."""
    fig, ax = plt.subplots(figsize=(5.4, 3.7))
    ax.plot(n, l_um, color=_PEN, lw=2)
    ax.set_xlabel("ALD cycle")
    ax.set_ylabel(r"Penetration depth per cycle ($\mu$m)")
    ax.set_title("Penetration depth vs cycle  (fixed per-cycle dose)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout(); return fig
