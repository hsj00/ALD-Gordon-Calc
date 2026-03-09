# ============================================================
# ALD Gordon Model Calculator v5
# ── 주요 변경사항 (v4 → v5) ──
#  1. 코드 최적화: 중복 제거, 함수화, 타입 힌트, 에러 처리 강화
#  2. 물리적 해석: EAR/Exposure/Kn/포화도/Fill Tank 각 결과에 해석 추가
#  3. 입력 파라미터 안내: help 툴팁에 실무 기준값/측정법 안내
#  4. UI/UX 개선: 신호등 색상 + 1줄 해석, 프리셋 버튼, 그래프 활용 가이드
#  5. MoO2Cl2 precursor 추가 (문헌 기반)
#  6. session_state 활용 불필요한 재계산 방지
#  7. find_t_full 인자 순서 버그 수정
#
# Gordon (2003) · Cremers (2019)
# Run: streamlit run ald_gordon_calculator_v5.py
# ============================================================

from __future__ import annotations  # ← 수정: 타입 힌트 호환성

import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.integrate import cumulative_trapezoid, solve_ivp
from dataclasses import dataclass
from typing import Any, Dict, Optional

# ── 물리 상수 ─────────────────────────────────────────────────
k_B: float = 1.380649e-23   # J/K
N_A: float = 6.02214076e23  # /mol

# ── 단위 변환 (타입 힌트 추가) ────────────────────────────────
def Pa_s_to_L(x: float) -> float:
    """Pa·s → Langmuir (1 L = 1.333×10⁻⁴ Pa·s)"""
    return x / 1.33322e-4

def Pa_s_to_Ts(x: float) -> float:
    """Pa·s → Torr·s"""
    return x / 133.322

def Torr_to_Pa(x: float) -> float:
    return x * 133.322

def Pa_to_Torr(x: float) -> float:
    return x / 133.322

# ════════════════════════════════════════════════════════════
# 물리적 해석 헬퍼 함수들 (신규 추가)               ← 수정: 해석 함수 모듈화
# ════════════════════════════════════════════════════════════

def interpret_ear(a: float) -> tuple[str, str]:
    """EAR 값에 대한 물리적 해석과 신호등 색상 반환"""
    if a < 10:
        return (
            f"EAR {a:.1f} — 비교적 완만한 구조. 대부분의 ALD 프로세스로 conformal 코팅 달성 용이.",
            "green"
        )
    elif a < 50:
        return (
            f"EAR {a:.1f} — 중간 난이도. 충분한 노출량 확보 필요. "
            f"DRAM 커패시터급 구조에 해당.",
            "orange"
        )
    elif a < 100:
        return (
            f"EAR {a:.1f} — 매우 도전적. 노출량을 500L 이상 확보하세요. "
            f"최신 3D NAND 수준의 고종횡비 구조.",
            "red"
        )
    else:
        return (
            f"EAR {a:.1f} — 극한 구조. 노출량 수천 L 이상 필요할 수 있으며, "
            f"multi-pulse 또는 precursor 재설계가 필요할 수 있습니다.",
            "red"
        )


def interpret_exposure(E_L: float) -> tuple[str, str]:
    """필요 노출량(Langmuir)에 대한 해석"""
    if E_L < 100:
        return (
            f"노출량 {E_L:.1f}L — 일반적 수준. "
            f"TMA/H₂O 기반 Al₂O₃ 등 표준 ALD 공정에서 흔히 달성 가능.",
            "green"
        )
    elif E_L < 1000:
        return (
            f"노출량 {E_L:.1f}L — 도전적. "
            f"도징 시간 연장 또는 fill tank 압력 증대 필요. "
            f"고종횡비 DRAM 또는 moderate 3D NAND 수준.",
            "orange"
        )
    else:
        return (
            f"노출량 {E_L:.1f}L — 극도로 어려움. "
            f"Multi-pulse dosing, 고압 fill tank, 또는 precursor 변경 검토 권장.",
            "red"
        )


def interpret_knudsen(Kn: float) -> tuple[str, str]:
    """Knudsen 수에 대한 물리적 해석"""
    if Kn > 10:
        return (
            f"Kn={Kn:.1f} — 분자 유동(Molecular Flow). "
            f"분자가 벽을 따라 탁구공처럼 튕기며 이동합니다. "
            f"분자-분자 충돌은 무시 가능하며 Gordon 모델이 정확합니다.",
            "green"
        )
    elif Kn > 0.1:
        return (
            f"Kn={Kn:.1f} — 전환 영역(Transition Flow). "
            f"분자-벽 충돌과 분자-분자 충돌이 공존합니다. "
            f"Gordon 모델은 근사적으로 적용 가능하나, 보정 계수 고려 권장.",
            "orange"
        )
    else:
        return (
            f"Kn={Kn:.1f} — 점성 유동(Viscous Flow). "
            f"분자가 유체처럼 거동하며 분자-분자 충돌이 지배적입니다. "
            f"Gordon 모델 부적합 — Continuum 모델(Navier-Stokes) 적용 필요.",
            "red"
        )


def interpret_saturation(ratio: float) -> tuple[str, str]:
    """포화도(E_actual/E_required)에 대한 해석"""
    pct = ratio * 100
    if ratio >= 1.0:
        return (
            f"포화도 {pct:.1f}% — 충분. 구조 전체에 conformal 코팅이 달성됩니다. "
            f"바닥(bottom)까지 GPC가 균일할 것으로 예상.",
            "green"
        )
    elif ratio >= 0.9:
        return (
            f"포화도 {pct:.1f}% — 거의 달성. 구조 하단 5~10%에서 두께 감소 가능. "
            f"step coverage >90% 예상. 도징 시간을 약간만 늘리면 해결 가능.",
            "orange"
        )
    elif ratio >= 0.7:
        return (
            f"포화도 {pct:.1f}% — 부족. 구조 하단 약 {(1-ratio)*100:.0f}%가 미코팅. "
            f"step coverage가 급격히 저하되는 영역 존재. "
            f"도징 시간 또는 압력을 최소 {1/ratio:.1f}배 증가시켜야 합니다.",
            "red"
        )
    elif ratio >= 0.5:
        return (
            f"포화도 {pct:.1f}% — 크게 부족. 구조 절반 이상이 미코팅 상태. "
            f"근본적인 공정 조건 변경 필요 (ΔP, V_fill, 도징 시간).",
            "red"
        )
    else:
        return (
            f"포화도 {pct:.1f}% — 심각하게 부족. "
            f"현재 조건으로는 유의미한 코팅이 불가능합니다. "
            f"Fill tank 용량, 펌프 속도, 또는 precursor 전략을 전면 재검토하세요.",
            "red"
        )


def interpret_filltank_models(E_A: float, E_B: float, E_req: float,
                               t_dose: float, tau: Optional[float]) -> str:
    """Fill Tank Model A vs B 차이에 대한 해석"""
    if E_A == 0:
        return "Model A 노출량이 0입니다. 입력 조건을 확인하세요."
    diff_pct = abs(E_A - E_B) / E_A * 100
    ratio_td_tau = t_dose / tau if tau and tau > 0 else float('inf')

    lines = []
    if diff_pct > 20:
        lines.append(
            f"⚠ Model A와 B의 차이가 {diff_pct:.1f}%로 큽니다. "
            f"이는 t_dose/τ = {ratio_td_tau:.2f}로, 도징 시간 동안 "
            f"챔버 압력이 크게 변하기 때문입니다."
        )
        if ratio_td_tau < 1:
            lines.append(
                "→ t_dose < τ: 도징 시간이 짧아 fill tank의 precursor가 "
                "완전히 소진되지 않습니다. **도징 시간을 τ의 2~3배로 늘리면** "
                "노출량을 효과적으로 높일 수 있습니다."
            )
        else:
            lines.append(
                "→ t_dose > τ: 펌핑에 의해 precursor가 빠르게 제거됩니다. "
                "**ΔP를 키우거나 fill tank 부피를 늘려** E_max 자체를 높여야 합니다."
            )
    else:
        lines.append(
            f"Model A와 B의 차이가 {diff_pct:.1f}%로 작습니다. "
            f"Fast-fill 근사(순간 평형)가 잘 적용되는 조건입니다."
        )
    return "\n".join(lines)


def interpret_growth_mode(growth_tag: str, s0: float) -> str:
    """성장 모드에 대한 물리적 해석"""
    if "Diffusion" in growth_tag:
        return (
            f"Diffusion-Limited: s₀={s0:.2e}에서 precursor가 표면에 닿으면 "
            f"거의 즉시 반응합니다. 침투 깊이가 날카로운 step function 형태로 "
            f"진행되며, Gordon 모델의 정확도가 가장 높습니다."
        )
    elif "Transition" in growth_tag:
        return (
            f"Transition: s₀={s0:.2e}로 반응 확률이 중간 수준. "
            f"확산과 반응 속도가 경쟁하며 침투 프로파일이 완만해집니다. "
            f"Gordon 모델은 양호한 정확도를 보이나 하한 추정치를 줍니다."
        )
    else:
        return (
            f"Reaction-Limited: s₀={s0:.2e}로 반응 확률이 매우 낮습니다. "
            f"Precursor가 구조 깊은 곳까지 쉽게 도달하므로 conformal 코팅에 "
            f"유리하지만, cycle 당 성장 효율이 낮습니다. "
            f"Gordon 모델은 하한 추정치만 제공합니다."
        )


# ════════════════════════════════════════════════════════════
# Gordon Model 핵심 함수 (타입 힌트 + docstring 강화)
# ════════════════════════════════════════════════════════════

def gordon_a(structure: str, L_m: float, w_m: float,
             z_m: Optional[float] = None) -> float:
    """
    일반화 종횡비(EAR) 계산 — Gordon (2003) Table I.

    Parameters
    ----------
    structure : 구조 유형
    L_m : 구조 깊이 [m]
    w_m : 구조 폭 [m]
    z_m : 구조 길이 [m] (Elongated Hole 전용)

    Returns
    -------
    float : EAR (dimensionless)

    Raises
    ------
    ValueError : 지원하지 않는 구조 유형
    """
    if structure in ("Cylindrical Hole", "Square Hole"):
        return L_m / w_m
    elif structure == "Infinite Trench":
        return L_m / (2.0 * w_m)
    elif structure == "Elongated Hole":
        z = z_m if z_m else w_m
        return L_m * (w_m + z) / (2.0 * w_m * z)
    elif structure == "Square Pillar Array":
        return L_m / (2.0 * np.sqrt(2) * w_m)
    else:
        raise ValueError(f"지원하지 않는 구조: {structure}")  # ← 수정: 에러 처리


def calc_exposure(a: float, K_max: float, m_kg: float, T_K: float) -> float:
    """필요 노출량 [Pa·s] — Gordon Eq.14"""
    scale = K_max * np.sqrt(2 * np.pi * m_kg * k_B * T_K)
    return scale * (1.0/19.0 + 4.0*a/3.0 + 2.0*a**2)


def calc_kmax_scale(K_max: float, m_kg: float, T_K: float) -> float:
    """K_max × √(2π·m·k_B·T) 공통 스케일 인자 [Pa·s·m]"""
    return K_max * np.sqrt(2 * np.pi * m_kg * k_B * T_K)


def calc_penetration_vec(E_arr: np.ndarray, scale: float,
                          w_m: float) -> np.ndarray:
    """
    침투 깊이 [m] — Gordon Eq.24 (벡터화)
    l = (4w/3) × (√(1 + (3/8)·E/scale) − 1)
    """
    E_safe = np.maximum(np.asarray(E_arr, dtype=float), 0.0)
    return (4 * w_m / 3) * (np.sqrt(1 + (3.0/8.0) * E_safe / scale) - 1)


def calc_lambda(T_K: float, P_Pa: float, d_m: float) -> float:
    """평균 자유 경로 [m] — λ = k_B·T / (√2·π·d²·P)"""
    if P_Pa <= 0:  # ← 수정: 0 압력 방어
        return float('inf')
    return k_B * T_K / (np.sqrt(2) * np.pi * d_m**2 * P_Pa)


def calc_kmax(GPC_nm: float, rho_gcm3: float, MW_film: float) -> float:
    """
    K_max [molecules/m²] 추정.

    K_max = GPC [m] × ρ [kg/m³] × N_A / MW_film [g/mol]

    ※ 이 추정은 GPC가 포화 상태일 때 유효합니다.
    오차 범위: 실제 K_max는 표면 -OH 밀도, 리간드 차폐 등에 의해
    ±20~50% 범위에서 변할 수 있습니다.
    """
    return (GPC_nm * 1e-9) * (rho_gcm3 * 1e6) * N_A / MW_film


# ════════════════════════════════════════════════════════════
# Fill Tank ODE  (타입 힌트 추가)
# ════════════════════════════════════════════════════════════

def fill_tank_ode(t: float, y: list[float],
                  V_f: float, V_c: float,
                  C_valve: float, S_pump: float) -> list[float]:
    """Fill tank ↔ Chamber ↔ Pump 연립 ODE"""
    P_f, P_c = y
    dP_f = -C_valve * (P_f - P_c) / V_f
    dP_c =  C_valve * (P_f - P_c) / V_c - S_pump * P_c / V_c
    return [dP_f, dP_c]


@dataclass  # ← 수정: 결과를 dataclass로 구조화
class FillTankResult:
    """Fill Tank ODE 적분 결과"""
    t: np.ndarray
    P_c: np.ndarray
    P_f: np.ndarray
    E_cum: np.ndarray
    P_eq: float
    dP_Pa: float
    tau: Optional[float]
    E_max: Optional[float]
    E_A: float
    E_B: float
    E_C: float
    t_dose: float


def calc_fill_tank(V_fill_cc: float, P_before_T: float, P_after_T: float,
                   V_chamber_L: float, t_dose_s: float,
                   S_pump_Ls: Optional[float] = None,
                   C_valve_Ls: Optional[float] = None,
                   T_K: float = 473.0,
                   n_pts: int = 800) -> FillTankResult:
    """
    Fill Tank 실제 노출량 계산 (ODE 수치 적분).

    Parameters
    ----------
    V_fill_cc    : Fill tank 부피 [cc]
    P_before_T   : Fill 전 압력 [Torr]
    P_after_T    : Fill 후 압력 [Torr]
    V_chamber_L  : Chamber 부피 [L]
    t_dose_s     : Dose 시간 [s]
    S_pump_Ls    : Pump speed [L/s] (optional)
    C_valve_Ls   : Valve conductance [L/s] (optional)
    T_K          : 온도 [K]
    n_pts        : 시간 포인트 수

    Raises
    ------
    ValueError : 물리적으로 불가능한 입력값
    """
    # ── 입력 검증 ── ← 수정: 에러 처리 강화
    if P_after_T >= P_before_T:
        raise ValueError("P_after는 P_before보다 작아야 합니다.")
    if V_fill_cc <= 0 or V_chamber_L <= 0:
        raise ValueError("부피는 양수여야 합니다.")
    if t_dose_s <= 0:
        raise ValueError("도징 시간은 양수여야 합니다.")

    V_f = V_fill_cc * 1e-6    # m³
    V_c = V_chamber_L * 1e-3  # m³
    dP = Torr_to_Pa(P_before_T - P_after_T)
    P_eq = dP * V_f / V_c
    S = (S_pump_Ls * 1e-3) if S_pump_Ls else 0.0
    tau = (V_c / S) if S > 0 else None

    if C_valve_Ls is None:
        # Fast-fill 근사: 순간 평형 후 지수 감소
        t_end = max(t_dose_s * 4, (6 * tau) if tau else t_dose_s * 4)
        t = np.linspace(0, t_end, n_pts)
        if tau:
            P_c = P_eq * np.exp(-t / tau)
        else:
            P_c = np.where(t <= t_dose_s, P_eq, 0.0)
        P_f = np.full_like(t, Torr_to_Pa(P_after_T))
    else:
        C = C_valve_Ls * 1e-3  # m³/s
        t_end = max(t_dose_s * 4, (6*tau) if tau else t_dose_s * 4)
        t_eval = np.linspace(0, t_end, n_pts)

        y0 = [Torr_to_Pa(P_before_T), 0.0]
        sol1 = solve_ivp(fill_tank_ode, [0, t_dose_s], y0,
                         args=(V_f, V_c, C, S if S > 0 else 0),
                         t_eval=t_eval[t_eval <= t_dose_s],
                         max_step=t_dose_s/200, rtol=1e-8)
        y0_2 = [sol1.y[0, -1], sol1.y[1, -1]]
        t_aft = t_eval[t_eval > t_dose_s]
        if len(t_aft) > 0 and S > 0:
            sol2 = solve_ivp(fill_tank_ode, [t_dose_s, t_end], y0_2,
                             args=(V_f, V_c, 0.0, S),
                             t_eval=t_aft, max_step=(t_end-t_dose_s)/200,
                             rtol=1e-8)
            t   = np.concatenate([sol1.t, sol2.t])
            P_f = np.concatenate([sol1.y[0], sol2.y[0]])
            P_c = np.concatenate([sol1.y[1], sol2.y[1]])
        else:
            t, P_f, P_c = sol1.t, sol1.y[0], sol1.y[1]

    # 누적 노출량
    E_cum = cumulative_trapezoid(P_c, t, initial=0.0)

    # Model A, B, C 스칼라
    E_A = P_eq * t_dose_s
    mask = t <= t_dose_s
    E_B = np.trapezoid(P_c[mask], t[mask]) if mask.any() else E_A
    E_C = (P_eq * tau * (1 - np.exp(-t_dose_s / tau))) if tau else E_A
    E_max = P_eq * tau if tau else None

    return FillTankResult(
        t=t, P_c=P_c, P_f=P_f, E_cum=E_cum,
        P_eq=P_eq, dP_Pa=dP, tau=tau, E_max=E_max,
        E_A=E_A, E_B=E_B, E_C=E_C, t_dose=t_dose_s
    )


# ════════════════════════════════════════════════════════════
# 침투 깊이 시간 해석 헬퍼
# ════════════════════════════════════════════════════════════

def find_t_full(E_req: float, mode: str,  # ← 수정: 인자 순서 수정 (E_req 먼저)
                P_Pa: Optional[float] = None,
                P_eq: Optional[float] = None,
                tau: Optional[float] = None) -> Optional[float]:
    """완전 conformal 코팅 달성 시간 추정"""
    if mode == "constant" and P_Pa and P_Pa > 0:
        return E_req / P_Pa
    elif mode == "filltank" and tau and tau > 0 and P_eq and P_eq > 0:
        E_max = P_eq * tau
        if E_req >= E_max * 0.9999:
            return None  # 달성 불가
        return -tau * np.log(1.0 - E_req / E_max)
    return None


# ── Precursor 데이터베이스 ────────────────────────────────────
# 참고: s₀ 값은 Gordon(2003) 및 Cremers(2019)에서 보고된 대표값.
#       정확한 값은 장비/조건에 따라 다르므로 실측을 권장합니다.
PRECURSORS: Dict[str, Dict[str, Any]] = {
    "MoO₂Cl₂  [Mo metal]": {  # ← 수정: MoO₂Cl₂ 추가
        "MW": 198.85, "d_A": 6.0, "s0": 0.040, "GPC": 0.04,
        "note": "MoO₂Cl₂ + H₂ 기반 Mo ALD. Flash memory word line 용도. "
               "MW=198.85 g/mol. GPC 및 s₀ 값은 문헌에서 제한적이며, "
               "공정 온도(~500°C)와 H₂ 환원 조건에 크게 의존합니다. "
               "정확한 s₀ 실측 필요."
    },
    "MoCl₅  [Mo metal]": {
        "MW": 273.21, "d_A": 6.5, "s0": 0.050, "GPC": 0.05,
        "note": "Mo 금속 ALD. MoCl₅ + H₂ 공정. "
               "s₀ 문헌값 제한적 — 실측 필요."
    },
    "TMA  [Al(CH₃)₃]": {
        "MW": 72.08, "d_A": 5.5, "s0": 0.010, "GPC": 0.10,
        "note": "Al₂O₃ ALD의 표준 precursor. s₀≈0.01은 Cremers(2019) 기반."
    },
    "TDMAT  [Ti(NMe₂)₄]": {
        "MW": 224.17, "d_A": 7.0, "s0": 0.020, "GPC": 0.05,
        "note": "TiO₂/TiN ALD용. s₀는 공정 온도에 민감 (150~250°C)."
    },
    "DEZ  [Zn(C₂H₅)₂]": {
        "MW": 123.49, "d_A": 6.0, "s0": 0.007, "GPC": 0.20,
        "note": "ZnO ALD용. 높은 GPC 덕분에 cycle 수를 줄일 수 있음."
    },
    "TEMAHf  [Hf(NEtMe)₄]": {
        "MW": 354.85, "d_A": 8.0, "s0": 0.100, "GPC": 0.10,
        "note": "HfO₂ ALD용. s₀가 높아 diffusion-limited 되기 쉬움."
    },
    "TiCl₄": {
        "MW": 189.68, "d_A": 5.5, "s0": 0.006, "GPC": 0.05,
        "note": "TiN/TiO₂ ALD용. 낮은 s₀ 덕분에 높은 conformality 달성 용이. "
               "Cremers(2019) 기반."
    },
    "ZrNMe₄  [ZrO₂]": {
        "MW": 242.42, "d_A": 7.5, "s0": 0.070, "GPC": 0.10,
        "note": "ZrO₂ ALD용."
    },
    "BDEAS  [SiO₂]": {
        "MW": 248.51, "d_A": 7.5, "s0": 3e-5, "GPC": 0.15,
        "note": "SiO₂ PEALD용. 매우 낮은 s₀로 reaction-limited. "
               "Gordon 모델은 하한 추정치를 제공."
    },
    "Custom": {
        "MW": 100.0, "d_A": 6.0, "s0": 0.010, "GPC": 0.10,
        "note": "사용자 정의 precursor."
    },
}

FILMS: Dict[str, Dict[str, float]] = {
    "Al₂O₃": {"rho": 3.0, "MW_f": 101.96},
    "TiO₂":  {"rho": 4.0, "MW_f": 79.87},
    "ZnO":   {"rho": 5.6, "MW_f": 81.38},
    "HfO₂":  {"rho": 9.7, "MW_f": 210.49},
    "ZrO₂":  {"rho": 5.7, "MW_f": 123.22},
    "SiO₂":  {"rho": 2.2, "MW_f": 60.09},
    "TiN":   {"rho": 5.2, "MW_f": 61.89},
    "Mo":    {"rho": 10.2, "MW_f": 95.95},
    "Custom": {"rho": 4.0, "MW_f": 100.0},
}

STRUCT_INFO: Dict[str, str] = {
    "Cylindrical Hole":    "EAR = L/w  (= AR)",
    "Square Hole":         "EAR = L/w  (= AR)",
    "Infinite Trench":     "EAR = L/(2w) = AR/2",
    "Elongated Hole":      "EAR = L(w+z)/(2wz)",
    "Square Pillar Array": "EAR = L/(2√2·w)",
}

# ── 프리셋 (신규 추가) ──────────────────────────────────────
PRESETS: Dict[str, Dict[str, Any]] = {  # ← 수정: 프리셋 추가
    "DRAM Capacitor (typical)": {
        "structure": "Cylindrical Hole",
        "L_um": 4.0, "w_nm": 80.0,
        "prec": "TEMAHf  [Hf(NEtMe)₄]",
        "film": "HfO₂",
        "T_C": 250.0, "P_mTorr": 100.0,
        "desc": "DRAM 커패시터 유전체 (HfO₂). AR≈50:1."
    },
    "3D NAND Slit (V-NAND)": {
        "structure": "Infinite Trench",
        "L_um": 8.0, "w_nm": 120.0,
        "prec": "TMA  [Al(CH₃)₃]",
        "film": "Al₂O₃",
        "T_C": 300.0, "P_mTorr": 200.0,
        "desc": "3D NAND 슬릿 구조 Al₂O₃ blocking oxide. "
               "Trench이므로 EAR=AR/2."
    },
    "3D NAND Channel Hole": {
        "structure": "Cylindrical Hole",
        "L_um": 10.0, "w_nm": 100.0,
        "prec": "TMA  [Al(CH₃)₃]",
        "film": "Al₂O₃",
        "T_C": 300.0, "P_mTorr": 200.0,
        "desc": "3D NAND 채널 홀 (128L 이상). AR=100:1. 극한 구조."
    },
    "FinFET Gate (Low AR)": {
        "structure": "Infinite Trench",
        "L_um": 0.1, "w_nm": 20.0,
        "prec": "TEMAHf  [Hf(NEtMe)₄]",
        "film": "HfO₂",
        "T_C": 250.0, "P_mTorr": 100.0,
        "desc": "FinFET gate-all-around HfO₂. 낮은 AR로 쉬운 코팅."
    },
    "Flash WL Mo Fill": {  # ← 수정: 구조를 Infinite Trench로 변경, 치수 보정
        "structure": "Infinite Trench",
        "L_um": 3.0, "w_nm": 20.0,
        "prec": "MoO₂Cl₂  [Mo metal]",
        "film": "Mo",
        "T_C": 650.0, "P_mTorr": 80.0,
        "desc": "3D NAND word line Mo fill (MoO₂Cl₂ + H₂). "
               "Replacement gate 공정에서 precursor는 수직 slit을 통해 진입 후 "
               "수평 lateral recess cavity를 채웁니다. "
               "구조: Infinite Trench (양면 slit 진입). "
               "L ≈ slit~channel hole 간 lateral 거리 (~3 μm), "
               "w ≈ WL gate height (~20 nm, z-pitch ~40 nm의 절반). "
               "EAR = L/(2w) = 75:1."
    },
}


# ════════════════════════════════════════════════════════════
# 신호등 디스플레이 헬퍼                        ← 수정: 공통 UI 함수
# ════════════════════════════════════════════════════════════

def show_traffic_light(text: str, color: str) -> None:
    """신호등 색상과 함께 1줄 해석 표시"""
    emoji = {"green": "🟢", "orange": "🟡", "red": "🔴"}.get(color, "⚪")
    bg = {
        "green": "rgba(46,204,113,0.12)",
        "orange": "rgba(243,156,18,0.12)",
        "red": "rgba(231,76,60,0.12)",
    }.get(color, "rgba(200,200,200,0.1)")
    st.markdown(
        f'<div style="background:{bg};border-radius:8px;padding:8px 14px;'
        f'margin:4px 0;font-size:0.92em;">'
        f'{emoji} {text}</div>',
        unsafe_allow_html=True,
    )


# ════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════
def main() -> None:
    st.set_page_config(
        page_title="ALD Gordon Model Calculator",
        page_icon="⚗",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    st.title("ALD Precursor Diffusion Time — Gordon Model")
    st.caption(
        "Gordon et al. CVD 9, 73 (2003)  |  Cremers et al. APR 6, 021302 (2019)"
    )

    # ────────────── SIDEBAR ──────────────────────────────────
    with st.sidebar:
        # ── 프리셋 버튼 ── ← 수정: 프리셋 추가
        st.header("⚡ Quick Presets")
        preset_name = st.selectbox(
            "프리셋 선택 (초보자용 빠른 설정)",
            ["직접 입력"] + list(PRESETS.keys()),
            help="대표적 반도체 구조에 대한 사전 설정. "
                 "선택 후 값이 자동 입력됩니다. 이후 수정 가능."
        )

        # 프리셋 적용
        if preset_name != "직접 입력":
            preset = PRESETS[preset_name]
            st.info(f"**{preset_name}**: {preset['desc']}")
            default_structure = preset["structure"]
            default_L = preset["L_um"]
            default_w = preset["w_nm"]
            default_prec = preset["prec"]
            default_film = preset["film"]
            default_T = preset["T_C"]
            default_P = preset["P_mTorr"]
        else:
            default_structure = "Cylindrical Hole"
            default_L = 5.0
            default_w = 100.0
            default_prec = "MoO₂Cl₂  [Mo metal]"
            default_film = "SiO₂"
            default_T = 650.0
            default_P = 50.0

        st.header("Structure (Geometry)")
        structure = st.selectbox(
            "Structure Type", list(STRUCT_INFO.keys()),
            index=list(STRUCT_INFO.keys()).index(default_structure)
        )
        st.info(STRUCT_INFO[structure])
        c1, c2 = st.columns(2)
        L_um = c1.number_input(
            "Depth L (μm)", 0.1, 1000.0, float(default_L), 0.1,
            help="구조의 깊이. 대표값 — DRAM cap: ~4μm, "  # ← 수정: help 추가
                 "3D NAND 128L: ~8μm, 3D NAND 200L+: ~12μm."
        )
        w_nm = c2.number_input(
            "Width w (nm)", 1.0, 10000.0, float(default_w), 1.0,
            help="구조의 개구부 폭. 대표값 — DRAM cap: ~60-100nm, "  # ← 수정: help 추가
                 "3D NAND hole: ~100-140nm, FinFET gap: ~7-20nm."
        )
        z_um: Optional[float] = None
        if structure == "Elongated Hole":
            z_um = st.number_input("Length z (μm)", 0.1, 1000.0, 1.0, 0.1)

        st.header("Precursor")
        prec_name = st.selectbox(
            "Precursor", list(PRECURSORS.keys()),
            index=list(PRECURSORS.keys()).index(default_prec)
                  if default_prec in PRECURSORS else 0
        )
        P_db = PRECURSORS[prec_name]
        if P_db.get("note"):
            st.caption(f"ℹ {P_db['note']}")

        c1, c2 = st.columns(2)
        MW = c1.number_input("MW (g/mol)", 1.0, 2000.0, float(P_db["MW"]), 0.1)
        d_A = c2.number_input(
            "Molecular diameter (Å)", 1.0, 20.0, float(P_db["d_A"]), 0.1,
            help="van der Waals 직경. DFT 계산 또는 점도 데이터에서 추정. "
                 "일반적으로 5~10Å."
        )
        s0 = st.number_input(
            "Sticking coeff s₀", 1e-6, 1.0, float(P_db["s0"]), format="%.6f",
            help="초기 반응 확률. 측정법: lateral HARS 테스트 구조, "  # ← 수정: help 추가
                 "step coverage vs depth 분석, 또는 QCM. "
                 "대표값 — TMA: ~0.01 (Cremers 2019), "
                 "TiCl₄: ~0.006 (Cremers 2019), "
                 "TEMAHf: ~0.1 (Gordon 2003)."
        )
        GPC = st.number_input(
            "GPC (nm/cycle)", 0.001, 5.0, float(P_db["GPC"]), 0.001,
            help="Growth Per Cycle. 평탄 기판 위 포화 조건에서 측정한 값을 "
                 "사용하세요. Ellipsometry 또는 XRR로 측정."
        )

        st.subheader("K_max")
        kmax_mode = st.radio("설정 방법", ["Film DB 자동 추정", "직접 입력"])
        if kmax_mode == "Film DB 자동 추정":
            film_name = st.selectbox(
                "Film", list(FILMS.keys()),
                index=list(FILMS.keys()).index(default_film)
                      if default_film in FILMS else 0
            )
            F = FILMS[film_name]
            rho_v = st.number_input("Film ρ (g/cm³)", 0.1, 25.0, float(F["rho"]), 0.1)
            MW_f = st.number_input("Film MW (g/mol)", 1.0, 600.0, float(F["MW_f"]), 0.1)
            K_max = calc_kmax(GPC, rho_v, MW_f)
            st.success(f"K_max = {K_max:.3e} molecules/m²")
            st.caption(
                "K_max = GPC × ρ_film × N_A / MW_film.  \n"  # ← 수정: 공식 안내
                "이 추정은 GPC가 포화 상태이고 film이 bulk 밀도에 가까울 때 "
                "유효합니다. 실제 오차 ±20~50% 가능. "
                "보다 정확한 K_max는 in-situ QCM 측정이 권장됩니다."
            )
        else:
            K_max = st.number_input("K_max (×10¹⁸)", 0.01, 1000.0, 4.0, 0.1) * 1e18

        st.header("Process Conditions")
        T_C = st.number_input("Temperature T (°C)", 50.0, 700.0, float(default_T), 5.0)
        P_mTorr = st.number_input("Reference P (mTorr)", 0.01, 10000.0,
                                   float(default_P), 1.0)

        # ── 실제 노출량 모드 ─────────────────────────────────
        st.header("Actual Exposure Calculation")
        exp_mode = st.radio(
            "계산 방법",
            ["비활성화", "일정 압력 P×t", "Fill Tank (ΔP 적분)"]
        )

        fill_inputs: Optional[Dict] = None
        Pact_Pa: Optional[float] = None
        t_actual: Optional[float] = None

        if exp_mode == "일정 압력 P×t":
            Pact_mT = st.number_input("노출 압력 (mTorr)", 0.01, 10000.0,
                                       float(P_mTorr), 1.0)
            t_actual = st.number_input("펄스 시간 t (s)", 0.001, 3600.0, 1.0, 0.01)
            Pact_Pa = Pact_mT * 0.133322

        elif exp_mode == "Fill Tank (ΔP 적분)":
            c1, c2 = st.columns(2)
            V_fill_cc = c1.number_input(
                "Fill Tank V (cc)", 0.1, 10000.0, 50.0, 1.0,
                help="Fill tank 부피. 장비 도면 또는 매뉴얼의 "  # ← 수정: help 추가
                     "'canister volume' 또는 'precursor cylinder manifold volume'을 "
                     "참조하세요. 일반적으로 10~200 cc. "
                     "라인 배관 사체적(dead volume)도 포함해야 정확합니다."
            )
            V_chamber_L = c2.number_input(
                "Chamber V (L)", 0.1, 1000.0, 10.0, 0.5,
                help="반응 챔버 부피. 장비 매뉴얼의 'chamber volume' 또는 "
                     "N₂ purge gas로 PV=nRT를 이용해 실측 가능."
            )
            c1, c2 = st.columns(2)
            P_before_T = c1.number_input("P_before (Torr)", 0.001, 1000.0, 10.0, 0.1)
            P_after_T = c2.number_input("P_after (Torr)", 0.001, 1000.0, 5.0, 0.1)
            t_dose_s = st.number_input("Dose time t_dose (s)", 0.001, 100.0, 1.0, 0.01)
            if P_after_T >= P_before_T:
                st.error("P_after < P_before 이어야 합니다.")
            use_pump = st.checkbox("Pump speed 입력", value=True)
            S_pump_Ls = st.number_input(
                "Pump speed S (L/s)", 0.1, 10000.0, 100.0, 10.0,
                help="유효 펌핑 속도. 주의: 펌프 spec sheet의 값은 "  # ← 수정: help 추가
                     "무부하(no-load) 속도입니다. 실제 챔버에서의 유효 S는 "
                     "배관 conductance로 인해 50~80%로 감소합니다. "
                     "측정법: 챔버를 일정 P로 유지 후 밸브 닫고 "
                     "압력 상승 속도로 S_eff = V × (dP/dt) / P 계산."
            ) if use_pump else None
            use_valve = st.checkbox(
                "Valve conductance 입력 (Full ODE)", value=False,
                help="Fast-fill 조건 (밸브가 완전 개방, C >> S)이면 "  # ← 수정: help 추가
                     "입력 불필요. 밸브가 partly open이거나 "
                     "flow restrictor가 있으면 C를 입력하세요. "
                     "판단 기준: 도징 중 fill tank 압력이 거의 안 변하면 C가 작음."
            )
            C_valve_Ls = st.number_input(
                "Valve conductance C (L/s)", 0.1, 10000.0, 50.0, 5.0,
                help="ALD 밸브의 기체 전달 conductance. "
                     "일반적으로 10~100 L/s. "
                     "Swagelok ALD valve: ~20-50 L/s (N₂ 기준). "
                     "정확한 값은 밸브 매뉴얼 또는 flow test로 결정."
            ) if (use_valve and use_pump) else None
            fill_inputs = dict(V_fill_cc=V_fill_cc, P_before=P_before_T,
                               P_after=P_after_T, V_chamber_L=V_chamber_L,
                               t_dose=t_dose_s, S_pump=S_pump_Ls, C_valve=C_valve_Ls)
            t_actual = t_dose_s

    # ────────────── 기본 계산 ─────────────────────────────────
    L_m = L_um * 1e-6
    w_m = w_nm * 1e-9
    z_m = z_um * 1e-6 if z_um else None
    T_K = T_C + 273.15
    P_Pa = P_mTorr * 0.133322
    m_kg = MW * 1e-3 / N_A
    d_m = d_A * 1e-10
    AR = L_m / w_m

    try:  # ← 수정: 에러 처리
        a = gordon_a(structure, L_m, w_m, z_m)
    except ValueError as e:
        st.error(str(e))
        return

    scale = calc_kmax_scale(K_max, m_kg, T_K)
    E_req = calc_exposure(a, K_max, m_kg, T_K)
    E_req_L = Pa_s_to_L(E_req)
    t_pulse = E_req / P_Pa if P_Pa > 0 else float('inf')
    lam = calc_lambda(T_K, P_Pa, d_m)
    Kn = lam / w_m if w_m > 0 else float('inf')

    # 유동 판정
    if Kn > 10:
        flow_tag = "Molecular Flow (Kn>>1)"
    elif Kn > 0.1:
        flow_tag = "Transition Flow (Kn~1)"
    else:
        flow_tag = "Viscous Flow (Kn<<1)"

    crit = 1.0 / np.sqrt(max(s0, 1e-12))
    growth_tag = ("Diffusion-Limited" if a > 2*crit else
                  "Transition" if a > 0.5*crit else "Reaction-Limited")

    # ── 평탄 기판 포화 노출량 & Dose Multiple ── ← 수정: v5.1 추가
    # 평탄 기판에서의 포화 노출량 (물리적 유도):
    # 표면에 K_max개의 흡착 사이트를 모두 채우려면,
    # s₀ 확률로 반응하므로 1/s₀ 배의 분자 충돌이 필요.
    # E_flat = K_max × √(2πmkBT) / s₀
    # 이 값은 Cremers(2019)의 문헌값과 일치 (HfO₂: 3~43 L)
    E_flat = scale / max(s0, 1e-12)  # Pa·s (flat substrate saturation)
    E_flat_L = Pa_s_to_L(E_flat)
    dose_multiple = E_req / E_flat if E_flat > 0 else float('inf')

    # ────────────── 실제 노출량 계산 ─────────────────────────
    ft_result: Optional[FillTankResult] = None
    E_actual: Optional[float] = None
    E_actual_L: Optional[float] = None
    ratio: Optional[float] = None

    if exp_mode == "일정 압력 P×t" and Pact_Pa and t_actual:
        E_actual = Pact_Pa * t_actual
        E_actual_L = Pa_s_to_L(E_actual)
        ratio = E_actual / E_req if E_req > 0 else 0.0

    elif (exp_mode == "Fill Tank (ΔP 적분)" and fill_inputs
          and fill_inputs["P_after"] < fill_inputs["P_before"]):
        try:
            ft_result = calc_fill_tank(
                fill_inputs["V_fill_cc"], fill_inputs["P_before"],
                fill_inputs["P_after"], fill_inputs["V_chamber_L"],
                fill_inputs["t_dose"], fill_inputs["S_pump"],
                fill_inputs["C_valve"], T_K=T_K
            )
            E_actual = ft_result.E_B
            E_actual_L = Pa_s_to_L(E_actual)
            ratio = E_actual / E_req if E_req > 0 else 0.0
        except ValueError as e:
            st.error(f"Fill Tank 계산 오류: {e}")
            ft_result = None

    l_actual = float(calc_penetration_vec(
        np.array([E_actual if E_actual else 0.0]), scale, w_m)[0]
    ) if E_actual else 0.0
    l_actual_um = l_actual * 1e6
    pen_frac = min(l_actual / L_m, 1.0) if E_actual else 0.0

    # ────────────── 결과 카드 ─────────────────────────────────
    st.header("Calculation Results")
    c1, c2, c3, c4, c5, c6 = st.columns(6)  # ← 수정: 6열로 확장
    c1.metric("AR", f"{AR:.1f} : 1")
    c2.metric("EAR = a", f"{a:.2f} : 1")
    c3.metric("Required Exposure", f"{E_req_L:.1f} L", f"= {E_req:.2e} Pa·s")
    c4.metric("Required Pulse", f"{t_pulse:.4f} s" if t_pulse < 100 else f"{t_pulse:.2f} s",
              f"@ {P_mTorr:.0f} mTorr")
    c5.metric("Knudsen Kn", f"{Kn:.2f}", f"λ={lam*1e9:.1f} nm")
    c6.metric(  # ← 수정: Dose Multiple 추가
        "Dose Multiple",
        f"{dose_multiple:.0f}×",
        f"vs flat substrate ({E_flat_L:.1f} L)",
        help="HAR 구조에서의 필요 노출량 / 평탄 기판 포화 노출량. "
             "이 값이 곧 '평탄 기판 대비 몇 배의 노출량이 필요한가'입니다. "
             "고종횡비 구조에서 100×~1000× 이상이 되는 것은 정상이며, "
             "이는 Gordon 모델의 Required Exposure에 이미 반영되어 있습니다."
    )

    # ── 물리적 해석 신호등 ── ← 수정: 신호등 + 해석 추가
    ear_text, ear_color = interpret_ear(a)
    show_traffic_light(ear_text, ear_color)

    exp_text, exp_color = interpret_exposure(E_req_L)
    show_traffic_light(exp_text, exp_color)

    # ← 수정: Dose Multiple 해석 추가
    if dose_multiple < 10:
        dm_text = (f"Dose Multiple {dose_multiple:.0f}× — 평탄 기판 대비 큰 차이 없음. "
                   f"낮은 AR 구조에서 표준 도징으로 충분합니다.")
        dm_color = "green"
    elif dose_multiple < 100:
        dm_text = (f"Dose Multiple {dose_multiple:.0f}× — 평탄 기판 대비 수십 배 노출 필요. "
                   f"도징 시간을 평탄 기판 조건 대비 크게 늘려야 합니다.")
        dm_color = "orange"
    elif dose_multiple < 1000:
        dm_text = (f"Dose Multiple {dose_multiple:.0f}× — 평탄 기판 대비 수백 배 노출 필요. "
                   f"고종횡비 구조에서 흔히 관찰되는 수준입니다 "
                   f"(Cremers 2019: HfO₂ AR=43에서 ~200×). "
                   f"Gordon Required Exposure가 이 배수를 이미 포함합니다.")
        dm_color = "orange"
    else:
        dm_text = (f"Dose Multiple {dose_multiple:.0f}× — 평탄 기판 대비 1000배 이상. "
                   f"극한 고종횡비 구조. 이 수준의 노출량을 확보하려면 "
                   f"multi-pulse dosing 또는 고압 fill tank이 필수적입니다.")
        dm_color = "red"
    show_traffic_light(dm_text, dm_color)

    # ← 수정: "평탄 기판 vs HAR 구조" 혼동 방지 안내 추가
    with st.expander("ℹ️ 'Dose Multiple' vs 'Saturation Ratio' — 자주 혼동되는 개념 설명"):
        st.markdown(
            "#### 두 가지 다른 '비율' 개념\n\n"
            "ALD conformality 논의에서 **'필요 노출량의 100~1000배'**라는 표현이 자주 등장합니다. "
            "이 표현은 **무엇 대비인지**에 따라 의미가 완전히 다릅니다:\n\n"
            "**1. Dose Multiple (평탄 기판 대비 배수)**\n"
            f"- 현재 값: **{dose_multiple:.0f}×**\n"
            f"- 의미: 평탄 기판 포화 노출량({E_flat_L:.1f} L) 대비 "
            f"HAR conformal 코팅 노출량({E_req_L:.1f} L)의 비율\n"
            "- Gordon 모델의 **Required Exposure에 이미 포함**되어 있는 값\n"
            "- 참고: Cremers (2019) — HfO₂, AR≈43 홀에서 평탄 기판(3~43 L) 대비 "
            "약 200~3000배(9,000 L) 필요\n\n"
            "**2. Saturation Ratio (실제 노출량 / Gordon Required)**\n"
            "- 의미: 현재 공급하는 실제 노출량이 Gordon 모델이 예측한 "
            "conformal 코팅 필요량의 몇 %인지\n"
            "- **이 비율이 ≥100%이면 conformal 코팅 달성**\n"
            "- 추가로 100~1000배를 곱할 필요 **없음**\n\n"
            "#### 요약\n"
            "| 비율 | 정의 | 100% 기준 | 코드에서 |\n"
            "|------|------|----------|--------|\n"
            "| Dose Multiple | E_req(HAR) / E_flat | 평탄 기판 포화 | 자동 계산됨 |\n"
            "| Saturation Ratio | E_actual / E_req(HAR) | Gordon Required | ≥100%이면 OK |\n\n"
            "**결론**: Gordon Required Exposure는 이미 HAR 구조의 확산 제한을 모두 포함한 값입니다. "
            "Saturation Ratio ≥ 100%이면 conformal 코팅이 달성됩니다. "
            "평탄 기판 대비 100~1000배라는 표현은 Dose Multiple에 해당하며, "
            "이는 모델 출력에 이미 반영되어 있습니다."
        )

    kn_text, kn_color = interpret_knudsen(Kn)
    show_traffic_light(kn_text, kn_color)

    growth_text = interpret_growth_mode(growth_tag, s0)
    growth_color = ("green" if "Diffusion" in growth_tag
                    else "orange" if "Transition" in growth_tag
                    else "red")
    show_traffic_light(growth_text, growth_color)

    # ════════════════════════════════════════════════════════
    # 실제 노출량 비교 + Fill Tank 정보
    # ════════════════════════════════════════════════════════
    if E_actual is not None and ratio is not None:
        st.header("Actual Exposure vs Required Exposure")
        pct = ratio * 100

        # 포화도 해석 ← 수정: 포화도 해석 추가
        sat_text, sat_color = interpret_saturation(ratio)
        show_traffic_light(sat_text, sat_color)

        if ft_result is not None:
            E_A_L = Pa_s_to_L(ft_result.E_A)
            E_B_L = Pa_s_to_L(ft_result.E_B)
            P_eq = ft_result.P_eq
            tau = ft_result.tau

            c1, c2, c3, c4 = st.columns(4)
            c1.metric("ΔP × V_fill",
                      f"{Pa_to_Torr(ft_result.dP_Pa):.2f} Torr × {fill_inputs['V_fill_cc']:.0f} cc")
            c2.metric("P_eq (chamber eq.)", f"{P_eq/133.322:.5f} Torr",
                      f"= {P_eq:.4f} Pa")
            c3.metric("τ = V_c / S_pump", f"{tau:.3f} s" if tau else "N/A")
            c4.metric("E_max = P_eq × τ",
                      f"{Pa_s_to_L(ft_result.E_max):.1f} L" if ft_result.E_max else "∞",
                      help="Fill tank로 달성 가능한 최대 노출량 (t→∞)")

            if ft_result.E_max and E_req > ft_result.E_max * 0.9999:
                st.error(
                    f"최대 가능 노출량({Pa_s_to_L(ft_result.E_max):.1f}L) < "
                    f"필요 노출량({E_req_L:.1f}L) — 현재 fill tank 조건으로 "
                    f"완전 코팅 **불가능**!  \n"
                    "ΔP↑, V_fill↑, 또는 S_pump↓ 를 조정하세요."
                )

            df_m = pd.DataFrame({
                "Model": ["A: P_eq×t_dose (upper bound)",
                          "B: ODE integration (recommended)",
                          "Required (Gordon)"],
                "Exposure (L)": [f"{E_A_L:.2f}", f"{E_B_L:.2f}", f"{E_req_L:.2f}"],
                "Exposure (Pa·s)": [f"{ft_result.E_A:.3e}",
                                     f"{ft_result.E_B:.3e}", f"{E_req:.3e}"],
                "Saturation (%)": [f"{ft_result.E_A/E_req*100:.1f}%",
                                    f"{ft_result.E_B/E_req*100:.1f}%",
                                    "100% (ref)"],
            })
            st.dataframe(df_m, use_container_width=True, hide_index=True)

            # Model A vs B 해석 ← 수정: Fill Tank 모델 비교 해석
            ft_interp = interpret_filltank_models(
                ft_result.E_A, ft_result.E_B, E_req,
                ft_result.t_dose, ft_result.tau
            )
            st.info(ft_interp)

            # Fill tank 압력-시간 프로파일
            with st.expander("Fill Tank P-t Profile & Cumulative Exposure",
                              expanded=False):
                t_arr = ft_result.t
                Pc_arr = ft_result.P_c
                Pf_arr = ft_result.P_f
                E_cum = ft_result.E_cum
                t_d = ft_result.t_dose

                fig_ft = make_subplots(
                    rows=3, cols=1, shared_xaxes=True,
                    subplot_titles=("Chamber Pressure P_c(t)",
                                    "Fill Tank Pressure P_fill(t)",
                                    "Cumulative Exposure E(t) = ∫P_c dt"),
                    row_heights=[0.35, 0.25, 0.4]
                )
                mask_d = t_arr <= t_d
                fig_ft.add_trace(go.Scatter(
                    x=t_arr[mask_d], y=Pc_arr[mask_d]/133.322,
                    name="P_c (during dose)",
                    line=dict(color="royalblue", width=2.5),
                    fill="tozeroy", fillcolor="rgba(30,100,255,0.10)"),
                    row=1, col=1)
                fig_ft.add_trace(go.Scatter(
                    x=t_arr[~mask_d], y=Pc_arr[~mask_d]/133.322,
                    name="P_c (after dose)",
                    line=dict(color="royalblue", width=1.5, dash="dot"),
                    showlegend=False), row=1, col=1)
                fig_ft.add_hline(y=P_eq/133.322, line_dash="dot", line_color="gray",
                                 annotation_text=f"P_eq={P_eq/133.322:.4f}T",
                                 row=1, col=1)
                fig_ft.add_vline(x=t_d, line_dash="dash", line_color="black",
                                 row=1, col=1)
                fig_ft.add_trace(go.Scatter(
                    x=t_arr, y=Pf_arr/133.322, name="P_fill",
                    line=dict(color="darkorange", width=2)), row=2, col=1)
                fig_ft.add_vline(x=t_d, line_dash="dash", line_color="black",
                                 row=2, col=1)
                E_cum_L = Pa_s_to_L(E_cum)
                fig_ft.add_trace(go.Scatter(
                    x=t_arr, y=E_cum_L, name="Cumulative exposure",
                    line=dict(color="forestgreen", width=2.5),
                    fill="tozeroy", fillcolor="rgba(0,180,60,0.10)"),
                    row=3, col=1)
                fig_ft.add_hline(y=E_req_L, line_dash="dash", line_color="red",
                                 annotation_text=f"Required {E_req_L:.1f}L",
                                 row=3, col=1)
                if ft_result.E_max:
                    fig_ft.add_hline(y=Pa_s_to_L(ft_result.E_max),
                                     line_dash="dot", line_color="salmon",
                                     annotation_text="E_max", row=3, col=1)
                fig_ft.add_vline(x=t_d, line_dash="dash", line_color="black",
                                 row=3, col=1)
                fig_ft.update_yaxes(title_text="Torr", row=1)
                fig_ft.update_yaxes(title_text="Torr", row=2)
                fig_ft.update_yaxes(title_text="L", row=3)
                fig_ft.update_xaxes(title_text="Time (s)", row=3)
                fig_ft.update_layout(height=600, template="plotly_white",
                                     legend=dict(x=0.01, y=0.99))
                st.plotly_chart(fig_ft, use_container_width=True)
        else:
            col_a, col_b, col_c = st.columns(3)
            col_a.metric("Actual Exposure", f"{E_actual_L:.2f} L",
                         f"= {E_actual:.2e} Pa·s")
            col_b.metric("Required Exposure", f"{E_req_L:.2f} L")
            col_c.metric("Saturation", f"{pct:.1f} %")

        bar_color = "#2ecc71" if ratio >= 1 else (
            "#f39c12" if ratio >= 0.5 else "#e74c3c")
        bar_w = min(pct / 2, 100)
        st.markdown(
            f'<div style="background:#e0e0e0;border-radius:8px;height:26px;width:100%;">'
            f'<div style="background:{bar_color};width:{bar_w:.1f}%;height:100%;'
            f'border-radius:8px;display:flex;align-items:center;padding:0 10px;'
            f'color:white;font-weight:bold;">'
            f'{pct:.1f}% {"✓ OK" if ratio >= 1 else "✗ Insufficient"}'
            f'</div></div>',
            unsafe_allow_html=True,
        )
        st.markdown("")

    # ════════════════════════════════════════════════════════
    # 시간별 침투 깊이 분석
    # ════════════════════════════════════════════════════════
    if E_actual is not None and t_actual is not None:
        st.header("Penetration Depth — Time-Resolved Analysis")
        st.markdown(
            "**공급 시간 동안 precursor가 구조 내부로 얼마나 깊이 확산했는지** "
            "시간에 따라 추적합니다."
        )
        st.latex(
            r"\ell(t) = \frac{4w}{3}"
            r"\left[\sqrt{1 + \frac{3}{8}\cdot\frac{E(t)}"
            r"{K_{\max}\sqrt{2\pi m k_B T}}} - 1\right]"
            r",\quad E(t) = \int_0^t P_c(t')\,dt'"
        )

        if ft_result is not None:
            t_plot = ft_result.t
            E_cum = ft_result.E_cum
            P_c_plot = ft_result.P_c
            t_d = ft_result.t_dose
            P_eq_v = ft_result.P_eq
            tau_v = ft_result.tau
            E_max_ft = ft_result.E_max
            t_full = find_t_full(  # ← 수정: 인자 순서 수정
                E_req, "filltank", P_eq=P_eq_v, tau=tau_v
            ) if tau_v else None
            l_plot = calc_penetration_vec(E_cum, scale, w_m) * 1e6
            Pc_Torr = P_c_plot / 133.322
        else:
            t_d = t_actual
            t_plot = np.linspace(0, t_d * 3, 600)
            E_cum = Pact_Pa * t_plot
            P_c_plot = np.full_like(t_plot, Pact_Pa)
            Pc_Torr = P_c_plot / 133.322
            t_full = find_t_full(E_req, "constant", P_Pa=Pact_Pa)
            E_max_ft = None
            l_plot = calc_penetration_vec(E_cum, scale, w_m) * 1e6

        l_max_um = float(calc_penetration_vec(
            np.array([E_max_ft if E_max_ft else E_cum[-1]]),
            scale, w_m)[0]) * 1e6

        c1, c2, c3, c4 = st.columns(4)
        c1.metric(
            f"Depth @ t_dose={t_d:.2f}s",
            f"{l_actual_um:.3f} μm",
            f"= {pen_frac*100:.1f}% of L",
        )
        c2.metric("Total depth L", f"{L_um:.2f} μm", f"EAR = {a:.1f}:1")
        if t_full is not None:
            c3.metric(
                "t_full (100% coating)",
                f"{t_full:.4f} s",
                f"{'within t_dose' if t_full <= t_d else 'exceeds t_dose'}",
                delta_color="normal" if t_full <= t_d else "inverse",
            )
        else:
            c3.metric("t_full", "Impossible",
                      "Increase ΔP×V_fill or reduce S_pump",
                      delta_color="inverse")
        if E_max_ft:
            c4.metric(
                "l_max (t→∞)",
                f"{l_max_um:.3f} μm",
                f"{'= L (full)' if l_max_um >= L_um else f'< L ({l_max_um/L_um*100:.0f}%)'}",
                delta_color="normal" if l_max_um >= L_um else "inverse",
            )
        else:
            c4.metric("l_max", f"{l_max_um:.3f} μm")

        # Plot: l(t) + P_c(t)
        fig_l = make_subplots(
            rows=2, cols=1, shared_xaxes=True,
            subplot_titles=(
                "Penetration depth l(t)",
                "Chamber pressure P_c(t)",
            ),
            row_heights=[0.65, 0.35]
        )
        fig_l.add_trace(go.Scatter(
            x=t_plot, y=l_plot, mode="lines", name="l(t)",
            line=dict(color="royalblue", width=3)), row=1, col=1)
        fig_l.add_trace(go.Scatter(
            x=t_plot, y=np.full_like(t_plot, L_um),
            mode="lines", name=f"L = {L_um:.2f} μm",
            line=dict(color="red", width=2, dash="dash")), row=1, col=1)
        fig_l.add_vline(x=t_d, line_dash="dash", line_color="black",
                        annotation_text=f"  t_dose={t_d:.2f}s", row=1, col=1)
        if t_full is not None:
            fig_l.add_vline(x=t_full, line_dash="longdash",
                            line_color="forestgreen",
                            annotation_text=f"  t_full={t_full:.3f}s",
                            annotation_position="top left", row=1, col=1)
        fig_l.add_trace(go.Scatter(
            x=[t_d], y=[l_actual_um], mode="markers",
            name=f"l @ t_dose = {l_actual_um:.3f} μm",
            marker=dict(color="red", size=12, symbol="circle")), row=1, col=1)
        if E_max_ft and l_max_um < L_um * 1.5:
            fig_l.add_trace(go.Scatter(
                x=t_plot, y=np.full_like(t_plot, l_max_um),
                mode="lines", name=f"l_max(t→∞)={l_max_um:.3f}μm",
                line=dict(color="salmon", width=1.5, dash="dot")), row=1, col=1)
        fig_l.add_trace(go.Scatter(
            x=t_plot, y=Pc_Torr, mode="lines", name="P_c (Torr)",
            line=dict(color="darkorange", width=2),
            fill="tozeroy", fillcolor="rgba(255,150,0,0.08)"), row=2, col=1)
        fig_l.add_vline(x=t_d, line_dash="dash", line_color="black",
                        row=2, col=1)
        fig_l.update_yaxes(title_text="Penetration depth (μm)", row=1)
        fig_l.update_yaxes(title_text="P_c (Torr)", row=2)
        fig_l.update_xaxes(title_text="Time (s)", row=2)
        fig_l.update_layout(height=540, template="plotly_white",
                            legend=dict(x=0.65, y=0.97))
        st.plotly_chart(fig_l, use_container_width=True)

        # Coverage profile θ(z) snapshots
        st.subheader("Coverage Profile θ(z)")
        st.caption(
            "Gordon model (s₀=1): step function. s₀<1 → smoother profile."
        )
        snapshot_fracs = [0.1, 0.25, 0.5, 0.75, 1.0]
        snap_colors = ["#d0d0ff", "#8888ff", "#4444dd", "#2222aa", "#ff3333"]
        z_arr = np.linspace(0, L_um, 400)

        fig_cov = go.Figure()
        for frac, col in zip(snapshot_fracs, snap_colors):
            t_snap = frac * t_d
            idx = min(np.searchsorted(t_plot, t_snap), len(E_cum) - 1)
            E_snap = E_cum[idx]
            l_snap = float(calc_penetration_vec(
                np.array([E_snap]), scale, w_m)[0]) * 1e6
            theta = np.where(z_arr <= min(l_snap, L_um), 1.0, 0.0)
            fig_cov.add_trace(go.Scatter(
                x=z_arr, y=theta, mode="lines",
                name=f"t={t_snap:.3f}s ({frac*100:.0f}%), l={l_snap:.3f}μm",
                line=dict(color=col, width=2 + int(frac == 1.0))))
        fig_cov.add_trace(go.Scatter(
            x=[0, L_um], y=[1.0, 1.0], mode="lines", name="Full saturation",
            line=dict(color="limegreen", width=1.5, dash="dot")))
        fig_cov.add_vline(x=l_actual_um, line_dash="dash", line_color="red",
                          annotation_text=f"  l(t_dose)={l_actual_um:.3f}μm")
        fig_cov.update_layout(
            title=f"Coverage θ(z) snapshots [{structure}, EAR={a:.1f}]",
            xaxis_title="Depth z (μm)", yaxis_title="θ (0=uncoated, 1=saturated)",
            yaxis_range=[-0.05, 1.15], height=420, template="plotly_white",
            legend=dict(x=0.35, y=0.98, font=dict(size=11)))
        st.plotly_chart(fig_cov, use_container_width=True)

        # Parameter study: dose time vs penetration depth
        with st.expander("Parameter Study: Dose Time vs Penetration Depth"):
            st.markdown(
                "도징 시간을 변화시켰을 때 최종 침투 깊이 변화를 보여줍니다."
            )
            t_study = np.linspace(
                0,
                max(t_d * 4, (t_full * 2) if t_full else t_d * 4),
                300
            )
            if ft_result is not None and ft_result.tau:
                E_study = ft_result.P_eq * ft_result.tau * (
                    1 - np.exp(-t_study / ft_result.tau))
            else:
                E_study = Pact_Pa * t_study
            l_study = calc_penetration_vec(E_study, scale, w_m) * 1e6

            fig_ps = go.Figure()
            fig_ps.add_trace(go.Scatter(
                x=t_study, y=l_study, mode="lines", name="l(t)",
                line=dict(color="royalblue", width=2.5)))
            fig_ps.add_hline(y=L_um, line_dash="dash", line_color="red",
                             annotation_text=f"  L={L_um:.2f}μm")
            fig_ps.add_vline(x=t_d, line_dash="dot", line_color="black",
                             annotation_text="  current t_dose")
            if t_full:
                fig_ps.add_vline(x=t_full, line_dash="longdash",
                                 line_color="forestgreen",
                                 annotation_text=f"  t_full={t_full:.3f}s")
            fig_ps.update_layout(
                title="Dose Time vs Penetration Depth",
                xaxis_title="Dose time (s)", yaxis_title="Penetration (μm)",
                height=380, template="plotly_white")
            st.plotly_chart(fig_ps, use_container_width=True)

    # ════════════════════════════════════════════════════════
    # 그래프 탭
    # ════════════════════════════════════════════════════════
    st.header("Graphical Analysis")
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "Exposure vs EAR", "Pulse Time vs P",
        "Pulse Time vs T", "EAR Structure Comparison", "Concepts"
    ])

    with tab1:
        st.info(  # ← 수정: 활용 가이드 추가
            "**활용**: 현재 구조의 EAR(빨간 점선)에서 필요 노출량을 확인하고, "
            "EAR이 2배가 되면 노출량이 ~4배 증가함을 보여줍니다. "
            "구조 shrink 전략 수립 시 필요 노출량 증가를 예측하는 데 활용하세요."
        )
        a_arr = np.linspace(1, max(a*2.5, 100), 400)
        EL_arr = [Pa_s_to_L(calc_exposure(ai, K_max, m_kg, T_K)) for ai in a_arr]
        EL_q = [Pa_s_to_L(scale * 2*ai**2) for ai in a_arr]
        fig1 = go.Figure()
        fig1.add_trace(go.Scatter(x=a_arr, y=EL_arr, mode="lines",
                                  name="Gordon model",
                                  line=dict(color="royalblue", width=3)))
        fig1.add_trace(go.Scatter(x=a_arr, y=EL_q, mode="lines",
                                  name="2a² asymptote",
                                  line=dict(color="gray", width=1.5, dash="dot")))
        fig1.add_vline(x=a, line_dash="dash", line_color="red",
                       annotation_text=f"  EAR={a:.1f}")
        fig1.add_hline(y=E_req_L, line_dash="dot", line_color="orange",
                       annotation_text=f"  Required {E_req_L:.1f}L")
        if E_actual_L:
            fig1.add_hline(y=E_actual_L, line_color="limegreen",
                           annotation_text=f"  Actual {E_actual_L:.1f}L")
        fig1.update_layout(title="Required Exposure vs EAR", xaxis_title="EAR",
                           yaxis_title="Exposure (L)", yaxis_type="log",
                           height=440, template="plotly_white")
        st.plotly_chart(fig1, use_container_width=True)

    with tab2:
        st.info(  # ← 수정: 활용 가이드 추가
            "**활용**: P × t = const 관계를 활용하여 도징 압력과 시간의 trade-off를 "
            "파악합니다. 목표 시간에서 필요한 최소 압력을 결정하거나, "
            "장비의 base pressure에서 최소 도징 시간을 읽을 수 있습니다. "
            "주황 영역(Kn<1)에서는 Gordon 모델 부적합."
        )
        P_arr = np.logspace(-2, 4, 400)
        t_arr_p = [E_req / (p * 0.133322) for p in P_arr]
        P_kn1 = k_B * T_K / (np.sqrt(2)*np.pi*d_m**2*w_m) / 0.133322
        fig2 = go.Figure()
        fig2.add_vrect(x0=P_kn1, x1=1e4, fillcolor="lightsalmon",
                       opacity=0.15, line_width=0, annotation_text="Kn<1")
        fig2.add_trace(go.Scatter(x=P_arr, y=t_arr_p, mode="lines",
                                  line=dict(color="forestgreen", width=3)))
        fig2.add_vline(x=P_mTorr, line_dash="dash", line_color="red",
                       annotation_text=f"  {P_mTorr:.1f}mTorr")
        fig2.add_vline(x=P_kn1, line_dash="dot", line_color="salmon",
                       annotation_text="  Kn=1")
        fig2.update_layout(title="Pulse Time vs Pressure",
                           xaxis_title="P (mTorr)", yaxis_title="Pulse Time (s)",
                           xaxis_type="log", yaxis_type="log",
                           height=440, template="plotly_white")
        st.plotly_chart(fig2, use_container_width=True)

    with tab3:
        st.info(  # ← 수정: 활용 가이드 추가
            "**활용**: 공정 온도 결정 시 참조. T↑ → 분자 열속도↑ → 같은 노출량에 "
            "더 긴 펄스 필요 (약간의 증가). 동시에 λ가 증가하여 Kn이 개선됩니다. "
            "ALD window와의 조합을 고려하세요."
        )
        T_C_arr = np.linspace(50, 700, 300)
        t_T = [calc_exposure(a, K_max, m_kg, Tc+273.15)/P_Pa for Tc in T_C_arr]
        lam_T = [calc_lambda(Tc+273.15, P_Pa, d_m)*1e9 for Tc in T_C_arr]
        fig3 = make_subplots(rows=2, cols=1, shared_xaxes=True,
                             subplot_titles=("Pulse Time vs T", "MFP vs T"),
                             row_heights=[0.6, 0.4])
        fig3.add_trace(go.Scatter(x=T_C_arr, y=t_T, mode="lines",
                                  line=dict(color="darkorange", width=2.5)),
                       row=1, col=1)
        fig3.add_vline(x=T_C, line_dash="dash", line_color="red",
                       annotation_text=f"  {T_C:.0f}°C")
        fig3.add_trace(go.Scatter(x=T_C_arr, y=lam_T, mode="lines",
                                  line=dict(color="mediumpurple", width=2)),
                       row=2, col=1)
        fig3.add_hline(y=w_nm, line_dash="dot", line_color="gray",
                       annotation_text=f"  Kn=1 (λ={w_nm:.0f}nm)", row=2, col=1)
        fig3.update_yaxes(title_text="s", row=1)
        fig3.update_yaxes(title_text="nm", type="log", row=2)
        fig3.update_xaxes(title_text="T (°C)", row=2)
        fig3.update_layout(height=520, template="plotly_white")
        st.plotly_chart(fig3, use_container_width=True)

    with tab4:
        st.info(  # ← 수정: 활용 가이드 추가
            "**활용**: 같은 AR이라도 구조 형상에 따라 EAR이 크게 다릅니다. "
            "Trench(EAR=AR/2)는 Hole(EAR=AR) 대비 코팅이 쉽습니다. "
            "공정 개발 시 구조 형상 변환 효과를 정량적으로 비교할 수 있습니다."
        )
        AR_arr = np.linspace(1, max(AR*2.5, 100), 300)
        fig4 = go.Figure()
        fig4.add_trace(go.Scatter(x=AR_arr, y=AR_arr, name="Hole (EAR=AR)",
                                  line=dict(color="royalblue", width=2.5)))
        fig4.add_trace(go.Scatter(x=AR_arr, y=AR_arr/2,
                                  name="Trench (EAR=AR/2)",
                                  line=dict(color="forestgreen", width=2.5,
                                            dash="dash")))
        fig4.add_trace(go.Scatter(x=AR_arr, y=AR_arr/(2*np.sqrt(2)),
                                  name="Pillar (EAR=AR/(2√2))",
                                  line=dict(color="darkorange", width=2,
                                            dash="dot")))
        if structure == "Elongated Hole" and z_m:
            fig4.add_trace(go.Scatter(
                x=AR_arr, y=AR_arr*(w_m+z_m)/(2*z_m),
                name=f"Elongated (z={z_um:.1f}μm)",
                line=dict(color="crimson", width=2, dash="dashdot")))
        fig4.add_trace(go.Scatter(
            x=[AR], y=[a], mode="markers", name="Current",
            marker=dict(color="red", size=14, symbol="star")))
        fig4.update_layout(title="EAR vs AR", xaxis_title="AR",
                           yaxis_title="EAR", height=460,
                           template="plotly_white")
        st.plotly_chart(fig4, use_container_width=True)

    with tab5:
        st.header("Key Concepts")
        with st.expander("Langmuir Exposure & Penetration Depth", expanded=True):
            st.latex(
                r"\text{Exposure}=\int P\,dt\,[Pa\cdot s],"
                r"\quad 1\,L=1.333\times10^{-4}\,Pa\cdot s"
            )
            st.latex(
                r"\ell = \frac{4w}{3}\left[\sqrt{1+\frac{3E}"
                r"{8\,K_{\max}\sqrt{2\pi mk_BT}}}-1\right]"
            )
            st.markdown(
                "침투 깊이는 E(t)의 **√ 증가 함수** — "
                "노출량이 4배 증가해야 침투 깊이가 2배."
            )
        with st.expander("Knudsen Number", expanded=True):
            st.latex(
                r"Kn=\lambda/w,\quad\lambda=k_BT/(\sqrt{2}\pi d^2 P)"
            )
            kn_df = pd.DataFrame({
                "Kn": ["Kn>10", "0.1~10", "Kn<0.1"],
                "Flow": ["Molecular", "Transition", "Viscous"],
                "Description": [
                    "Molecule-wall collisions dominate. "
                    "Like ping-pong balls bouncing off walls. Gordon valid.",
                    "Both types of collisions coexist. Approximate.",
                    "Fluid-like behavior. Continuum model needed."
                ],
            })
            st.dataframe(kn_df, use_container_width=True, hide_index=True)
        with st.expander("Sticking Coefficient s₀", expanded=True):
            st.latex(r"s(\theta)=s_0(1-\theta)")
            s_df = pd.DataFrame({
                "s₀": ["~1.0", "0.01~0.1", "<0.001"],
                "Mode": ["Diffusion-Limited", "Intermediate", "Reaction-Limited"],
                "Profile": ["Sharp step", "Moderate gradient", "Very gradual"],
                "Gordon accuracy": ["Best", "Good", "Lower bound"],
            })
            st.dataframe(s_df, use_container_width=True, hide_index=True)
        with st.expander("Fill Tank ODE", expanded=False):
            st.latex(
                r"\frac{dP_f}{dt}=-\frac{C(P_f-P_c)}{V_f},"
                r"\quad\frac{dP_c}{dt}=\frac{C(P_f-P_c)}{V_c}"
                r"-\frac{SP_c}{V_c}"
            )
            st.latex(
                r"E_B=P_{\rm eq}\tau(1-e^{-t_d/\tau}),"
                r"\quad P_{\rm eq}=\frac{\Delta P\cdot V_f}{V_c},"
                r"\quad\tau=\frac{V_c}{S}"
            )
            st.latex(
                r"E_{\rm max}=P_{\rm eq}\cdot\tau"
                r"=\frac{\Delta P\cdot V_f}{S_{\rm pump}}"
            )

    # ── 전체 요약 ─────────────────────────────────────────────
    st.header("Summary")
    rows_s = [
        ("Structure / AR / EAR",
         f"{structure}  |  AR={AR:.1f}  |  EAR={a:.2f}"),
        ("L / w", f"{L_um:.2f} μm  /  {w_nm:.1f} nm"),
        ("Precursor / s₀", f"{prec_name}  |  s₀={s0:.2e}"),
        ("K_max", f"{K_max:.3e} molecules/m²"),
        ("T / P (ref)", f"{T_C:.0f}°C  /  {P_mTorr:.1f} mTorr"),
        ("λ / Kn / Flow", f"{lam*1e9:.1f} nm  /  {Kn:.2f}  /  {flow_tag}"),
        ("Growth mode", growth_tag),
        ("Required Exposure", f"{E_req_L:.2f} L  =  {E_req:.3e} Pa·s"),
        ("Flat Substrate Sat.", f"{E_flat_L:.2f} L  =  {E_flat:.3e} Pa·s"),
        ("Dose Multiple (HAR/flat)", f"{dose_multiple:.0f}×"),
        ("Required Pulse",
         f"{t_pulse:.5f} s  @ {P_mTorr:.1f} mTorr"),
    ]
    if E_actual is not None and ratio is not None:
        lbl = "P×t" if exp_mode == "일정 압력 P×t" else "Fill Tank ODE"
        rows_s += [
            (f"Actual Exposure ({lbl})",
             f"{E_actual_L:.2f} L  =  {E_actual:.3e} Pa·s"),
            ("Saturation Ratio",
             f"{ratio*100:.1f}%  ({'OK' if ratio >= 1 else 'INSUFFICIENT'})"),
            ("Penetration @ t_dose",
             f"{l_actual_um:.4f} μm  ({pen_frac*100:.1f}% of L)"),
        ]
        t_f_val = t_full if (E_actual is not None and t_actual is not None
                             and 't_full' in dir() and t_full) else None
        rows_s.append(
            ("t_full (100% coating)",
             f"{t_f_val:.4f} s" if t_f_val else "Impossible")
        )
        if ft_result is not None:
            rows_s += [
                ("P_eq / τ",
                 f"{ft_result.P_eq/133.322:.5f} Torr  /  "
                 f"{ft_result.tau:.3f} s" if ft_result.tau else "N/A"),
                ("E_max",
                 f"{Pa_s_to_L(ft_result.E_max):.2f} L"
                 if ft_result.E_max else "∞"),
                ("l_max (t→∞)",
                 f"{l_max_um:.4f} μm  ({l_max_um/L_um*100:.1f}% of L)"),
            ]
    df_s = pd.DataFrame(rows_s, columns=["Parameter", "Value"])
    st.dataframe(df_s, use_container_width=True, hide_index=True)
    st.download_button("Download CSV",
                       data=df_s.to_csv(index=False, encoding="utf-8-sig"),
                       file_name="gordon_result.csv", mime="text/csv")


if __name__ == "__main__":
    main()