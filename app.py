"""app.py — Gordon 모델 ALD 컨포멀리티 계산기 (Phase 5).

실행:  streamlit run app.py

Phase 5 추가:
  - 결과 요약 테이블 + CSV 다운로드
  - 그래프별 PNG 다운로드 + 곡선 데이터 CSV 다운로드
  - 폴리시(레이아웃·푸터·인용)

근거: Cremers et al., Appl. Phys. Rev. 6, 021302 (2019) / Gordon et al. (2003).
물리 계산은 physics.py(순수 함수)에 위임. 본 파일은 UI 전용.
"""
import streamlit as st
import matplotlib.pyplot as plt

import physics as phys
import units as u
import plots
import export
import multicycle
from presets import PRESETS, FILM_PRESETS, MO_PRECURSORS

APP_VERSION = "0.5"

st.set_page_config(page_title="Gordon ALD 컨포멀리티 계산기", layout="centered")

GEOMETRY_LABELS = {
    "원형 홀 (circular hole)": "circular_hole",
    "사각 홀 (square hole)": "square_hole",
    "트렌치 (trench)": "trench",
    "길쭉한 홀 (elongated hole)": "elongated_hole",
    "사각 기둥 배열 (square pillar)": "square_pillar",
}

# ===========================================================================
# 헤더
# ===========================================================================
st.title("Gordon 모델 ALD 컨포멀리티 계산기")
st.caption(
    "근거: Cremers et al., Appl. Phys. Rev. 6, 021302 (2019) / Gordon et al. (2003). "
    "**가정: 열 ALD · 분자 흐름 · 확산 제한(s=1) · 비가역 반응.** "
    "PE/오존(재결합), 점성 흐름, 강한 반응 제한에는 부적용."
)

with st.expander("📖 입력 가이드 — 수식 의미 · 입력값 정하는 기준 (처음이면 펼쳐보세요)"):
    st.markdown(
        "#### 수식 한눈에 (term 의미)\n"
        "- **(Pt)_flat = K_max·√(2πm·k_BT)** — 평탄면 한 겹을 덮는 최소 노출량. "
        "`K_max`=덮을 분자 수, `√(2πmk_BT)`=플럭스 분모(무거운 분자·고온일수록 노출량↑).\n"
        "- **Pt = (Pt)_flat·(1 + 19/4·a + 3/2·a²)** [Eq.14] — 필요 노출량. "
        "`1`=평탄면, `19/4·a`=입구/벽 보정, `3/2·a²`=깊이 수송 비용. "
        "**고AR에서 a² 지배 → 노출량 ∝ EAR²** (EAR 2배 ≈ 노출 4배).\n"
        "- **l = (4w/3)(√(1+3/8·E*)−1)** [Eq.24] — 노출량으로 완전 피복되는 깊이. "
        "`E*=Pt/(Pt)_flat`(무차원 노출비). **l ∝ √E* → 깊이 2배 ≈ 노출 4배**, `l ∝ w`.\n"
        "- **Kn = λ/w, λ=k_BT/(√2πd²P)** — 분자 흐름 점검. Kn≥10 ✅ 적용 / Kn<1 🔴 부적용.\n"
        "- **EAR a** — 같은 L/w라도 트렌치는 홀의 ½, pillar는 1/(2√2)배 (코팅 더 쉬움).\n"
        "\n"
        "ℹ️ '필요 노출량'(Eq.14)과 '침투 깊이'(Eq.24)는 완전 역함수가 아니라 저~중 AR에서 "
        "≈1+0.75a 차이가 날 수 있습니다(고AR은 거의 일치)."
    )
    st.markdown(
        "#### 입력값 어디서 정하나\n"
        "| 그룹 | 입력 | 정하는 기준 |\n"
        "|---|---|---|\n"
        "| ① 구조 | geometry·L·w·z | **단면 TEM·SEM 실측** 우선(없으면 도면 CD/depth) |\n"
        "| ② 공정 | 온도 T | 레시피 **스텝별 척 온도** |\n"
        "| ② 공정 | 분압 P | 챔버압 × 전구체 비율(MFC/ampoule). 고체 전구체는 불확실 |\n"
        "| ③ 측정 | **K_max** | 측정 GPC=**XRR 두께÷사이클**, ρ=**XRR**, 상=**XRD** |\n"
        "| ③ 측정 | 분자 지름 d | ~6–7 Å 근사, 대개 그대로 |\n"
        "| ③ 측정 | sticking s | (선택) s<1이면 반응제한 상한 ×(1/s) 표시 |\n"
        "\n"
        "🔴 **K_max가 정확도 최약점.** GPC 기반은 '흡착 사이트≠증착 원자' 단순화라 차수 해석 — "
        "가능하면 본인 측정 GPC·밀도·상을 쓰세요. 자세한 설명은 README 참조."
    )

# ===========================================================================
# 사이드바 — 입력
# ===========================================================================
with st.sidebar:
    st.header("입력")
    mode = st.radio(
        "계산 모드",
        ["정방향 — 필요 노출량 / 펄스시간",
         "역방향 — 노출량 예산 → 코팅 가능 EAR / 침투깊이"],
        help="정방향: 목표 EAR을 코팅하는 데 필요한 노출량·펄스시간. "
             "역방향: 가진 노출량 예산으로 코팅 가능한 EAR·침투깊이.")
    forward = mode.startswith("정방향")

    # 한 번 정하고 잘 안 바꾸는 항목은 접어둠
    with st.expander("⚙️ 고급 설정 (단위 · 분자지름)", expanded=False):
        temp_unit = st.selectbox("온도 단위", ["°C", "K"], index=0)
        pres_unit = st.selectbox("압력 단위", ["Pa", "Torr", "mTorr", "mbar"], index=0)
        len_unit = st.selectbox("길이 단위", ["nm", "µm", "mm"], index=0)
        exp_unit = st.selectbox("노출량 단위", ["L", "Pa·s", "Torr·s"], index=0)
        d_A = st.number_input(
            "분자 지름 d [Å]", value=7.0, min_value=0.5, step=0.5,
            help="Knudsen 점검 전용 (λ=k_BT/√2πd²P). ~6–7 Å 근사 — 대개 그대로 두면 됩니다.")
        s_stick = st.number_input(
            "sticking s (반응제한 브래킷용)", value=1.0, min_value=0.001, max_value=1.0,
            step=0.05, format="%.3f",
            help="Gordon 기본은 s=1(하한). s<1을 넣으면 반응제한 상한 ≈ ×(1/s)을 함께 표시합니다. "
                 "실제값은 둘 사이 — 고AR(확산제한)은 하한, 저AR(반응제한)은 상한에 근접(Fig.23). "
                 "참고 s값은 ③ 전구체 프리셋에 표시됨. 1이면 브래킷 미표시.")
    d = u.to_m(d_A, "Å")

    # ── ① 구조 · 단면 TEM/도면에서 ─────────────────────────
    st.subheader("① 구조  ·  TEM/도면")
    geom_label = st.selectbox(
        "구조 (geometry)", list(GEOMETRY_LABELS),
        help="EAR 환산식을 결정. 출처: 소자 구조 / 단면 TEM 형상.")
    geom = GEOMETRY_LABELS[geom_label]
    L_in = st.number_input(
        f"깊이/높이 L [{len_unit}]", value=u.from_m(7300e-9, len_unit), min_value=0.0,
        help="피처 깊이/높이. 출처: 단면 TEM·SEM 실측 또는 도면 depth.")
    w_hint = "기둥 사이 갭" if geom == "square_pillar" else "피처 폭/직경"
    w_in = st.number_input(
        f"폭/갭 w [{len_unit}]", value=u.from_m(170e-9, len_unit), min_value=1e-9,
        help=f"{w_hint}. 출처: 단면 TEM·SEM 실측 또는 도면 CD.")
    z_in = None
    if geom == "elongated_hole":
        z_in = st.number_input(
            f"홀 길이 z [{len_unit}]", value=u.from_m(170e-9, len_unit), min_value=1e-9,
            help="elongated hole의 긴 방향 길이. 출처: 평면 TEM/도면.")
    L = u.to_m(L_in, len_unit)
    w = u.to_m(w_in, len_unit)
    z = u.to_m(z_in, len_unit) if z_in is not None else None

    # ── ② 공정 조건 · 레시피/장비에서 ──────────────────────
    st.subheader("② 공정 조건  ·  레시피")
    T_in = st.number_input(
        f"온도 T [{temp_unit}]", value=200.0 if temp_unit == "°C" else 473.15, step=10.0,
        help="기판(척) 온도. flux의 √(2πm·k_BT) 항. 출처: 레시피 스텝별 척온도.")
    T = u.celsius_to_kelvin(T_in) if temp_unit == "°C" else T_in
    P_in = st.number_input(
        f"분압 P [{pres_unit}]", value=u.from_pa(0.27, pres_unit), min_value=0.0,
        help="reactant 분압. 펄스시간 t=Pt/P·Knudsen 점검에 사용. "
             "출처: 챔버압×전구체 비율(MFC/ampoule). 고체 전구체는 불확실.")
    P = u.to_pa(P_in, pres_unit)
    if not forward:
        dose_in = st.number_input(
            f"노출량 예산 [{exp_unit}]",
            value=u.from_pa_s(u.to_pa_s(9000.0, "L"), exp_unit), min_value=0.0,
            help="쓸 수 있는 총 노출량. 출처: 레시피 펄스×분압, 또는 목표치.")
        dose_pa_s = u.to_pa_s(dose_in, exp_unit)

    # ── ③ 전구체·막 · 본인 측정값에서 ──────────────────────
    st.subheader("③ 전구체·막  ·  측정")
    preset_name = st.selectbox(
        "전구체 프리셋", list(PRESETS),
        help="전구체 선택 → 몰질량 M 자동. 참고 sticking s 표시(Gordon은 s=1, 계산 미반영).")
    preset = PRESETS[preset_name]
    s_note = "" if preset["s"] in ("—", "미공개") else f"  ·  참고 s≈{preset['s']}"
    M = st.number_input(
        "전구체 몰질량 M [g/mol]", value=float(preset["M"]), min_value=1.0, step=1.0,
        key=f"M_{preset_name}",
        help="전구체 1분자 몰질량(flux의 √m 항). 프리셋이 자동 채움." + s_note)
    if preset_name in MO_PRECURSORS:
        with st.popover("🧪 Mo 공정 주의 (펼치기)"):
            st.warning(
                "**Mo 염화물 공정** — co-reactant가 **thermal H₂/NH₃**이면 Gordon(열·확산제한) "
                "적용 가능. 단, ① HCl 부산물 사이트 점유로 깊은 곳 GPC 저하, ② 고온에서 CVD 성분, "
                "③ 다단계(2-step Mo / 3-step MoN seed) H₂·NH₃ 율속이면 실제 요구 노출 상이, "
                "④ 고체 precursor 증기압 불안정으로 P·t 불확실. 본 계산은 **MoO₂Cl₂(reactant A)** "
                "노출 기준 — 결과는 상한·차수로 해석.")

    kmax_mode = st.radio(
        "K_max 입력", ["직접 입력 [/nm²]", "측정 GPC로 추정"],
        help="K_max=포화 면적당 분자수 → (Pt)_flat=K_max·√(2πm·k_BT). "
             "본인 측정 GPC가 있으면 'GPC로 추정' 권장.")
    if kmax_mode.startswith("직접"):
        K_max_nm2 = st.number_input(
            "포화 면밀도 K_max [/nm²]", value=2.77, min_value=1e-6, step=0.1,
            help="⚠️ 흡착 사이트 수 ≠ 증착 원자 수(논문 지적). 차수 수준 해석.")
    else:
        film_name = st.selectbox(
            "막 (film)", list(FILM_PRESETS),
            help="막 선택 시 ρ·M_film 자동. 상(phase) 출처: XRD.")
        film = FILM_PRESETS[film_name]
        gpc_nm = st.number_input(
            "측정 GPC [nm/cycle]", value=float(film["gpc"]), min_value=1e-6,
            step=0.005, format="%.3f", key=f"gpc_{film_name}",
            help="사이클당 두께. 출처: XRR 두께 ÷ 사이클 수 (본인 측정).")
        rho = st.number_input(
            "막 밀도 ρ [g/cm³]", value=float(film["rho"]), min_value=0.1, step=0.1,
            key=f"rho_{film_name}",
            help="출처: XRR 측정 권장. Mo계 기본값은 격자 기반 추정.")
        M_film = st.number_input(
            "막 몰질량 M_film [g/mol]", value=float(film["M"]), min_value=1.0, step=1.0,
            key=f"Mf_{film_name}", help="막 화학식단위 몰질량. 상(phase)에 따라 다름.")
        K_max_nm2 = u.areal_density_per_m2_to_per_nm2(
            phys.kmax_from_gpc(gpc_nm * 1e-9, rho * 1000.0, M_film / 1000.0))
        st.info(f"추정 K_max = **{K_max_nm2:.2f} /nm²**"
                + (f"  ·  {film['note']}" if film['note'] else ""))

# ===========================================================================
# 공통 계산
# ===========================================================================
m = u.molar_mass_to_molecule_mass(M)
K_max = u.areal_density_per_nm2_to_per_m2(K_max_nm2)
a = phys.aspect_ratio(geom, L, w, z)
is_hole = geom in ("circular_hole", "square_hole")
ar_simple = phys.aspect_ratio_simple(L, w)
flat_pa_s = phys.flat_saturation_exposure(K_max, m, T)
Kn = phys.knudsen_number(phys.mean_free_path(T, P, d), w) if P > 0 else None
t_reach = phys.feeding_time_for_depth(L, P, w, K_max, m, T) if P > 0 else None


def fmt_exp(pa_s):
    return f"{u.from_pa_s(pa_s, exp_unit):,.4g} {exp_unit}"


# 모드별 출력값 일괄 계산
if forward:
    Pt = phys.required_exposure(a, K_max, m, T)
    pulse_t = phys.pulse_time(Pt, P) if P > 0 else None
else:
    Pt = dose_pa_s
    max_a = phys.max_aspect_ratio(Pt, K_max, m, T)
    l = phys.penetration_depth(Pt, w, K_max, m, T)
    coated_a = phys.aspect_ratio(geom, max(l, 0.0), w, z) if l > 0 else 0.0
    spend_t = phys.pulse_time(Pt, P) if P > 0 else None

# ===========================================================================
# 결과 — 구조 + 신뢰도 배지
# ===========================================================================
st.subheader("구조")
c1, c2, c3 = st.columns(3)
c1.metric("AR = L/w (단면)", f"{ar_simple:.1f}")
c2.metric("EAR (등가, 기하반영)", f"{a:.1f}")
c3.metric("평탄면 포화 (Pt)_flat", fmt_exp(flat_pa_s))
if geom == "square_pillar":
    st.caption("⚠️ pillar EAR식 `L/(2√2·w)`는 MC 결과 — `w/w_pillar=3`, `L/w=5–50`에서만 유효(그 밖은 외삽).")

if P > 0:
    msg = (f"Knudsen 수 Kn = λ/w = **{Kn:.3g}**  "
           f"(λ={u.from_m(phys.mean_free_path(T, P, d),'µm'):.3g} µm, w={w_in:g} {len_unit})")
    if Kn >= 10:
        st.success("✅ 분자 흐름 영역 — Gordon 모델 적용 가능.  " + msg)
    elif Kn >= 1:
        st.warning("⚠️ 전이 영역 (1 ≤ Kn < 10) — 정량 신뢰도 저하.  " + msg)
    else:
        st.error("🔴 점성 흐름 (Kn < 1) — Gordon(분자 흐름) 모델 부적용.  " + msg)
else:
    st.caption("분압 P를 입력하면 분자 흐름(Knudsen) 점검을 표시합니다.")

if a < 30:
    st.info("ℹ️ EAR < 30 영역에서는 sticking probability 의존성이 큽니다. "
            "실제 reactant는 s<1 이므로 필요 노출량이 본 추정(s=1)보다 클 수 있습니다 (cf. Fig.23).")

st.divider()

# ===========================================================================
# 결과 — 모드별
# ===========================================================================
if forward:
    st.subheader("정방향 결과 — 필요 노출량 [Eq.14]")
    c1, c2 = st.columns(2)
    c1.metric("필요 노출량", fmt_exp(Pt))
    c2.metric("필요 노출량 (SI)", f"{Pt:.3e} Pa·s")
    if pulse_t is not None:
        t_str = (f"{pulse_t:,.1f} s" if s_stick >= 1.0
                 else f"{pulse_t:,.1f} ~ {pulse_t / s_stick:,.1f} s")
        st.metric(f"필요 펄스 시간 (P={P_in:g} {pres_unit}, Eq.14 기준)", t_str)
    st.caption(f"Pt = (Pt)_flat × (1 + 19/4·a + 3/2·a²) = {fmt_exp(flat_pa_s)} × "
               f"(1 + {19/4*a:.1f} + {1.5*a*a:.1f})")
    if s_stick < 1.0:
        st.warning(
            f"🔁 **반응제한 브래킷 (s={s_stick:g})**: 필요 노출량 "
            f"**{fmt_exp(Pt)} ~ {fmt_exp(Pt / s_stick)}** "
            f"(하한=Gordon s=1 / 상한≈×1/s). 실제값은 둘 사이 — "
            f"고AR(확산제한)은 하한, 저AR(반응제한)은 상한에 근접 (Fig.23).")
else:
    st.subheader("역방향 결과 — 노출량 예산으로 가능한 한계")
    c1, c2 = st.columns(2)
    c1.metric("코팅 가능 최대 EAR [Eq.14 역산]", f"{max_a:.1f}")
    c2.metric("이 구조 EAR (목표)", f"{a:.1f}",
              delta="달성" if a <= max_a else "부족",
              delta_color="normal" if a <= max_a else "inverse")
    c3, c4 = st.columns(2)
    c3.metric("완전 피복 침투 깊이 l [Eq.24]", f"{u.from_m(l,'µm'):.2f} µm")
    c4.metric("→ 침투 깊이 기준 coated EAR", f"{coated_a:.1f}")
    if not is_hole:
        st.caption("⚠️ 침투 깊이(Eq.24, 4w/3)는 **원형 홀 기준 유도** — 트렌치·pillar에서는 근사입니다.")
    if spend_t is not None:
        st.metric(f"예산 소진 시간 (P={P_in:g} {pres_unit})", f"{spend_t:,.1f} s")
    if max_a == 0.0:
        st.warning("노출량 예산이 평탄면 포화(Pt_flat)에도 못 미칩니다.")
    st.info("ℹ️ 두 식은 완전한 역함수가 아닙니다(§2.5). '코팅 가능 EAR'=Eq.14 역산, "
            "'침투 깊이'=Eq.24. 저~중간 AR에서 차이(≈1+0.75a) 가능, 고AR은 거의 일치.")
    if s_stick < 1.0:
        max_a_lo = phys.max_aspect_ratio(dose_pa_s * s_stick, K_max, m, T)
        l_lo = phys.penetration_depth(dose_pa_s * s_stick, w, K_max, m, T)
        st.warning(
            f"🔁 **반응제한 브래킷 (s={s_stick:g})**: 같은 예산으로 코팅 가능 EAR "
            f"**{max_a_lo:.1f} ~ {max_a:.1f}**, 침투깊이 "
            f"**{u.from_m(l_lo,'µm'):.2f} ~ {u.from_m(l,'µm'):.2f} µm** "
            f"(낮은 쪽=반응제한 s={s_stick:g} / 높은 쪽=Gordon s=1). 고AR은 s=1 쪽에 근접.")

# ===========================================================================
# 결과 요약 + CSV 다운로드
# ===========================================================================
st.divider()
st.subheader("📋 결과 요약 & 다운로드")

rows = [
    ("생성시각", export.timestamp(), ""),
    ("모드", "정방향" if forward else "역방향", ""),
    ("온도 T", T_in, temp_unit),
    ("구조", geom_label, ""),
    ("깊이 L", L_in, len_unit),
    ("폭/갭 w", w_in, len_unit),
]
if z_in is not None:
    rows.append(("홀 길이 z", z_in, len_unit))
rows += [
    ("reactant", preset_name, ""),
    ("몰질량 M", M, "g/mol"),
    ("K_max", K_max_nm2, "/nm²"),
    ("분압 P", P_in, pres_unit),
    ("AR (L/w)", ar_simple, ""),
    ("EAR", a, ""),
    ("(Pt)_flat", u.from_pa_s(flat_pa_s, exp_unit), exp_unit),
]
if Kn is not None:
    rows.append(("Knudsen Kn", Kn, ""))
if s_stick < 1.0:
    rows.append(("sticking s (브래킷)", s_stick, ""))
if forward:
    rows += [
        ("필요 노출량", u.from_pa_s(Pt, exp_unit), exp_unit),
        ("필요 노출량(SI)", Pt, "Pa·s"),
    ]
    if s_stick < 1.0:
        rows.append(("필요 노출량 반응제한 상한(×1/s)", u.from_pa_s(Pt / s_stick, exp_unit), exp_unit))
    if pulse_t is not None:
        rows.append(("필요 펄스시간(Eq.14)", pulse_t, "s"))
    if t_reach is not None:
        rows.append(("목표깊이 도달 feeding time(Eq.24)", t_reach, "s"))
else:
    rows += [
        ("노출량 예산", dose_in, exp_unit),
        ("코팅 가능 최대 EAR(Eq.14)", max_a, ""),
        ("침투깊이 l(Eq.24)", u.from_m(l, "µm"), "µm"),
        ("coated EAR", coated_a, ""),
    ]
    if s_stick < 1.0:
        rows.append(("코팅 가능 EAR 반응제한 하한", phys.max_aspect_ratio(dose_pa_s * s_stick, K_max, m, T), ""))
    if spend_t is not None:
        rows.append(("예산 소진 시간", spend_t, "s"))

st.table([{"항목": p, "값": export._g(v), "단위": un} for p, v, un in rows])
st.download_button("⬇️ 결과 요약 CSV", data=export.summary_csv(rows),
                   file_name="ald_summary.csv", mime="text/csv", key="dl_summary")

# ===========================================================================
# 그래프 (4탭) + 그래프별 PNG·CSV 다운로드
# ===========================================================================
st.divider()
st.subheader("📈 그래프")

ear_max = max(100.0, a * 1.3)
if forward:
    dose_max_L = u.from_pa_s(Pt, "L") * 1.5
    budget_L = None
    cur_dose_L = u.from_pa_s(Pt, "L")
    cur_t = t_reach
else:
    cur_dose_L = u.from_pa_s(dose_pa_s, "L")
    budget_L = cur_dose_L
    dose_max_L = max(cur_dose_L, u.from_pa_s(phys.required_exposure(a, K_max, m, T), "L")) * 1.3
    cur_t = spend_t
t_max = max((t_reach or 1.0) * 3.0, (cur_t or 0.0) * 1.3, 1.0) if P > 0 else None


def dl_buttons(fig, colnames, columns, base):
    """그래프 아래 PNG/CSV 다운로드 버튼 한 쌍."""
    cc1, cc2 = st.columns(2)
    cc1.download_button("⬇️ PNG", data=export.fig_png_bytes(fig),
                        file_name=f"{base}.png", mime="image/png", key=f"png_{base}")
    cc2.download_button("⬇️ 데이터 CSV", data=export.series_csv(colnames, columns),
                        file_name=f"{base}.csv", mime="text/csv", key=f"csv_{base}")


tabs = st.tabs(["노출량 – EAR", "침투깊이 – 노출량", "침투깊이 – feeding time", "기하 비교 (Fig.17)"])

with tabs[0]:
    figA = plots.fig_exposure_vs_ear(K_max, m, T, current_a=a, ear_max=ear_max,
                                     budget_dose_L=budget_L, s_stick=s_stick)
    st.pyplot(figA)
    xa, ya = plots.curve_exposure_vs_ear(K_max, m, T, ear_max=ear_max)
    dl_buttons(figA, ["EAR", "Required_exposure_L"], [xa, ya], "exposure_vs_EAR")
    plt.close(figA)
    st.caption("Eq.(14): 노출량 ∝ EAR²(고AR). 빨간 점 = 현재 구조."
               + ("" if s_stick >= 1.0 else f"  음영=반응제한 범위(×1/s, s={s_stick:g})."))

with tabs[1]:
    figB = plots.fig_penetration_vs_dose(K_max, m, T, w, current_dose_L=cur_dose_L,
                                         dose_max_L=dose_max_L, target_depth_m=L, s_stick=s_stick)
    st.pyplot(figB)
    xb, yb = plots.curve_penetration_vs_dose(K_max, m, T, w, dose_max_L)
    dl_buttons(figB, ["Exposure_L", "Penetration_depth_um"], [xb, yb], "pendepth_vs_dose")
    plt.close(figB)
    st.caption("Eq.(24): 주어진 노출량에서 완전 피복 깊이. 점선 = 목표 깊이 L."
               + ("" if is_hole else "  ⚠️ 4w/3는 원형 홀 기준 — 트렌치/pillar에서는 근사."))

with tabs[2]:
    if P > 0:
        figD = plots.fig_penetration_vs_time(K_max, m, T, w, P, current_t=cur_t,
                                             t_max=t_max, target_depth_m=L, s_stick=s_stick)
        st.pyplot(figD)
        xd, yd, ed = plots.curve_penetration_vs_time(K_max, m, T, w, P, t_max)
        dl_buttons(figD, ["Feeding_time_s", "Penetration_depth_um", "coated_EAR_hole"],
                   [xd, yd, ed], "pendepth_vs_time")
        plt.close(figD)
        tr_str = (f"{t_reach:,.1f} s" if s_stick >= 1.0
                  else f"{t_reach:,.1f} ~ {t_reach / s_stick:,.1f} s")
        st.metric(f"목표 깊이 L={u.from_m(L,'µm'):.2f} µm 도달 feeding time (P={P_in:g} {pres_unit})",
                  tr_str)
        st.caption("Eq.(24)를 Pt=P·t 로 치환. l ∝ √t (깊이 2배 ≈ 시간 4배). "
                   "⚠️ P 일정 가정 — 반응기 고갈·펄스 rise/fall 무시, 실측 시간은 더 클 수 있음."
                   + ("" if is_hole else " 침투 깊이 4w/3는 원형 홀 기준(트렌치/pillar 근사)."))
    else:
        st.warning("분압 P > 0 을 입력하면 feeding time 그래프가 표시됩니다.")

with tabs[3]:
    figC = plots.fig_geometry_comparison(K_max, m, T)
    st.pyplot(figC)
    lw, hL, tL, pL = plots.curve_geometry_comparison(K_max, m, T)
    dl_buttons(figC, ["L_over_w", "hole_L", "trench_L", "pillar_L"], [lw, hL, tL, pL], "geometry_comparison")
    plt.close(figC)
    st.caption("같은 L/w라도 트렌치는 ~1/4, pillar는 ~1/8 노출량으로 코팅(고AR 근사).")

# ===========================================================================
# 다중 사이클 EAR 변화 (확장)
# ===========================================================================
st.divider()
st.subheader("🔁 다중 사이클 EAR 변화 (확장)")
st.caption("매 사이클 벽에 GPC가 쌓여 폭이 줄고 EAR이 증가(Gordon 36→43 / Perez). "
           "고정 per-cycle 노출에서 침투 깊이는 사이클마다 감소. ⚠️ 등각·일정 GPC 가정, "
           "침투 깊이는 Eq.(24)(s=1) 기반.")

if st.checkbox("다중 사이클 분석 활성화", value=False):
    m1, m2, m3 = st.columns(3)
    gpc_mc = m1.number_input("GPC [nm/cycle]", value=0.10, min_value=1e-4, step=0.01, key="mc_gpc")
    n_cyc = m2.number_input("사이클 수 N", value=200, min_value=1, step=10, key="mc_n")
    t_pulse_mc = m3.number_input("사이클당 펄스시간 [s]", value=5.0, min_value=0.0, step=0.5, key="mc_t")
    gpc_m = gpc_mc * 1e-9

    n_arr, w_arr, ear_arr, film_arr = multicycle.evolution(geom, L, w, gpc_m, int(n_cyc), z)
    if len(n_arr) == 0:
        st.error("입력 폭/GPC로는 첫 사이클 전에 이미 닫힘 상태입니다. 값을 확인하세요.")
    else:
        close_n = multicycle.cycles_to_close(w, gpc_m)
        cA, cB, cC = st.columns(3)
        cA.metric("최종 EAR", f"{ear_arr[-1]:.1f}", delta=f"{ear_arr[-1]-ear_arr[0]:+.1f}")
        cB.metric("최종 폭", f"{w_arr[-1]*1e9:.1f} nm",
                  delta=f"{(w_arr[-1]-w_arr[0])*1e9:+.1f} nm")
        cC.metric("클로깅까지", f"{close_n:,.0f} cycle")
        if close_n <= n_cyc:
            st.warning(f"⚠️ 약 {close_n:.0f} 사이클에서 피처가 닫힙니다(클로깅). "
                       f"그래프는 닫히기 전(N={int(n_arr[-1])})까지만 표시합니다.")

        mtab1, mtab2 = st.tabs(["EAR · 폭 vs 사이클", "침투깊이 vs 사이클"])
        with mtab1:
            figE = plots.fig_multicycle_ear(n_arr, ear_arr, w_arr * 1e9)
            st.pyplot(figE)
            cc1, cc2 = st.columns(2)
            cc1.download_button("⬇️ PNG", data=export.fig_png_bytes(figE),
                                file_name="multicycle_ear.png", mime="image/png", key="png_mc_ear")
            cc2.download_button("⬇️ 데이터 CSV",
                                data=export.series_csv(["cycle", "pore_width_nm", "EAR", "film_nm"],
                                                       [n_arr, w_arr*1e9, ear_arr, film_arr*1e9]),
                                file_name="multicycle_ear.csv", mime="text/csv", key="csv_mc_ear")
            plt.close(figE)
        with mtab2:
            if P > 0 and t_pulse_mc > 0:
                Pt_cycle = P * t_pulse_mc
                n_p, l_p = multicycle.penetration_per_cycle(geom, L, w, gpc_m, Pt_cycle,
                                                            K_max, m, T, int(n_cyc), z)
                if s_stick < 1.0:
                    _, l_p_lo = multicycle.penetration_per_cycle(
                        geom, L, w, gpc_m, Pt_cycle * s_stick, K_max, m, T, int(n_cyc), z)
                    figP = plots.fig_multicycle_penetration(n_p, l_p * 1e6, l_lo_um=l_p_lo * 1e6,
                                                            s_stick=s_stick)
                    csv_cols = (["cycle", "penetration_um_s1", "penetration_um_reaction_limited"],
                                [n_p, l_p * 1e6, l_p_lo * 1e6])
                else:
                    figP = plots.fig_multicycle_penetration(n_p, l_p * 1e6)
                    csv_cols = (["cycle", "penetration_um"], [n_p, l_p * 1e6])
                st.pyplot(figP)
                cc1, cc2 = st.columns(2)
                cc1.download_button("⬇️ PNG", data=export.fig_png_bytes(figP),
                                    file_name="multicycle_penetration.png", mime="image/png",
                                    key="png_mc_pen")
                cc2.download_button("⬇️ 데이터 CSV",
                                    data=export.series_csv(csv_cols[0], csv_cols[1]),
                                    file_name="multicycle_penetration.csv", mime="text/csv",
                                    key="csv_mc_pen")
                plt.close(figP)
                st.caption(f"per-cycle 노출 Pt = P·t = {P_in:g} {pres_unit} × {t_pulse_mc:g} s. "
                           "l ∝ w 이므로 폭이 줄며 침투 깊이 감소(Perez/Gordon)."
                           + ("" if s_stick >= 1.0 else f" 음영=반응제한(s={s_stick:g}) 범위."))
            else:
                st.warning("분압 P > 0 및 사이클당 펄스시간 > 0 을 입력하면 표시됩니다.")

# ===========================================================================
# 푸터
# ===========================================================================
with st.expander("사용한 입력값 (SI) 펼쳐보기 — 검증/디버깅용"):
    st.write({
        "T [K]": T, "m [kg]": m, "K_max [/m²]": K_max,
        "L [m]": L, "w [m]": w, "z [m]": z, "d [m]": d, "P [Pa]": P,
        "geometry": geom, "EAR (a)": a, "(Pt)_flat [Pa·s]": flat_pa_s,
    })

st.divider()
st.caption(
    f"Gordon ALD Conformality Calculator v{APP_VERSION} · 모델: Gordon et al., "
    "Chem. Vap. Depos. 9(2), 73–78 (2003) · 리뷰: Cremers et al., Appl. Phys. Rev. 6, 021302 (2019). "
    "결과는 모델 가정(s=1·분자흐름·확산제한)에 의존하며 절대값은 차수 수준으로 해석하십시오."
)
