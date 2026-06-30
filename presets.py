"""presets.py — reactant / 막(film) 프리셋 라이브러리.

방침(불확실성 표기):
  - precursor M, 막 M: 화학식 기반 → 신뢰 가능.
  - s(참고 sticking probability): 논문 Table V 인용. **계산 미사용**(Gordon은 s=1).
  - 막 밀도 ρ: 결정상/측정법 의존. ALD 막은 XRR 측정값 사용 권장(프리셋은 참고/추정).
출처: Cremers et al., Appl. Phys. Rev. 6, 021302 (2019), Table V;
      Mo: Vos et al. JVST A 41, 052402 (2023) 등(본문 검토 참조).
"""

# --- precursor 프리셋: M(몰질량) 및 참고 s ---
PRESETS = {
    "사용자 정의 (custom)":        {"M": 354.80, "s": "—",         "note": "직접 입력"},
    "Hf(NMe2)4 → HfO2":           {"M": 354.80, "s": "0.03–0.6",  "note": "Table V (Hf(NEtMe)4, MC)"},
    "TEMAHf Hf(NEtMe)4 → HfO2":   {"M": 410.91, "s": "0.03–0.6",  "note": "Table V (MC)"},
    "TMA Al(CH3)3 → Al2O3":       {"M": 72.09,  "s": "≈0.002–0.1","note": "Table V (방법 의존)"},
    "TiCl4 → TiN/TiO2":           {"M": 189.68, "s": "0.006–0.1", "note": "Table V (QCM~0.006, 연속체~0.1)"},
    "TDMAT Ti(NMe2)4 → TiO2":     {"M": 224.18, "s": "≈0.02",     "note": "Table V (MC)"},
    "DEZ ZnEt2 → ZnO":            {"M": 123.50, "s": "≈0.007",    "note": "Table V (MC)"},
    "Zr(NMe2)4 → ZrO2":           {"M": 267.53, "s": "≈0.07",     "note": "Table V (QCM)"},
    # --- Mo 염화물 (SJ 공정) ---
    "MoO2Cl2 (Mo/MoN, 염화물)":   {"M": 198.85, "s": "미공개",     "note": "thermal H2/NH3; HCl 부산물·고온 CVD 주의"},
    "MoCl5 (MoN, 염화물)":        {"M": 273.20, "s": "미공개",     "note": "이량화(Mo2Cl10) 가능; HCl 주의"},
}

# Mo 염화물 precursor (선택 시 Mo 경고 배너 표시)
MO_PRECURSORS = {"MoO2Cl2 (Mo/MoN, 염화물)", "MoCl5 (MoN, 염화물)"}

# --- 막(film) 프리셋: GPC→K_max 변환용 (M_film, ρ) ---
#   Mo계 밀도는 격자 기반 추정 — 실제는 XRR 측정값 사용 권장.
FILM_PRESETS = {
    "사용자 정의 (custom)":   {"M": 210.49, "rho": 9.68, "gpc": 0.10, "note": "직접 입력 (기본: HfO2)"},
    "Mo (metal, BCC)":        {"M": 95.95,  "rho": 10.22, "gpc": 0.022, "note": "벌크 Mo. 참고 GPC≈0.022 nm (Vos 2023)"},
    "γ-Mo2N (cubic)":         {"M": 205.91, "rho": 9.5,  "gpc": 0.03, "note": "격자 기반 추정 ρ — XRR 측정 권장"},
    "δ-MoN (hex, WC형)":      {"M": 109.96, "rho": 9.3,  "gpc": 0.03, "note": "격자 기반 추정 ρ — XRR 측정 권장"},
    "HfO2":                   {"M": 210.49, "rho": 9.68, "gpc": 0.10, "note": ""},
    "Al2O3 (ALD, 비정질)":    {"M": 101.96, "rho": 3.0,  "gpc": 0.11, "note": "비정질 ALD ρ≈3.0"},
    "TiO2 (anatase)":         {"M": 79.87,  "rho": 3.9,  "gpc": 0.05, "note": ""},
    "ZnO":                    {"M": 81.38,  "rho": 5.6,  "gpc": 0.20, "note": ""},
}
